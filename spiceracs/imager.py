#!/usr/bin/env python3
"""SPICE-RACS imager"""
import hashlib
import logging
import multiprocessing as mp
import os
import pickle
import shutil
import traceback
from argparse import Namespace
from glob import glob
from typing import List, Union, Dict
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import mad_std
from astropy.wcs import WCS
from casatasks import casalog, vishead
from dask import compute, delayed
from dask.delayed import Delayed
from dask.distributed import Client, LocalCluster, get_client
from dask_mpi import initialize
from IPython import embed
from racs_tools import beamcon_2D
from radio_beam import Beam
from schwimmbad import SerialPool
from spython.main import Client as sclient
from tqdm.auto import tqdm

from spiceracs import fix_ms_dir
from spiceracs.logger import logger
from spiceracs.utils import (
    beam_from_ms,
    chunk_dask,
    field_idx_from_ms,
    inspect_client,
    wsclean,
)


def get_wsclean(wsclean: Union[Path, str]) -> Path:
    """Pull wsclean image from dockerhub (or wherver).

    Args:
        version (str, optional): wsclean image tag. Defaults to "3.1".

    Returns:
        Path: Path to wsclean image.
    """
    sclient.load(str(wsclean))
    if isinstance(wsclean, str):
        return Path(sclient.pull(wsclean))
    return wsclean


def get_prefix(
    ms: Path,
    out_dir: Path,
) -> Path:
    """Get prefix for output files"""
    idx = field_idx_from_ms(ms.resolve(strict=True).as_posix())
    field = vishead(vis=ms.resolve(strict=True).as_posix(), mode="list")["field"][0][idx]
    beam = beam_from_ms(ms.resolve(strict=True).as_posix())
    prefix = f"image.{field}.contcube.beam{beam:02}"
    return out_dir / prefix


@delayed()
def image_beam(
    ms: Path,
    field_idx: int,
    out_dir: Path,
    prefix: str,
    simage: Path,
    pols: str = "IQU",
    nchan: int = 36,
    scale: u.Quantity = 2.5 * u.arcsec,
    npix: int = 4096,
    join_polarizations: bool = True,
    join_channels: bool = True,
    squared_channel_joining: bool = True,
    mgain: float = 0.7,
    niter: int = 100_000,
    auto_mask: float = 3,
    force_mask_rounds: Union[int, None] = None,
    auto_threshold: float = 1,
    gridder: Union[str, None] = None,
    robust: float = -0.5,
    mem: float = 90,
    taper: float = None,
    reimage: bool = False,
    minuv_l: float = 0.0,
    parallel_deconvolution: Union[int, None]=None,
    nmiter: Union[int, None]=None,
    local_rms: bool = False,
    local_rms_window: Union[float, None] = None,
):
    """Image a single beam"""
    logger = logging.getLogger(__name__)
    if not reimage:
        # Look for existing images
        checkvals = np.array(
            [f"{prefix}-{i:04}-{s}-image.fits" for s in pols for i in range(nchan)]
        )
        checks = np.array([os.path.exists(f) for f in checkvals])
        if all(checks):
            logger.critical("Images already exist, skipping imaging")
            return True
        else:
            logger.critical(f"At least some images do not exist: {checkvals[~checks]}")
            logger.critical("starting imaging...")

    commands = []
    # Do any I cleaning separately
    do_stokes_I = "I" in pols
    if do_stokes_I:
        command = wsclean(
            mslist=[ms.resolve(strict=True).as_posix()],
            use_mpi=False,
            name=prefix,
            pol="I",
            verbose=True,
            channels_out=nchan,
            scale=f"{scale.to(u.arcsec).value}asec",
            size=f"{npix} {npix}",
            join_polarizations=join_polarizations,
            join_channels=join_channels,
            squared_channel_joining=False,
            mgain=mgain,
            niter=niter,
            auto_mask=auto_mask,
            force_mask_rounds=force_mask_rounds,
            auto_threshold=auto_threshold,
            gridder=gridder,
            weight=f"briggs {robust}",
            log_time=False,
            mem=mem,
            taper_gaussian=f"{taper}asec" if taper else None,
            field=field_idx,
            parallel_deconvolution=parallel_deconvolution,
            minuv_l=minuv_l,
            nmiter=nmiter,
            local_rms=local_rms,
            local_rms_window=local_rms_window,
        )

    if squared_channel_joining:
        logger.info("Using squared channel joining")
        logger.info("Reducing mask and threshold by sqrt(2) to account for this")
        auto_mask /= np.sqrt(2)
        auto_threshold /= np.sqrt(2)
        logger.debug(f"auto_mask = {auto_mask}")
        logger.debug(f"auto_threshold = {auto_threshold}")


    command = wsclean(
        mslist=[ms.resolve(strict=True).as_posix()],
        use_mpi=False,
        name=prefix,
        pol=pols,
        verbose=True,
        channels_out=nchan,
        scale=f"{scale.to(u.arcsec).value}asec",
        size=f"{npix} {npix}",
        join_polarizations=join_polarizations,
        join_channels=join_channels,
        squared_channel_joining=squared_channel_joining,
        mgain=mgain,
        niter=niter,
        auto_mask=auto_mask,
        force_mask_rounds=force_mask_rounds,
        auto_threshold=auto_threshold,
        gridder=gridder,
        weight=f"briggs {robust}",
        log_time=False,
        mem=mem,
        taper_gaussian=f"{taper}asec" if taper else None,
        field=field_idx,
        parallel_deconvolution=parallel_deconvolution,
        minuv_l=minuv_l,
        nmiter=nmiter,
        local_rms=local_rms,
        local_rms_window=local_rms_window,
    )
    commands.append(command)

    root_dir = ms.parent

    for command in commands:
        logger.info(f"Running wsclean with command: {command}")
        output = sclient.execute(
            image=simage.resolve(strict=True).as_posix(),
            command=command.split(),
            bind=f"{out_dir}:{out_dir}, {root_dir.resolve(strict=True).as_posix()}:{root_dir.resolve(strict=True).as_posix()}",
            return_result=True,
            quiet=False,
            stream=True,
        )
        for line in output:
            logger.info(line)

    # Check rms of image to check for divergence
    if do_stokes_I:
        pols += "I"
    for pol in pols:
        mfs_image = f"{prefix}-MFS-image.fits" if pol == "I" else f"{prefix}-MFS-{pol}-image.fits"
        rms = mad_std(fits.getdata(mfs_image), ignore_nan=True)
        if rms > 1:
            raise ValueError(f"RMS of {rms} is too high in image {mfs_image}, try imaging with lower mgain {mgain - 0.1}")
    return True


def get_images(image_done: bool, pol: str, prefix: Path) -> List[Path]:
    # image_lists = {s: sorted(glob(f"{prefix}*[0-9]-{s}-image.fits")) for s in pols}
    if not image_done:
        raise ValueError("Imaging must be done")

    imglob = "*[0-9]-image.fits" if pol=="I" else f"*[0-9]-{pol}-image.fits"
    image_list = sorted(prefix.glob(imglob))
    return image_list



@delayed()
def get_aux(image_done: bool, pol: str, prefix: Path) -> Dict[str, List[Path]]:
    if not image_done:
        raise ValueError("Imaging must be done")
    aux_lists = {}
    for aux in ["model", "psf", "residual", "dirty"]:
        if pol == "I" or aux == "psf":
            aux_lists[aux] = sorted(prefix.glob(f"*[0-9]-{aux}.fits"))
        else:
            aux_lists[aux] = sorted(prefix.glob(f"*[0-9]-{pol}-{aux}.fits"))
    return aux_lists



@delayed(nout=2)
def make_cube(
    pol: str,
    image_list: list,
    common_beam_pkl: str,
) -> tuple:
    """Make a cube from the images"""
    # First combine images into cubes
    freqs = []
    rmss = []
    for chan, image in enumerate(
        tqdm(
            image_list,
            desc="Reading channel image",
            leave=False,
        )
    ):
        # init cube
        if chan == 0:
            old_name = image
            old_header = fits.getheader(old_name)
            wcs = WCS(old_header)
            idx = 0
            for j, t in enumerate(
                wcs.axis_type_names[::-1]
            ):  # Reverse to match index order
                if t == "FREQ":
                    idx = j
                    break

            plane_shape = list(fits.getdata(old_name).shape)
            cube_shape = plane_shape.copy()
            cube_shape[idx] = len(image_list)

            data_cube = np.zeros(cube_shape)

            out_dir = os.path.dirname(old_name)
            old_base = os.path.basename(old_name)
            new_base = old_base
            b_idx = new_base.find("beam") + len("beam") + 2
            sub = new_base[b_idx:]
            new_base = new_base.replace(sub, ".fits")
            new_base = new_base.replace("image", f"image.restored.{pol.lower()}")
            new_name = os.path.join(out_dir, new_base)

        plane = fits.getdata(image) / 2  # Divide by 2 because of ASKAP Stokes
        plane_rms = mad_std(plane, ignore_nan=True)
        rmss.append(plane_rms)
        data_cube[:, chan] = plane
        freq = WCS(image).spectral.pixel_to_world(0)
        freqs.append(freq.to(u.Hz).value)
    # Write out cubes
    freqs = np.array(freqs) * u.Hz
    rmss_arr = np.array(rmss) * u.Jy / u.beam
    assert np.diff(freqs).std() < 1e-6 * u.Hz, "Frequencies are not evenly spaced"
    new_header = old_header.copy()
    new_header["NAXIS"] = len(cube_shape)
    new_header["NAXIS3"] = len(freqs)
    new_header["CRPIX3"] = 1
    new_header["CRVAL3"] = freqs[0].value
    new_header["CDELT3"] = np.diff(freqs).mean().value
    new_header["CUNIT3"] = "Hz"
    # Deserialise beam
    with open(common_beam_pkl, "rb") as f:
        common_beam = pickle.load(f)
    new_header = common_beam.attach_to_header(new_header)
    fits.writeto(new_name, data_cube, new_header, overwrite=True)
    logger.info(f"Written {new_name}")

    # Copy image cube
    new_w_name = new_name.replace("image.restored", "weights").replace(".fits", ".txt")
    np.savetxt(new_w_name, rmss_arr.value, fmt="%s")

    return new_name, new_w_name


@delayed(nout=2)
def get_beam(ms_dict, pols, cutoff=None):
    # convert dict to list
    image_list = []
    for ms in ms_dict.keys():
        for s in pols:
            image_lists = ms_dict[ms]["image_lists"]
            image_list.extend(image_lists[s])
    # Create a unique hash for the beam log filename
    beam_log = f"beam_{hashlib.md5(''.join(image_list).encode()).hexdigest()}.log"
    with SerialPool() as pool:
        common_beam = beamcon_2D.main(
            pool=pool,
            infile=image_list,
            cutoff=cutoff,
            dryrun=True,
            log=os.path.join("/tmp", beam_log),
        )
    # serialise the beam
    common_beam_pkl = os.path.abspath(
        f"beam_{hashlib.md5(''.join(image_list).encode()).hexdigest()}.pkl"
    )

    with open(common_beam_pkl, "wb") as f:
        pickle.dump(common_beam, f)

    return common_beam_pkl, beam_log


@delayed()
def smooth_image(image, common_beam_pkl, cutoff=None):
    # Smooth image
    # Deserialise the beam
    with open(common_beam_pkl, "rb") as f:
        common_beam = pickle.load(f)
    with SerialPool() as pool:
        _ = beamcon_2D.main(
            pool=pool,
            infile=[image],
            suffix="conv",
            bmaj=common_beam.major.to(u.arcsec).value,
            bmin=common_beam.minor.to(u.arcsec).value,
            bpa=common_beam.pa.to(u.deg).value,
            cutoff=cutoff,
        )
    sm_image = image.replace(".fits", ".conv.fits")
    return sm_image


@delayed()
def cleanup(
    purge: bool,
    im_cube_name,
    w_cube_name,
    image_list,
    aux_list,
    sm_image_list,
):
    if not purge:
        logger.info("Not purging intermediate files")
        return
    for image in image_list:
        logger.critical(f"Removing {image}")
        os.remove(image)
    for sm_image in sm_image_list:
        logger.critical(f"Removing {sm_image}")
        os.remove(sm_image)
    for aux in aux_list.keys():
        for aux_image in aux_list[aux]:
            try:
                logger.critical(f"Removing {aux_image}")
                os.remove(aux_image)
            except FileNotFoundError:
                logger.error(f"Could not find {aux_image}")
                logger.error(f"aux_lists: {aux_list}")
    return


@delayed()
def fix_ms(ms: Path) -> Path:
    fix_ms_dir.main(ms.resolve(strict=True).as_posix())
    return ms


def main(
    msdir: Path,
    out_dir: Path,
    cutoff: Union[float, None] = None,
    robust: float = -0.5,
    pols: str = "IQU",
    nchan: int = 36,
    size: int = 6074,
    scale: u.Quantity = 2.5 * u.arcsec,
    mgain: float = 0.8,
    niter: int = 100_000,
    auto_mask: float = 3,
    force_mask_rounds: Union[int, None] = None,
    auto_threshold: float = 1,
    taper: Union[float, None] = None,
    reimage: bool = False,
    purge: bool = False,
    minuv: float = 0.0,
    parallel_deconvolution: Union[int, None] = None,
    gridder: Union[str, None] = None,
    nmiter: Union[int, None] = None,
    local_rms: bool = False,
    local_rms_window: Union[float, None] = None,
    wsclean_path: Union[Path, str] = "docker://alecthomson/wsclean:latest"
):
    simage = get_wsclean(wsclean=wsclean_path)
    get_image_task = delayed(get_images, nout=nchan)

    mslist = sorted(msdir.glob("scienceData*_averaged_cal.leakage.ms"))

    assert (len(mslist) > 0) & (
        len(mslist) == 36
    ), f"Incorrect number of MS files found: {len(mslist)} / 36"

    logger.info(f"Will image {len(mslist)} MS files in {msdir} to {out_dir}")
    cleans = []  # type: List[Delayed]

    # Do this in serial since CASA gets upset
    prefixs = {}
    field_idxs = {}
    for ms in tqdm(mslist, "Getting metadata"):
        prefix = get_prefix(ms, out_dir)
        prefixs[ms] = prefix
        field_idxs[ms] = field_idx_from_ms(ms.resolve(strict=True).as_posix())

    ms_dict = {}
    for ms in mslist:
        logger.info(f"Imaging {ms}")
        # Apply Emil's fix for MSs feed centre
        ms_fix = fix_ms(ms)
        # Image with wsclean
        image_done = image_beam(
            ms=ms_fix,
            field_idx=field_idxs[ms],
            out_dir=out_dir,
            prefix=prefixs[ms].resolve(strict=False).as_posix(),
            simage=simage.resolve(strict=True),
            robust=robust,
            pols=pols,
            nchan=nchan,
            scale=scale,
            npix=size,
            mgain=mgain,
            niter=niter,
            auto_mask=auto_mask,
            force_mask_rounds=force_mask_rounds,
            auto_threshold=auto_threshold,
            taper=taper,
            reimage=reimage,
            minuv_l=minuv,
            parallel_deconvolution=parallel_deconvolution,
            gridder=gridder,
            nmiter=nmiter,
            local_rms=local_rms,
            local_rms_window=local_rms_window,
        )
        # Get images
        image_lists = {}
        aux_lists = {}
        for s in pols:
            image_list = get_image_task(
                image_done=image_done,
                pol=s,
                prefix=prefixs[ms],
            )
            image_lists[s] = image_list
            aux_list = get_aux(
                image_done=image_done,
                pol=s,
                prefix=prefixs[ms],
            )
            aux_lists[s] = aux_list
        ms_dict[ms] = {
            "image_lists": image_lists,
            "aux_lists": aux_lists,
        }
        # Smooth images
    common_beam_pkl, beam_log = get_beam(
        ms_dict=ms_dict,
        pols=pols,
        cutoff=cutoff,
    )

    for ms in mslist:
        sm_image_lists = {}
        image_lists = ms_dict[ms]["image_lists"]
        aux_lists = ms_dict[ms]["aux_lists"]

        for s in pols:
            sm_image_list = []
            image_list = image_lists[s]
            aux_list = aux_lists[s]
            for image in image_list:
                sm_image = smooth_image(
                    image,
                    common_beam_pkl=common_beam_pkl,
                    cutoff=cutoff,
                )
                sm_image_list.append(sm_image)

            # Make a cube
            im_cube_name, w_cube_name = make_cube(
                pol=s,
                image_list=sm_image_list,
                common_beam_pkl=common_beam_pkl,
            )
            sm_image_lists[s] = sm_image_list
            # Clean up
            clean = cleanup(
                purge=purge,
                im_cube_name=im_cube_name,
                w_cube_name=w_cube_name,
                image_list=image_list,
                sm_image_list=sm_image_list,
                aux_list=aux_list,
            )
            cleans.append(clean)

    futures = chunk_dask(
        outputs=cleans,
        task_name="Image and cube",
        progress_text="Imaging",
        verbose=True,
        batch_size=1,
    )

    return [f.result() for f in futures]


def cli():
    import argparse

    """Command-line interface"""
    # Help string to be shown using the -h option
    logostr = """
     mmm   mmm   mmm   mmm   mmm
     )-(   )-(   )-(   )-(   )-(
    ( S ) ( P ) ( I ) ( C ) ( E )
    |   | |   | |   | |   | |   |
    |___| |___| |___| |___| |___|
     mmm     mmm     mmm     mmm
     )-(     )-(     )-(     )-(
    ( R )   ( A )   ( C )   ( S )
    |   |   |   |   |   |   |   |
    |___|   |___|   |___|   |___|

    """

    # Help string to be shown using the -h option
    descStr = f"""
    {logostr}
    SPICE-RACS Stage X:
    Image calibrated visibilities

    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "msdir",
        type=Path,
        help="Directory containing MS files",
    )
    parser.add_argument(
        "outdir",
        type=Path,
        help="Directory to output images",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        help="Cutoff for smoothing",
    )
    parser.add_argument(
        "--robust",
        type=float,
        default=-0.5,
    )
    parser.add_argument(
        "--nchan",
        type=int,
        default=36,
    )
    parser.add_argument(
        "--pols",
        type=str,
        default="IQU",
    )
    parser.add_argument(
        "--size",
        type=int,
        default=4096,
    )
    parser.add_argument(
        "--scale",
        type=u.Quantity,
        default=2.5,
    )
    parser.add_argument(
        "--mgain",
        type=float,
        default=0.8,
    )
    parser.add_argument(
        "--niter",
        type=int,
        default=100_000,
    )
    parser.add_argument(
        "--nmiter",
        type=int,
        default=None,
    )
    parser.add_argument(
        "--auto_mask",
        type=float,
        default=3.0,
    )
    parser.add_argument(
        "--auto_threshold",
        type=float,
        default=1.0,
    )
    parser.add_argument(
        "--local-rms",
        action="store_true",
    )
    parser.add_argument(
        "--local-rms-window",
        type=float,
        default=None,
    )
    parser.add_argument(
        "--force-mask-rounds",
        type=int,
        default=None,
    )
    parser.add_argument(
        "--gridder",
        type=str,
        default=None,
        choices=["direct-ft", "idg", "wgridder", "tuned-wgridder", "wstacking"],
    )
    parser.add_argument(
        "--taper",
        type=float,
        default=None,
    )
    parser.add_argument(
        "--minuv",
        type=float,
        default=0.0,
    )
    parser.add_argument(
        "--parallel",
        type=int,
        default=None,
    )
    parser.add_argument(
        "--purge",
        action="store_true",
        help="Purge intermediate files",
    )
    parser.add_argument(
        "--mpi",
        action="store_true",
        help="Use MPI",
    )
    parser.add_argument(
        "--reimage",
        action="store_true",
        help="Force a new round of imaging. Otherwise, will skip if images already exist.",
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--hosted-wsclean",
        type=str,
        default="docker://alecthomson/wsclean:latest",
        help="Docker or Singularity image for wsclean [docker://alecthomson/wsclean:latest]",
    )
    group.add_argument(
        "--local-wsclean",
        type=Path,
        default=None,
        help="Path to local wsclean Singularity image",
    )


    args = parser.parse_args()

    if args.mpi:
        initialize(interface="ipogif0")
        cluster = None

    else:
        cluster = LocalCluster(
            n_workers=1,
            threads_per_worker=1,
            # processes=False,
        )

    with Client(cluster) as client:
        logger.debug(f"{cluster=}")
        logger.debug(f"{client=}")
        main(
            msdir=args.msdir,
            out_dir=args.outdir,
            cutoff=args.cutoff,
            robust=args.robust,
            pols=args.pols,
            nchan=args.nchan,
            size=args.size,
            scale=args.scale,
            mgain=args.mgain,
            niter=args.niter,
            auto_mask=args.auto_mask,
            force_mask_rounds=args.force_mask_rounds,
            auto_threshold=args.auto_threshold,
            minuv=args.minuv,
            purge=args.purge,
            taper=args.taper,
            reimage=args.reimage,
            parallel_deconvolution=args.parallel,
            gridder=args.gridder,
            wsclean_path=Path(args.local_wsclean) if args.local_wsclean else args.hosted_wsclean,
        )


if __name__ == "__main__":
    cli()
