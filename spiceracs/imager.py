#!/usr/bin/env python3
"""SPICE-RACS imager"""
import hashlib
import logging as log
import multiprocessing as mp
import os
import pickle
import shutil
from argparse import Namespace
from glob import glob
from typing import List
import traceback

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
from dask.distributed import Client, get_client, LocalCluster
from dask_mpi import initialize
from IPython import embed
from racs_tools import beamcon_2D
from radio_beam import Beam
from schwimmbad import SerialPool
from spython.main import Client as sclient
from tqdm.auto import tqdm

from spiceracs import fix_ms_dir
from spiceracs.utils import beam_from_ms, chunk_dask, inspect_client, wsclean, field_idx_from_ms


def get_wsclean(hub_address: str = "docker://alecthomson/wsclean", tag="3.1") -> str:
    """Pull wsclean image from dockerhub (or wherver).

    Args:
        version (str, optional): wsclean image tag. Defaults to "3.1".

    Returns:
        str: Path to yandasoft image.
    """
    sclient.load(f"{hub_address}:{tag}")
    image = os.path.abspath(sclient.pull())
    return image


def get_prefix(
    ms: str,
    out_dir: str,
) -> str:
    """Get prefix for output files"""
    idx = field_idx_from_ms(ms)
    field = vishead(vis=ms, mode="list")["field"][0][idx]
    beam = beam_from_ms(ms)
    prefix = f"image.{field}.contcube.beam{beam:02}"
    return os.path.join(out_dir, prefix)


@delayed
def image_beam(
    ms: str,
    field_idx: int,
    out_dir: str,
    prefix: str,
    simage: str,
    pols: str = "IQU",
    nchan: int = 36,
    scale: u.Quantity = 2.5 * u.arcsec,
    npix: int = 4096,
    join_polarizations: bool = True,
    join_channels: bool = True,
    squared_channel_joining: bool = True,
    mgain: float = 0.8,
    niter: int = 100_000,
    auto_mask: float = 3,
    auto_threshold: float = 1,
    use_wgridder: bool = True,
    robust: float = 0.0,
    mem: float = 90,
    taper: float = None,
):
    """Image a single beam"""
    import logging
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )
    log = logging.getLogger(__name__)
    # (
    #     addr,
    #     nworkers,
    #     nthreads,
    #     memory,
    #     threads_per_worker,
    #     memory_per_worker,
    # ) = inspect_client()
    # abs_mem = float(memory_per_worker.to(u.gigabyte).value * mem / 100)
    # log.debug(f"Using {abs_mem} GB of memory")

    command = wsclean(
        mslist=[ms],
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
        auto_threshold=auto_threshold,
        use_wgridder=use_wgridder,
        weight=f"briggs {robust}",
        log_time=False,
        # abs_mem=abs_mem,
        mem=mem,
        # no_mf_weighting=True,
        # j=1,
        taper_gaussian=f"{taper}asec" if taper else None,
        field=field_idx,
        parallel_deconvolution=2048,
    )

    root_dir = os.path.dirname(ms)

    # for command in (command_1, command_2):
    log.info(f"Running wsclean with command: {command}")
    output = sclient.execute(
        image=simage,
        command=command.split(),
        bind=f"{out_dir}:{out_dir}, {root_dir}:{root_dir}",
        return_result=True,
        quiet=False,
        stream=True,
    )
    for line in output:
        log.info(line)
    return True


def get_images(image_done: bool, pol, prefix):
    # image_lists = {s: sorted(glob(f"{prefix}*[0-9]-{s}-image.fits")) for s in pols}
    if image_done:
        image_list = sorted(glob(f"{prefix}*[0-9]-{pol}-image.fits"))
        return image_list
    else:
        raise ValueError("Imaging must be done")


@delayed
def get_aux(image_done: bool, pol, prefix):
    if image_done:
        aux_lists = {
            aux: sorted(
                glob(
                    f"{prefix}*[0-9]-{pol}-{aux}.fits" if aux != "psf" else f"{prefix}*[0-9]-{aux}.fits"
                )
            )
            for aux in ["model", "psf", "residual", "dirty"]
        }
        return aux_lists
    else:
        raise ValueError("Imaging must be done")


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
    log.info(f"Written {new_name}")

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


@delayed(nout=2)
def smooth_image(image, common_beam_pkl):
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
        )
    sm_image = image.replace(".fits", ".conv.fits")
    return sm_image


@delayed
def cleanup(
    im_cube_name,
    w_cube_name,
    image_list,
    aux_list,
    sm_image_list,
):
    for image in image_list:
        log.error(f"Removing {image}")
        os.remove(image)
    for sm_image in sm_image_list:
        log.error(f"Removing {sm_image}")
        os.remove(sm_image)
    for aux in aux_list.keys():
        for aux_image in aux_list[aux]:
            try:
                log.error(f"Removing {aux_image}")
                os.remove(aux_image)
            except FileNotFoundError:
                log.error(f"Could not find {aux_image}")
                log.error(f"aux_lists: {aux_list}")
    return


@delayed
def fix_ms(ms):
    fix_ms_dir.main(ms)
    return ms


def main(
    msdir: str,
    out_dir: str,
    cutoff: float = None,
    robust: float = -0.5,
    pols: str = "IQU",
    nchan: int = 36,
    # size: int = 4096,
    size: int = 6074,
    taper: float = None,
):
    simage = get_wsclean(tag="latest")
    msdir = os.path.abspath(msdir)
    out_dir = os.path.abspath(out_dir)
    get_image_task = delayed(get_images, nout=nchan)

    mslist = sorted(glob(os.path.join(msdir, "scienceData*_averaged_cal.leakage.ms")))

    assert (len(mslist) > 0) & (
        len(mslist) == 36
    ), f"Incorrect number of MS files found: {len(mslist)} / 36"

    mslist = [mslist[0]]

    log.info(f"Will image {len(mslist)} MS files in {msdir} to {out_dir}")
    cleans = []  # type: List[Delayed]

    # Do this in serial since CASA gets upset
    prefixs = {}
    field_idxs = {}
    for ms in tqdm(mslist, "Getting metadata"):
        prefix = get_prefix(ms, out_dir)
        prefixs[ms] = prefix
        field_idxs[ms] = field_idx_from_ms(ms)

    ms_dict = {}
    for ms in mslist:
        log.info(f"Imaging {ms}")
        # Apply Emil's fix for MSs feed centre
        ms_fix = fix_ms(ms)
        # Image with wsclean
        image_done = image_beam(
            ms=ms_fix,
            field_idx=field_idxs[ms],
            out_dir=out_dir,
            prefix=prefixs[ms],
            simage=simage,
            robust=robust,
            pols=pols,
            nchan=nchan,
            npix=size,
            taper=taper,
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
        type=str,
        help="Directory containing MS files",
    )
    parser.add_argument(
        "outdir",
        type=str,
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
        "--taper",
        type=float,
        default=None,
    )
    parser.add_argument(
        "--mpi",
        action="store_true",
        help="Use MPI",
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
        log.debug(f"{cluster=}")
        log.debug(f"{client=}")
        main(
            msdir=args.msdir,
            out_dir=args.outdir,
            cutoff=args.cutoff,
            robust=args.robust,
            nchan=args.nchan,
            pols=args.pols,
            size=args.size,
            taper=args.taper,
        )


if __name__ == "__main__":
    log.basicConfig(
        level=log.DEBUG,
        format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )
    cli()
