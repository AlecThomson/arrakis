#!/usr/bin/env python3
"""SPICE-RACS imager"""
import os
import shutil
import pickle
from glob import glob
import logging as log
from bleach import clean
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from spiceracs.utils import wsclean, beam_from_ms
from casatasks import vishead, casalog
from spiceracs import fix_ms_dir
import astropy.units as u
from spython.main import Client as sclient
from tqdm.auto import tqdm
from racs_tools import beamcon_2D
from schwimmbad import SerialPool
import multiprocessing as mp
from argparse import Namespace
from IPython import embed
from dask.distributed import get_client, Client
from dask import delayed, compute
import hashlib
from radio_beam import Beam

def get_wsclean(
    hub_address: str = "docker://alecthomson/wsclean",
    tag="3.1"
) -> str:
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
    field = vishead(vis=ms, mode="list")["field"][0][0]
    beam = beam_from_ms(ms)
    prefix = f"image.{field}.contcube.beam{beam:02}"
    return os.path.join(out_dir, prefix)

@delayed(nout=3)
def image_beam(
    ms: str,
    out_dir: str,
    prefix: str,
    image: str,
    pols: str = "IQU",
    nchan: int = 36,
    scale: u.Quantity = 2.5 * u.arcsec,
    npix: int = 4096,
    join_polarizations: bool = True,
    join_channels: bool = True,
    squared_channel_joining: bool = True,
    mgain: float = 0.8,
    niter: int = 50000,
    auto_mask: float = 5,
    auto_threshold: float = 3,
    use_wgridder: bool = True,
    robust: float = 0.0,
    mem: float = 90,
):
    """Image a single beam

    """
    client = get_client()
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
        temp_dir=out_dir,
        mem=mem,
        save_weights=True,
        no_mf_weighting=True,
        j=[i for i in client.nthreads().values()][0], # Set number of threads to match dask
        # j=1,
    )

    root_dir = os.path.dirname(ms)

    log.info(f"Running wsclean with command: {command}")

    output = sclient.execute(
        image=image,
        command=command.split(),
        bind=f"{out_dir}:{out_dir}, {root_dir}:{root_dir}",
        return_result=True,
        quiet=False,
        stream=True
    )
    for line in output:
        log.info(line)
    image_lists = {s: sorted(glob(f"{prefix}*[0-9]-{s}-image.fits")) for s in pols}
    weight_list = sorted(glob(f"{prefix}*[0-9]-weights.fits"))
    aux_lists = {
        aux: sorted(glob(f"{prefix}*[0-9]-{s}-{aux}.fits")) for aux in ["model", "psf", "residual", "dirty" ] for s in pols
    }

    return image_lists, weight_list, aux_lists

@delayed(nout=2)
def make_cube(
    image_lists: dict,
    weight_list: list,
    common_beam_pkl: str,
) -> tuple:
    """Make a cube from the images"""
    im_cube_names = []
    for s in tqdm("IQU", desc=f"Making Stokes image cube"):
        # First combine images into cubes
        freqs = []
        for chan, image in enumerate(
                tqdm(
                    image_lists[s],
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
                for j,t in enumerate(wcs.axis_type_names[::-1]): # Reverse to match index order
                    if t == "FREQ":
                        idx = j
                        break

                plane_shape = list(fits.getdata(old_name).shape)
                cube_shape = plane_shape.copy()
                cube_shape[idx] = len(image_lists[s])

                data_cube = np.zeros(cube_shape)

                out_dir = os.path.dirname(old_name)
                old_base = os.path.basename(old_name)
                new_base = old_base
                b_idx = new_base.find("beam") + len("beam") + 2
                sub = new_base[b_idx:]
                new_base = new_base.replace(sub, ".fits")
                new_base = new_base.replace("image", f"image.restored.{s.lower()}")
                new_name = os.path.join(out_dir, new_base)

            data_cube[:,chan] = fits.getdata(image)
            freq = WCS(image).spectral.pixel_to_world(0)
            freqs.append(freq.to(u.Hz).value)
        # Write out cubes
        freqs = np.array(freqs) * u.Hz
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
        im_cube_names.append(new_name)

    # Now combine weights into cubes
    weight_cube_names = []
    for i, s in enumerate(tqdm("IQU",desc=f"Making Stokes {s} weights cube")):
        # Copy image cube
        old_name = im_cube_names[i]
        new_name = old_name.replace("image.restored", "weights")
        if os.path.exists(new_name):
            os.remove(new_name)
        shutil.copy(old_name, new_name)
        weight_cube_names.append(new_name)
        with fits.open(new_name, mode="update") as hdul:
            for chan, weight in enumerate(
                    tqdm(
                        weight_list,
                        desc="Reading channel weights",
                        leave=False,
                    )
                ):
                hdul[0].data[:,chan] = fits.getdata(weight)

    return im_cube_names, weight_cube_names

@delayed(nout=2)
def get_beam(image_lists, cutoff=None):
    # convert dict to list
    image_list = []
    for s in image_lists.keys():
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
    common_beam_pkl = os.path.join(
        "/tmp",
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
def smooth_images(image_lists, common_beam_pkl):
    sm_image_lists = {}
    for s in image_lists.keys():
        sm_image_lists[s] = []
        for image in image_lists[s]:
            sm_image = smooth_image(
                image,
                common_beam_pkl=common_beam_pkl,
            )
            sm_image_lists[s].append(sm_image)
    # return delayed(sm_image_lists)
    return sm_image_lists

@delayed
def cleanup(
    im_cube_name,
    w_cube_name,
    image_lists,
    weight_list,
    aux_lists,
    sm_image_lists,
):
    for s in image_lists.keys():
        for image in image_lists[s]:
            os.remove(image)
        for image in sm_image_lists[s]:
            os.remove(image)
    for weight in weight_list:
        os.remove(weight)
    for aux in aux_lists.keys():
        for image in aux_lists[aux]:
            os.remove(image)
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
):
    image = get_wsclean(tag="latest")
    msdir = os.path.abspath(msdir)
    out_dir = os.path.abspath(out_dir)

    mslist = glob(
        os.path.join(msdir, "scienceData*_averaged_cal.leakage.ms")
    )
    assert (len(mslist) > 0) & (len(mslist) == 36), f"Incorrect number of MS files found: {len(mslist)}"

    log.info(
        f"Will image {len(mslist)} MS files in {msdir} to {out_dir}"
    )
    cleans = []

    # Do this in serial since
    prefixs = {}
    for ms in mslist:
        prefix = get_prefix(ms, out_dir)
        prefixs[ms] = prefix

    for ms in mslist:
        log.info(f"Imaging {ms}")
        # Apply Emil's fix for MSs
        ms_fix = fix_ms(ms)
        # Image with wsclean
        image_lists, weight_list, aux_lists = image_beam(
            ms=ms_fix,
            out_dir=out_dir,
            prefix=prefixs[ms],
            image=image,
            niter=1,
            robust=robust,
        )
        # Smooth images
        common_beam_pkl, beam_log = get_beam(image_lists, cutoff=cutoff)
        sm_image_lists = smooth_images(image_lists, common_beam_pkl=common_beam_pkl)

        # Make a cube
        im_cube_name, w_cube_name = make_cube(
            image_lists=sm_image_lists,
            weight_list=weight_list,
            common_beam_pkl=common_beam_pkl,
        )
        # Clean up
        clean = cleanup(
            im_cube_name=im_cube_name,
            w_cube_name=w_cube_name,
            image_lists=image_lists,
            sm_image_lists=sm_image_lists,
            weight_list=weight_list,
            aux_lists=aux_lists,
        )
        cleans.append(clean)

    compute(*cleans)


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

    args = parser.parse_args()
    main(
        msdir=args.msdir,
        out_dir=args.outdir,
        cutoff=args.cutoff,
        robust=args.robust,
    )

if __name__ == "__main__":
    log.basicConfig(
            level=log.DEBUG,
            format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            force=True,
    )
    client = Client(n_workers=4)
    cli()