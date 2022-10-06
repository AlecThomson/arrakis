#!/usr/bin/env python3
"""SPICE-RACS imager"""
import os
import shutil
from glob import glob
import logging as log
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from spiceracs.utils import wsclean, beam_from_ms
from spiceracs import fix_ms_dir
import astropy.units as u
from casatasks import vishead
from spython.main import Client as sclient
from tqdm.auto import tqdm
from racs_tools import beamcon_2D
from argparse import Namespace
from IPython import embed
from dask.distributed import get_client, Client
from dask import delayed, compute

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
    field = vishead(ms, "list")["field"][0][0]
    beam = beam_from_ms(ms)
    prefix = f"image.{field}.contcube.beam{beam}"
    return os.path.join(out_dir, prefix)

@delayed(nout=3)
def image_beam(
    ms: str,
    out_dir: str,
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
    prefix = get_prefix(ms, out_dir)
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
        aux: sorted(glob(f"{prefix}*[0-9]-{aux}.fits")) for aux in ["model", "psf", "residual", "dirty"]
    }

    return image_lists, weight_list, aux_lists

@delayed(nout=2)
def make_cube(
    image_lists: dict,
    weight_list: list,
) -> tuple:
    """Make a cube from the images"""
    im_cube_names = []
    for s in tqdm("IQU",desc=f"Making Stokes {s} image cube"):
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
                hdul[0].data = fits.getdata(weight)

    return im_cube_names, weight_cube_names

@delayed
def smooth_images(image_lists):
    # convert dict to list
    out_lists = {}
    for s in image_lists.keys():
        # Smooth images
        beamcon_2D.main(
            pool=get_client(),
            infile=image_lists[s],
            suffix="conv",
            cutoff=25,
        )
        out_list  = []
        for image in image_lists[s]:
            out_list.append(image.replace(".fits", ".conv.fits"))
        out_lists[s] = out_list
    return out_lists

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
    for ms in mslist:
        log.info(f"Imaging {ms}")
        # Apply Emil's fix for MSs
        ms_fix = fix_ms(ms)
        # Image with wsclean
        image_lists, weight_list, aux_lists = image_beam(
            ms=ms_fix,
            out_dir=out_dir,
            image=image,
            niter=1,
        )
        # Smooth images
        sm_image_lists = smooth_images(image_lists)
        # Make a cube
        im_cube_name, w_cube_name = make_cube(
            image_lists=sm_image_lists,
            weight_list=weight_list,
        )
        # Clean up
        clean = cleanup(
            im_cube_name=im_cube_name,
            w_cube_name=w_cube_name,
            image_lists=image_lists,
            weight_list=weight_list,
            aux_lists=aux_lists,
        )
        cleans.append(clean)
        clean.compute()
        exit()
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

    args = parser.parse_args()
    main(
        msdir=args.msdir,
        out_dir=args.outdir,
    )

if __name__ == "__main__":
    log.basicConfig(
            level=log.INFO,
            format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            force=True,
    )
    client = Client()
    cli()