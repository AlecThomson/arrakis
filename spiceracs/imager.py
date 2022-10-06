#!/usr/bin/env python3
"""SPICE-RACS imager"""
import os
from glob import glob
from spiceracs.utils import wsclean, beam_from_ms
from spiceracs import fix_ms_dir.main as fix_ms
import astropy.units as u
from casatasks import vishead
import subprocess as sp

def get_prefix(
    ms: str,
    out_dir: str,
) -> str:
    """Get prefix for output files"""
    field = vishead(ms, "list")["field"][0][0]
    beam = beam_from_ms(ms)
    prefix = f"image.{field}.contcube.beam{beam}"
    return os.path.join(out_dir, prefix)

def image_beam(
    ms: str,
    out_dir: str,
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
        log_time=True,
        temp_dir=out_dir,
        mem=mem,
        save_weights=True,
    )

    sp.run(command.split(), check=True)
    image_list = sorted(glob(f"{prefix}*.fits"))
    weight_list = sorted(glob(f"{prefix}*.weight.fits"))

    return image_list, weight_list

def make_cube(
    image_list: list,
    weight_list: list,
):
    """Make a cube from the images"""
    im_cube_name = image_list[0].replace(".fits", ".cube.fits")
    w_cube_name = weight_list[0].replace(".fits", ".cube.fits")


    return im_cube_name, w_cube_name


def main(
    msdir: str,
    out_dir: str,
):
    mslist = glob(
        os.path.join(msdir, "scienceData*_averaged_cal.leakage.ms")
    )
    assert len(mslist) > 0, "No MS files found"
    assert len(mslist) == 36, f"Incorrect number of MS files found: {len(mslist)}"

    for ms in mslist:
        # Apply Emil's fix for MSs
        fixed_ms = fix_ms(ms)
        # Image with wsclean
        image_list, weight_list = image_beam(
            ms=fixed_ms,
            out_dir=out_dir,
        )
        # Make a cube
        im_cube_name, w_cube_name = make_cube(
            image_list=image_list,
            weight_list=weight_list,
        )
        # Clean up
        for image in image_list:
            os.remove(image)
        for weight in weight_list:
            os.remove(weight)



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

    main(
        msdir,
    )

if __name__ == "__main__":
    cli()