#!/usr/bin/env python3
"""SPICE-RACS imager"""
import os
from glob import glob
from spiceracs.utils import wsclean
import astropy.units as u
from casatasks import vishead
import subprocess as sp

def get_prefix(
    ms: str,
    out_dir: str,
) -> str:
    """Get prefix for output files"""
    field = vishead(ms, "list")["field"][0][0]
    beam = ms[ms.find("beam")+len("beam"):ms.find("beam")+len("beam")+2]
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

    return sorted(glob(f"{prefix}*.fits"))

def make_cube(
    msdir: str,
    out_dir: str,
):
    mslist = glob(
        os.path.join(msdir, "scienceData*_averaged_cal.leakage.ms")
    )
    assert len(mslist) > 0, "No MS files found"
    assert len(mslist) == 36, f"Incorrect number of MS files found: {len(mslist)}"

    for ms in mslist:
        fixed_ms = fix_ms(ms)
        image_list = image_beam(
            ms=fixed_ms,
            out_dir=out_dir,
        )

def main(
    mslist,
):
    pass

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
    SPICE-RACS Stage 5:
    Run RMsynthesis on cubelets.

    Note: Runs on brightest sources first.

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