"""
The interpolation works as follows:
Take pixels offsets x,y from reference pixel in input image, multiply by
axis increments to get offx and offy.

Then compute offset = arcsin(offx^2+offy^2) and angle=atan2(offx,offy),
which should be the angular offset on the sky of the pixel position.

For the leakage image the inverse is used.
Take the offset and angle and turn them into pixel positions on the leakage map:

x = sin(offset)*cos(angle)/incx + refx
y = sin(offset)*sin(angle)/incy + refy
"""
import os
import warnings
from glob import glob

import astropy
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import mad_std, sigma_clip
from astropy.wcs import WCS
from dask import delayed
from dask.distributed import Client, LocalCluster

from arrakis.linmos import gen_seps
from arrakis.logger import logger, logging
from arrakis.utils.database import get_db
from arrakis.utils.fitsutils import getfreq
from arrakis.utils.pipeline import chunk_dask, logo_str


def make_plot(data, comp, imfile):
    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(10, 10))
    fig.suptitle(f"{comp['Gaussian_ID']} leakage")
    for i, s in enumerate(["q", "u"]):
        ax = axs[i]
        for beam, dat in data.items():
            freq = dat["freq"]
            frac = dat[f"{s}_image"] / dat["i_image"]
            filt = sigma_clip(frac, sigma=5, stdfunc=mad_std)
            frac[filt.mask] = np.nan
            (line,) = ax.step(
                freq,
                frac,
                label=f"beam {beam} -- off={dat['offset']:0.3f}, ang={dat['angle']:0.3f}",
            )
            ax.plot(freq, dat[f"{s}_holo"], ":", color=line.get_color())
        ax.set_ylabel(f"Stokes {s} [fractional]")

    plt.legend()
    plt.xlabel("Frequency [Hz]")

    outname = os.path.join("./", f"{comp['Gaussian_ID']}_leakage.pdf")
    plt.savefig(outname)
    plt.close()
    return outname


@delayed
def interpolate(field, comp, beams, cutdir, septab, holofile, verbose=True):
    beam = beams["beams"][field]

    ra = comp["RA"]
    dec = comp["Dec"]
    coord = SkyCoord(ra * u.deg, dec * u.deg)

    wcs_holo = WCS(holofile)
    incx, incy = astropy.wcs.utils.proj_plane_pixel_scales(wcs_holo.celestial)
    refx = int(wcs_holo.celestial.to_header()["CRPIX1"] - 1)
    refy = int(wcs_holo.celestial.to_header()["CRPIX2"] - 1)
    holo_data = fits.getdata(holofile)
    data = {}
    for bm in list(set(beam["beam_list"])):  # Ensure list of beams is unique!
        data.update({bm: {}})
        # imfile = beam[f'i_beam{bm}_image_file']
        try:
            imfile = glob(
                os.path.join(cutdir, f"{comp['Source_ID']}*beam{bm:02d}.conv.fits")
            )[0]
        except:
            logger.critical(f"No image file for source {comp['Source_ID']} beam {bm}")
            return

        freq = getfreq(imfile)
        wcs = WCS(imfile)
        sep = septab[bm]
        beamcoord = SkyCoord(sep["BEAM_RA"], sep["BEAM_DEC"], unit=(u.hourangle, u.deg))

        x, y = np.array(wcs.celestial.world_to_pixel(coord)).round().astype(int)

        offset = beamcoord.separation(coord).to(u.deg)
        angle = beamcoord.position_angle(coord).to(u.deg)

        data[bm].update({"offset": offset})
        data[bm].update({"angle": angle})

        x_holo = int(np.round(np.sin(offset) * np.sin(angle) / incx + refx))
        y_holo = int(np.round(np.sin(offset) * np.cos(angle) / incy + refy))
        # ell = offset * np.sin(angle)
        # emm = offset * np.cos(angle)
        # x_holo, y_holo = wcs_holo.celestial.world_to_array_index(SkyCoord(ell, emm))
        for i, s in enumerate(["i", "q", "u"]):
            # imfile = beam[f'{s}_beam{bm}_image_file']
            # imfile = os.path.join(os.path.abspath(cutdir), imfile)
            imfile = glob(f"{cutdir}/{comp['Source_ID']}*.{s}.*beam{bm:02d}.conv.fits")[
                0
            ]
            imdata = np.squeeze(fits.getdata(imfile))
            im_spec = imdata[:, y, x]
            filt = sigma_clip(im_spec, sigma=5, stdfunc=mad_std)
            im_spec[filt.mask] = np.nan
            holo_spec = holo_data[bm, i, :, y_holo, x_holo]

            data[bm].update({f"{s}_holo": holo_spec})
            data[bm].update({f"{s}_image": im_spec})
            data[bm].update({f"freq": freq})

    try:
        outname = make_plot(data, comp, imfile)
        # plotdir = os.path.join(os.path.join(cutdir, 'plots'), os.path.basename(outname))
        # copyfile(outname, plotdir)
    except Exception as e:
        logger.warning(f"No plot made : {e}")
        return


def main(
    field,
    datadir,
    host,
    epoch: int,
    holofile,
    username=None,
    password=None,
    verbose=True,
    snr_cut=None,
):
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    beamseps = gen_seps(field)

    if datadir is not None:
        datadir = os.path.abspath(datadir)

    cutdir = os.path.abspath(os.path.join(datadir, "cutouts"))
    holofile = os.path.abspath(holofile)

    beams_col, island_col, comp_col = get_db(
        host=host, epoch=epoch, username=username, password=password
    )

    # Query the DB
    beam_query = {"$and": [{f"beams.{field}": {"$exists": True}}]}
    island_ids = sorted(beams_col.distinct("Source_ID", beam_query))
    isl_query = {"Source_ID": {"$in": island_ids}}
    beams = pd.DataFrame(list(beams_col.find(isl_query).sort("Source_ID")))
    beams.set_index("Source_ID", drop=False, inplace=True)
    components = pd.DataFrame(
        list(
            comp_col.find(
                isl_query,
                # Only get required values
                {
                    "Source_ID": 1,
                    "Gaussian_ID": 1,
                    "RA": 1,
                    "Dec": 1,
                    "Noise": 1,
                    "Total_flux_Gaussian": 1,
                },
            ).sort("Source_ID")
        )
    )
    components.set_index("Source_ID", drop=False, inplace=True)
    component_ids = list(components["Gaussian_ID"])
    assert len(set(beams.index)) == len(set(components.index))

    outputs = []
    for i, c in components.iterrows():
        if snr_cut is not None:
            noise = c.Noise
            signal = c.Total_flux_Gaussian
            snr_total = signal / noise
            if snr_total < snr_cut:
                continue
        out = interpolate(
            field=field,
            comp=c,
            beams=beams.loc[c.Source_ID],
            cutdir=cutdir,
            septab=beamseps,
            holofile=holofile,
        )
        outputs.append(out)
    futures = chunk_dask(
        outputs=outputs,
        task_name="leakage plots",
        progress_text="Making leakage plots",
        verbose=verbose,
    )

    logger.info("Comparing leakge done!")


def cli():
    """Command-line interface"""
    import argparse

    from astropy.utils.exceptions import AstropyWarning

    warnings.simplefilter("ignore", category=AstropyWarning)
    from astropy.io.fits.verify import VerifyWarning

    warnings.simplefilter("ignore", category=VerifyWarning)
    warnings.simplefilter("ignore", category=RuntimeWarning)
    # Help string to be shown using the -h option
    descStr = f"""
    {logo_str}
    Arrakis:
    Make leakage comparison plots.

    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "field", metavar="field", type=str, help="RACS field - e.g. 2132-50A."
    )

    parser.add_argument(
        "datadir",
        metavar="datadir",
        type=str,
        help="Directory containing data cutout direcory (in datadir/cutouts).",
    )

    parser.add_argument("--holofile", type=str, help="Path to holography image")

    parser.add_argument(
        "--host",
        type=str,
        help="Host of mongodb.",
    )

    parser.add_argument(
        "--username", type=str, default=None, help="Username of mongodb."
    )

    parser.add_argument(
        "--password", type=str, default=None, help="Password of mongodb."
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose output [False]."
    )

    parser.add_argument("--snr", type=float, default=None, help="SNR cut (full band).")

    args = parser.parse_args()

    cluster = LocalCluster(
        n_workers=10,
        threads_per_worker=1,
    )
    client = Client(cluster)

    if args.verbose:
        logger.setLevel(logging.INFO)

    main(
        field=args.field,
        datadir=args.datadir,
        host=args.host,
        holofile=args.holofile,
        username=args.username,
        password=args.password,
        verbose=args.verbose,
        snr_cut=args.snr,
    )

    client.close()
    cluster.close()


if __name__ == "__main__":
    cli()
