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
from IPython.core.pylabtools import figsize
import numpy as np
import warnings
from astropy.coordinates import SkyCoord
import astropy.units as u
import astropy
from astropy.io import fits
from astropy.wcs import WCS
from spiceracs.linmos import gen_seps
from spiceracs.utils import tqdm_dask, get_db, test_db, coord_to_string, getfreq, chunk_dask
from dask import delayed
from dask.distributed import LocalCluster, Client
import matplotlib.pyplot as plt
from shutil import copyfile
import time
from astropy.stats import sigma_clip, mad_std

def make_plot(data, comp, imfile):

    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(10,10))
    fig.suptitle(f"{comp['Gaussian_ID']} leakage")
    for i, s in enumerate(['q', 'u']):
        ax = axs[i]
        for beam, dat in data.items():
            freq = dat['freq']
            frac = dat[f'{s}_image'] / dat['i_image']
            filt = sigma_clip(frac, sigma=5, stdfunc=mad_std)
            frac[filt.mask] = np.nan
            line, = ax.step(freq, frac, label = f"beam {beam} -- off={dat['offset']:0.3f}, ang={dat['angle']:0.3f}")
            ax.plot(freq, dat[f'{s}_holo'], ':', color = line.get_color())
        ax.set_ylabel(f'Stokes {s} [fractional]')
        
    plt.legend()
    plt.xlabel('Frequency [Hz]')
    
    outname = os.path.join(os.path.dirname(imfile),f"{comp['Gaussian_ID']}_leakage.pdf")
    plt.savefig(outname)
    plt.close()
    return outname

@delayed
def interpolate(
    field,
    comp,
    beams,
    cutdir,
    septab,
    holofile,
    verbose=True
):
    beam = beams["beams"][field]

    ra = comp['RA']
    dec = comp['Dec']
    coord = SkyCoord(ra*u.deg, dec*u.deg)

    wcs_holo = WCS(holofile)
    incx, incy = astropy.wcs.utils.proj_plane_pixel_scales(wcs_holo.celestial)
    refx = int(wcs_holo.celestial.to_header()['CRPIX1'] - 1)
    refy = int(wcs_holo.celestial.to_header()['CRPIX2'] - 1)
    holo_data = fits.getdata(holofile)

    data = {}

    for bm in list(set(beam['beam_list'])):  # Ensure list of beams is unique!
        data.update({bm: {}})
        imfile = beam[f'i_beam{bm}_image_file']
        imfile = os.path.join(os.path.abspath(cutdir), imfile)

        freq = getfreq(imfile)
        wcs = WCS(imfile)
        sep = septab[bm]
        beamcoord = SkyCoord(sep['BEAM_RA'], sep['BEAM_DEC'], unit=(u.hourangle, u.deg))
        
        x, y = np.array(
            wcs.celestial.world_to_pixel(coord)
        ).round().astype(int)

        offset = coord.separation(beamcoord).to(u.deg)
        angle = coord.position_angle(beamcoord).to(u.deg)

        data[bm].update({"offset": offset})
        data[bm].update({"angle": angle})
        
        x_holo = int(np.round(np.sin(offset) * np.sin(angle) / incx + refx))
        y_holo = int(np.round(np.sin(offset) * np.cos(angle) / incy + refy))

        for i, s in enumerate(['i','q','u']):
            imfile = beam[f'{s}_beam{bm}_image_file']
            imfile = os.path.join(os.path.abspath(cutdir), imfile)
            imdata = np.squeeze(fits.getdata(imfile))
            im_spec = imdata[:, y, x]
            filt = sigma_clip(im_spec, sigma=5, stdfunc=mad_std)
            im_spec[filt.mask] = np.nan
            holo_spec = holo_data[bm, i, :, y_holo, x_holo]

            data[bm].update({f"{s}_holo": holo_spec})
            data[bm].update({f"{s}_image": im_spec})
            data[bm].update({f"freq": freq})

    outname = make_plot(data, comp, imfile)
    plotdir = os.path.join(os.path.join(cutdir, 'plots'), os.path.basename(outname))
    copyfile(outname, plotdir)

def main(
    field,
    datadir,
    client,
    host,
    holofile,
    username=None,
    password=None,
    verbose=True,
    snr_cut=None,
):
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    beamseps = gen_seps(field, scriptdir)

    if datadir is not None:
        datadir = os.path.abspath(datadir)

    cutdir = os.path.abspath(os.path.join(datadir, "cutouts"))
    holofile = os.path.abspath(holofile)

    beams_col, island_col, comp_col = get_db(
        host=host, username=username, password=password
    )

    # Query the DB
    query = {
        "$and": [{f"beams.{field}": {"$exists": True}}, {f"beams.{field}.DR1": True}]
    }

    island_ids = sorted(beams_col.distinct("Source_ID", query))
    big_beams = list(beams_col.find(
        {"Source_ID": {'$in': island_ids}}).sort("Source_ID"))
    # files = sorted([name for name in glob(f"{cutdir}/*") if os.path.isdir(os.path.join(cutdir, name))])
    big_comps = list(comp_col.find(
        {"Source_ID": {'$in': island_ids}}).sort("Source_ID"))
    comps = []
    for island_id in island_ids:
        _comps = []
        for c in big_comps:
            if c["Source_ID"] == island_id:
                _comps.append(c)
        comps.append(_comps)

    assert len(big_beams) == len(comps)

    outputs = []
    for beams, comp in zip(big_beams, comps):
        for c in comp:
            if snr_cut is not None:
                noise = c['Noise']
                signal = c['Total_flux_Gaussian']
                snr_total = signal / noise
                if snr_total < snr_cut:
                    continue
            out = interpolate(
                field=field,
                comp=c,
                beams=beams,
                cutdir=cutdir,
                septab=beamseps,
                holofile=holofile
            )
            outputs.append(out)

    futures = chunk_dask(
        outputs=outputs,
        client=client,
        task_name="leakage plots",
        progress_text="Making leakage plots",
        verbose=verbose,
    )

    print("Comparing leakge done!")

def cli():
    """Command-line interface
    """
    import argparse
    from astropy.utils.exceptions import AstropyWarning

    warnings.simplefilter("ignore", category=AstropyWarning)
    from astropy.io.fits.verify import VerifyWarning

    warnings.simplefilter("ignore", category=VerifyWarning)
    warnings.simplefilter("ignore", category=RuntimeWarning)
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
    SPICE-RACS:
    Make leakage comparison plots.

    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.RawTextHelpFormatter
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

    parser.add_argument(
        "--holofile",
        type=str,
        help="Path to holography image"
    )

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

    parser.add_argument(
        "--snr", type=float, default=None, help="SNR cut (full band)."
    )

    args = parser.parse_args()

    cluster = LocalCluster(
        n_workers=12, threads_per_worker=1,
    )
    client = Client(cluster)

    main(
        field=args.field,
        datadir=args.datadir,
        client=client,
        host=args.host,
        holofile=args.holofile,
        username=args.username,
        password=args.password,
        verbose=args.verbose,
        snr_cut=args.snr,
    )

if __name__ == "__main__":
    cli()