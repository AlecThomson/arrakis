#!/usr/bin/env python3
"""Make a SPICE-RACS catalogue"""
import numpy as np
import warnings
from astropy.table import QTable, Column
from astropy.io import fits
import pymongo
from tqdm import tqdm, trange
from spiceracs import columns_possum
from spiceracs.utils import get_db, test_db
import rmtable as RMT
import logging as log
from pprint import pformat
from scipy.stats import lognorm, norm


def lognorm_from_percentiles(x1, p1, x2, p2):
    """ Return a log-normal distribuion X parametrized by:

            P(X < p1) = x1
            P(X < p2) = x2
    """
    x1 = np.log(x1)
    x2 = np.log(x2)
    p1ppf = norm.ppf(p1)
    p2ppf = norm.ppf(p2)

    scale = (x2 - x1) / (p2ppf - p1ppf)
    mean = ((x1 * p2ppf) - (x2 * p1ppf)) / (p2ppf - p1ppf)

#     return lognorm(s=scale, scale=np.exp(mean))
    return scale, np.exp(mean)


def sigma_add_fix(tab):
    sigma_Q_low = np.array(tab['Sigma_add_Q'] - tab['E_Sigma_add_minus_Q'])
    sigma_Q_high = np.array(tab['Sigma_add_Q'] + tab['E_Sigma_add_plus_Q'])

    sigma_U_low = np.array(tab['Sigma_add_U'] - tab['E_Sigma_add_minus_U'])
    sigma_U_high = np.array(tab['Sigma_add_U'] + tab['E_Sigma_add_plus_U'])

    s_Q, scale_Q = lognorm_from_percentiles(
        sigma_Q_low,
        15.72/100,
        sigma_Q_high,
        84.27/100
    )

    s_U, scale_U = lognorm_from_percentiles(
        sigma_U_low,
        15.72/100,
        sigma_U_high,
        84.27/100
    )

    med, std = np.zeros_like(s_Q), np.zeros_like(s_Q)
    for i, (_s_Q, _scale_Q, _s_U, _scale_U) in tqdm(
        enumerate(zip(s_Q, scale_Q, s_U, scale_U)),
        total=len(s_Q)
    ):
        try:
            Q_dist = lognorm.rvs(s=_s_Q, scale=_scale_Q, size=(1000))
            U_dist = lognorm.rvs(s=_s_U, scale=_scale_U, size=(1000))
            P_dist = np.hypot(Q_dist, U_dist)
            med[i] = np.median(P_dist)
            std[i] = np.std(P_dist)
        except ValueError:
            med[i] = np.nan
            std[i] = np.nan

    tab.add_column(Column(data=med, name='Sigma_add_P'))
    tab.add_column(Column(data=std, name='E_Sigma_add_P'))
    tab.remove_columns(
        [
            'Sigma_add_Q',
            'Sigma_add_U',
            'E_Sigma_add_minus_Q',
            'E_Sigma_add_plus_Q',
            'E_Sigma_add_minus_U',
            'E_Sigma_add_plus_U'
        ]
    )

    return tab


def is_leakage(frac, sep, fit):
    """Determine if a source is leakage

    Args:
        frac (float): Polarised fraction
        sep (float): Separation from tile centre
        fit (function): Fitting function

    Returns:
        bool: True if source is leakage
    """
    fit_frac = fit(sep)
    return frac < fit_frac


def get_fit_func(tab, nbins=21, offset=0.002):
    """Fit an envelope to define leakage sources

    Args:
        tab (Table): Catalogue to fit
        nbins (int, optional): Number of bins along seperation axis. Defaults to 21.

    Returns:
        np.polynomial.Polynomial.fit: 3rd order polynomial fit.
    """
    # Select high SNR sources
    hi_snr = (tab['Peak_I_flux_Gaussian'] / tab['Noise_I']) > 100
    hi_i_tab = tab[hi_snr]
    # Get fractional pol
    frac_P = np.array(hi_i_tab['Fractional_polarisation']) * \
        hi_i_tab['Fractional_polarisation'].unit
    # Bin sources by separation from tile centre
    bins = np.histogram_bin_edges(
        hi_i_tab['Separation_Tile_Centre'], bins=nbins)
    bins_c = np.median(np.vstack([bins[0:-1], bins[1:]]), axis=0)
    # Compute the median and standard deviation of the fractional pol
    meds = np.zeros_like(bins_c)
    s1_ups = np.zeros_like(bins_c)
    s1_los = np.zeros_like(bins_c)
    s2_ups = np.zeros_like(bins_c)
    s2_los = np.zeros_like(bins_c)
    for i in range(len(bins)-1):
        idx = ((hi_i_tab['Separation_Tile_Centre'] < bins[i+1])
               & (hi_i_tab['Separation_Tile_Centre'] >= bins[i]))
        s2_los[i], s1_los[i], meds[i], s1_ups[i], s2_ups[i] = np.nanpercentile(
            frac_P[idx], [2.3, 16, 50, 84, 97.6]
        )
    # Fit to median with small offset
    fit = np.polynomial.Polynomial.fit(bins_c, meds+offset, deg=3, full=False)
    return fit


def cuts_and_flags(cat):
    """Cut out bad sources, and add flag columns

    A flag of 'True' means the source is bad.

    Args:
        cat (RMT): Catalogue to cut and flag
    """
    # SNR cut
    snr_cut = cat['SNR_Polarised_intensity'] > 5
    cat = cat[snr_cut]
    # SNR flag
    snr_flag = cat['SNR_Polarised_intensity'] < 8
    cat.add_column(Column(data=snr_flag, name='SNR_flag'))
    # Leakage flag
    fit = get_fit_func(cat)
    leakage_flag = is_leakage(
        cat['Fractional_polarisation'],
        cat['Separation_Tile_Centre'],
        fit
    )
    cat.add_column(Column(data=leakage_flag, name='Leakage_flag'))
    # Channel flag
    chan_flag = cat['Number_of_channels'] < 144
    cat.add_column(Column(data=chan_flag, name='Channel_flag'))
    # Fitting flag
    fit_flag = cat['Stokes_I_fit_flag'] >= 64
    cat.remove_column('Stokes_I_fit_flag')
    cat.add_column(Column(data=fit_flag, name='Stokes_I_fit_flag'))
    # Sigma_add flag
    sigma_flag = cat['Sigma_add_P'] > 1
    cat.add_column(Column(data=sigma_flag, name='Complex_sigma_add_flag'))
    # M2_CC flag
    m2_flag = cat['M2_CC'] > cat['RMSF_FWHM']
    cat.add_column(Column(data=m2_flag, name='Complex_M2_CC_flag'))

    return cat


def main(
    field: str,
    host: str,
    username: str = None,
    password: str = None,
    verbose=True,
    outfile: str = None,
    cat_format: str = None,
) -> None:
    """Main

    Args:
        field (str): RACS field name
        host (str): MongoDB host IP
        username (str, optional): Mongo username. Defaults to None.
        password (str, optional): Mongo password. Defaults to None.
        verbose (bool, optional): Verbose output. Defaults to True.
        outfile (str, optional): Output file name. Defaults to None.
        cat_format (str, optional): Type of catalogue .e.g. fits. Defaults to None.
    """
    # default connection (ie, local)
    beams_col, island_col, comp_col = get_db(
        host=host, username=username, password=password
    )
    query = {
        "$and": [{f"beams.{field}": {"$exists": True}}, {f"beams.{field}.DR1": True}]
    }

    beams = beams_col.find(query).sort("Source_ID")
    all_island_ids = sorted(beams_col.distinct("Source_ID", query))
    query = {
        "$and": [{"Source_ID": {"$in": all_island_ids}}, {"rmclean1d": True}]}

    fields = {}
    for n in columns_possum.input_names:
        fields.update({n: 1})
    for n in columns_possum.sourcefinder_columns:
        fields.update({n: 1})
    fields.update({"rmsynth_summary": 1})
    fields.update({"rmclean_summary": 1})
    fields.update({"header": 1})

    comps = list(
        comp_col.find(
            query,
            fields
        )
    )

    # tab = RMT.RMTable()
    tab = QTable()
    # Add items to main cat using RMtable standard
    for j, [name, typ, src, col, unit] in enumerate(
        tqdm(
            zip(
                columns_possum.output_cols,
                columns_possum.output_types,
                columns_possum.input_sources,
                columns_possum.input_names,
                columns_possum.output_units,
            ),
            total=len(columns_possum.output_cols),
            desc="Making table by column",
            disable=not verbose,
        ),
    ):
        data = []
        if src == "cat":
            for comp in comps:
                data += [comp[col]]
            new_col = Column(data=data, name=name, dtype=typ, unit=unit)
            tab.add_column(new_col)

        if src == "synth":
            for comp in comps:
                try:
                    data += [comp["rmclean_summary"][col]]
                except KeyError:
                    data += [comp["rmsynth_summary"][col]]
            new_col = Column(data=data, name=name, dtype=typ, unit=unit)
            tab.add_column(new_col)

        if src == "header":
            for comp in comps:
                data += [comp["header"][col]]
            new_col = Column(data=data, name=name, dtype=typ, unit=unit)
            tab.add_column(new_col)

    for selcol in tqdm(columns_possum.sourcefinder_columns, desc="Adding BDSF data"):
        data = []
        for comp in comps:
            data += [comp[selcol]]
        new_col = Column(data=data, name=selcol)
        tab.add_column(new_col)
    rmtab = RMT.from_table(tab)
    # Get Galatic coords
    rmtab["Gal_lon"], rmtab["Gal_lat"] = RMT.calculate_missing_coordinates_column(
        rmtab["RA"], rmtab["Dec"], to_galactic=True
    )

    # Add common columns
    rmtab["rm_method"] = "RM Synthesis"
    rmtab["telescope"] = "ASKAP"
    rmtab["pol_bias"] = "2012PASA...29..214G"

    # Verify table
    rmtab.verify_columns()
    # rmtab.verify_standard_strings()
    rmtab.verify_limits()
    rmtab.add_missing_columns()

    # Fix sigma_add
    rmtab = sigma_add_fix(rmtab)

    # Add flags
    rmtab = cuts_and_flags(rmtab)

    if outfile is None:
        log.info(pformat(rmtab))

    if outfile is not None:
        rmtab.table.write(outfile, format=cat_format, overwrite=True)
        log.info(f"{outfile} written to disk")

    log.info("Done!")


def cli():
    """Command-line interface"""
    import argparse
    from astropy.utils.exceptions import AstropyWarning

    warnings.simplefilter("ignore", category=AstropyWarning)
    from astropy.io.fits.verify import VerifyWarning

    warnings.simplefilter("ignore", category=VerifyWarning)
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
    SPICE-RACS Stage 7:
    Make RM catalogue.

    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "field", metavar="field", type=str, help="RACS field to mosaic - e.g. 2132-50A."
    )

    parser.add_argument(
        "host",
        metavar="host",
        type=str,
        help="Host of mongodb (probably $hostname -i).",
    )

    parser.add_argument(
        "--username", type=str, default=None, help="Username of mongodb."
    )

    parser.add_argument(
        "--password", type=str, default=None, help="Password of mongodb."
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="verbose output [False]."
    )

    parser.add_argument(
        "-w",
        "--write",
        dest="outfile",
        default=None,
        type=str,
        help="File to save table to [None].",
    )

    parser.add_argument(
        "-f",
        "--format",
        dest="format",
        default=None,
        type=str,
        help="Format for output file [None].",
    )

    args = parser.parse_args()
    if args.outfile and not args.format:
        parser.error("Please provide an output file format.")

    verbose = args.verbose

    if verbose:
        log.basicConfig(
            level=log.INFO,
            format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            force=True
        )
    else:
        log.basicConfig(
            format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            force=True
        )

    host = args.host
    test_db(
        host=args.host, username=args.username, password=args.password, verbose=verbose
    )

    main(
        field=args.field,
        host=host,
        username=args.username,
        password=args.password,
        verbose=verbose,
        outfile=args.outfile,
        cat_format=args.format,
    )


if __name__ == "__main__":
    cli()
