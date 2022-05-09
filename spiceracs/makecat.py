#!/usr/bin/env python3
"""Make a SPICE-RACS catalogue"""
import os
import time
import numpy as np
import warnings
from astropy.table import QTable, Column
from astropy.io import fits
from astropy.io import votable as vot
import astropy.units as u
from tqdm import tqdm, trange
from spiceracs import columns_possum
from spiceracs.utils import get_db, test_db, get_field_db, latexify
import rmtable as rmt
import logging as log
from pprint import pformat
from scipy.stats import lognorm, norm
import matplotlib.pyplot as plt
from IPython import embed



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

    return scale, np.exp(mean)


def sigma_add_fix(tab):
    sigma_Q_low = np.array(tab['sigma_add_Q'] - tab['sigma_add_Q_err_minus'])
    sigma_Q_high = np.array(tab['sigma_add_Q'] + tab['sigma_add_Q_err_plus'])

    sigma_U_low = np.array(tab['sigma_add_U'] - tab['sigma_add_U_err_minus'])
    sigma_U_high = np.array(tab['sigma_add_U'] + tab['sigma_add_U_err_plus'])

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
        total=len(s_Q),
        desc="Calculating sigma_add"
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

    tab.add_column(Column(data=med, name='sigma_add'))
    tab.add_column(Column(data=std, name='sigma_add_err'))
    tab.remove_columns(
        [
            'sigma_add_Q',
            'sigma_add_U',
            'sigma_add_Q_err_minus',
            'sigma_add_Q_err_plus',
            'sigma_add_U_err_minus',
            'sigma_add_U_err_plus'
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


def get_fit_func(tab, nbins=21, offset=0.002, degree=2, do_plot=False):
    """Fit an envelope to define leakage sources

    Args:
        tab (Table): Catalogue to fit
        nbins (int, optional): Number of bins along seperation axis. Defaults to 21.

    Returns:
        np.polynomial.Polynomial.fit: 3rd order polynomial fit.
    """
    # Select high SNR sources
    hi_snr = (tab['stokesI'].to(u.Jy/u.beam) / tab['stokesI_err'].to(u.Jy/u.beam)) > 100
    hi_i_tab = tab[hi_snr]
    # Get fractional pol
    frac_P = np.array(hi_i_tab['fracpol'].value)
    # Bin sources by separation from tile centre
    bins = np.histogram_bin_edges(
        hi_i_tab['beamdist'].to(u.deg).value, bins=nbins)
    bins_c = np.median(np.vstack([bins[0:-1], bins[1:]]), axis=0)
    # Compute the median and standard deviation of the fractional pol
    meds = np.zeros_like(bins_c)
    s1_ups = np.zeros_like(bins_c)
    s1_los = np.zeros_like(bins_c)
    s2_ups = np.zeros_like(bins_c)
    s2_los = np.zeros_like(bins_c)
    for i in range(len(bins)-1):
        idx = ((hi_i_tab['beamdist'].to(u.deg).value < bins[i+1])
               & (hi_i_tab['beamdist'].to(u.deg).value >= bins[i]))
        s2_los[i], s1_los[i], meds[i], s1_ups[i], s2_ups[i] = np.nanpercentile(
            frac_P[idx], [2.3, 16, 50, 84, 97.6]
        )
    # Fit to median with small offset
    fit = np.polynomial.Polynomial.fit(
        bins_c, 
        meds+offset, 
        deg=degree, 
        full=False
    )
    if not do_plot: 
        return fit
    
    # Plot the fit
    latexify(columns=1)
    figure = plt.figure(facecolor='w')
    fig = plt.figure(facecolor='w')
    color = 'tab:green'
    stoke = {
        "s2_los": s2_los,
        "s1_los": s1_los,
        "meds": meds,
        "s1_ups": s1_ups,
        "s2_ups" :s2_ups,
    }
    plt.scatter(
        hi_i_tab['beamdist'].to(u.deg).value,
        frac_P, 
        s=1, 
        alpha=0.2,
        marker='.',
        c='k',
        zorder=0,
    )
    plt.plot(
        bins_c,
        meds,
        alpha=1,
        c=color,
        label="Median"
    )
    for s, ls in zip((1,2), ("--",":")):
        for r in ("ups", "los"):
            plt.plot(
                bins_c,
                stoke[f"s{s}_{r}"],
                alpha=1,
                c=color,
                linestyle=ls,
                label=f"${s}\sigma$" if r=="ups" else ""
            )
    xx = np.linspace(0, 4.5, 100)
    plt.plot(xx, fit(xx), 'tab:orange', label="Leakage envelope")
    plt.legend(loc='upper left')
    plt.xlabel('Separation from tile centre [deg]')
    plt.ylabel(f'$P/I$ fraction')
    plt.ylim(0,+0.05)
    plt.grid()
    return fit, fig


def cuts_and_flags(cat):
    """Cut out bad sources, and add flag columns

    A flag of 'True' means the source is bad.

    Args:
        cat (rmt): Catalogue to cut and flag
    """
    # SNR flag
    snr_flag = cat['snr_polint'] < 8
    cat.add_column(Column(data=snr_flag, name='snr_flag'))
    # Leakage flag
    fit = get_fit_func(cat)
    leakage_flag = is_leakage(
        cat['fracpol'].value,
        cat['beamdist'].to(u.deg).value,
        fit
    )
    cat.add_column(Column(data=leakage_flag, name='leakage_flag'))
    # Channel flag
    chan_flag = cat['Nchan'] < 144
    cat.add_column(Column(data=chan_flag, name='channel_flag'))
    # Fitting flag
    fit_flag = cat['stokes_I_fit_flag'] >= 64
    cat.remove_column('stokes_I_fit_flag')
    cat.add_column(Column(data=fit_flag, name='stokes_I_fit_flag'))
    # sigma_add flag
    sigma_flag = cat['sigma_add'] > 1
    cat.add_column(Column(data=sigma_flag, name='complex_sigma_add_flag'))
    # M2_CC flag
    m2_flag = cat['rm_width'] > cat['rmsf_fwhm']
    cat.add_column(Column(data=m2_flag, name='complex_M2_CC_flag'))

    return cat, fit

def get_alpha(cat):
    coefs_str = cat["stokesI_model_coef"]
    alphas = []
    for c in coefs_str:
        coefs = c.split(",")
        alpha = float(coefs[-2]) # alpha is the 2nd last coefficient
        alphas.append(alpha)
    return np.array(alphas)

def get_integration_time(cat, field_col):
    field_names = list(cat['tile_id'])
    query = {
        "$and": [
            {"FIELD_NAME": {"$in": field_names}},
            {"SELECT": 1}
        ]
    }
    tint_dicts = list(field_col.find(query, {"_id":0,"SCAN_TINT": 1, "FIELD_NAME": 1}))
    tint_dict = {}
    for d in tint_dicts:
        tint_dict.update(
                {
                    d["FIELD_NAME"]: d["SCAN_TINT"]
                }
        )

    tints = []
    for name in field_names:
        tints.append(tint_dict[name])
    
    return np.array(tints) * u.s

# Stolen from GASKAP pipeline
# Credit to J. Dempsey
# https://github.com/GASKAP/GASKAP-HI-Absorption-Pipeline/
# https://github.com/GASKAP/GASKAP-HI-Absorption-Pipeline/blob/
# def add_col_metadata(vo_table, col_name, description, units=None, ucd=None, datatype=None):
#     """Add metadata to a VO table column.

#     Args:
#         vo_table (vot.): VO Table
#         col_name (str): Column name
#         description (str): Long description of the column
#         units (u.Unit, optional): Unit of column. Defaults to None.
#         ucd (str, optional): UCD string. Defaults to None.
#         datatype (_type_, optional): _description_. Defaults to None.
#     """    
#     col = vo_table.get_first_table().get_field_by_id(col_name)
#     col.description = description
#     if units:
#         col.unit = units
#     if ucd:
#         col.ucd = ucd
#     if datatype:
#         col.datatype = datatype

def add_metadata(vo_table: vot.tree.Table, filename: str):
    """Add metadata to VO Table for CASDA

    Args:
        vo_table (vot): VO Table object

    Returns:
        vot: VO Table object with metadata
    """    
    # Add RMtable metadata
    for col_name, desc in columns_possum.rmtab_column_descriptions.items():
        col = vo_table.get_first_table().get_field_by_id(col_name)
        col.description = desc
        if col_name == "ra":
            col.ucd = "pos.eq.ra;meta.main"
        elif col_name == "dec":
            col.ucd = "pos.eq.dec;meta.main"
        elif col_name == "cat_id":
            col.ucd = "meta.id;meta.main"
        elif col_name == "source_id":
            col.ucd = "meta.id"
    # Add extra metadata
    for col_name, desc in columns_possum.extra_column_descriptions.items():
        col = vo_table.get_first_table().get_field_by_id(col_name)
        col.description = desc

    # Add params for CASDA
    if len(vo_table.params) > 0:
        log.warning(f"{filename} already has params - not adding")
        return vo_table
    _ , ext = os.path.splitext(filename)
    cat_name = os.path.basename(filename).replace(ext, "")
    idx_fields = "ra,dec,cat_id,source_id"
    pri_fields = "ra,dec,cat_id,source_id,rm,polint,snr_polint,fracpol,stokesI,sigma_add"
    params = [
        vot.tree.Param(
            vo_table,
            ID="Catalogue_Name", 
            name="Catalogue Name", 
            value=cat_name
        ),
        vot.tree.Param(
            vo_table,
            ID="Indexed_Fields", 
            name="Indexed Fields", 
            value=idx_fields
        ),
        vot.tree.Param(
            vo_table,
            ID="Principal_Fields", 
            name="Principal Fields", 
            value=pri_fields
        ),
    ]
    vo_table.params.extend(params)

    return vo_table

def replace_nans(filename:str):
    """Replace NaNs in a XML table with a string

    Args:
        filename (str): File name
    """  
    pass
    # with open(filename, "r") as f:
    #     xml = f.read()
    # xml = xml.replace("NaN", "null")
    # with open(filename, "w") as f:
    #     f.write(xml)

def write_votable(rmtab:rmt.RMTable,outfile:str) -> None:
    # CASDA needs v1.3
    vo_table = vot.from_table(rmtab.table)
    vo_table.version = "1.3"
    vo_table = add_metadata(vo_table, outfile)
    vot.writeto(vo_table, outfile)
    # Fix NaNs for CASDA
    replace_nans(outfile)

def main(
    field: str,
    host: str,
    username: str = None,
    password: str = None,
    verbose=True,
    outfile: str = None,
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
    log.info("Starting beams collection query")
    tick = time.time()
    query = {
        "$and": [{f"beams.{field}": {"$exists": True}}, {f"beams.{field}.DR1": True}]
    }
    all_island_ids = sorted(beams_col.distinct("Source_ID", query))
    tock = time.time()
    log.info(f"Finished beams collection query - {tock-tick:.2f}s")

    log.info("Starting component collection query")
    tick = time.time()
    query = {
        "$and": [
            {"Source_ID": {"$in": all_island_ids}}, 
            {"rmsynth1d": True},
            {"rmclean1d": True},
        ]
    }

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
    tock = time.time()
    log.info(f"Finished component collection query - {tock-tick:.2f}s")

    # tab = rmt.rmtable()
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
    
    # Fix sigma_add
    tab = sigma_add_fix(tab)

    # Add flags
    tab, fit = cuts_and_flags(tab)

    # Add spectral index from fitted model
    alphas = get_alpha(tab)
    tab.add_column(Column(data=alphas, name='spectral_index'))

    # Add integration time
    field_col = get_field_db(
        host=host, username=username, password=password
    )
    tints = get_integration_time(tab, field_col)
    tab.add_column(Column(data=tints, name='int_time'))
    # Add epoch
    tab.add_column(Column(data=tab['start_time'] + (tints / 2), name='epoch'))

    # Convert to rmtab
    rmtab = rmt.from_table(tab)
    # Get Galatic coords
    glon, glat = rmt.calculate_missing_coordinates_column(
        rmtab["ra"], rmtab["dec"], to_galactic=True
    )
    rmtab.add_column(data=glon, name="l")
    rmtab.add_column(data=glat, name="b")
    rmtab.add_column(data=np.max([rmtab['ra_err'], rmtab['dec_err']]), name="pos_err")

    # Add common columns
    rmtab["rm_method"] = "RM Synthesis - Fractional polarization"
    rmtab["telescope"] = "ASKAP"
    rmtab["pol_bias"] = "2012PASA...29..214G"
    rmtab["catalog"] = "SPICE-RACS-DR1"
    rmtab["ionosphere"] = "FRion"
    rmtab["flux_type"] = "Peak"
    rmtab["aperture"] = 0*u.deg

    rmtab.add_column(
        data=fit(
            rmtab["separation_tile_centre"].to(u.deg).value,
        ),
        name="leakage"
    )

    rmtab.add_column(
        data=np.logical_or(
            rmtab['complex_sigma_add_flag'], 
            rmtab['complex_M2_CC_flag']
        ), 
        name='complex_flag'
    )

    # Verify table
    rmtab.verify_columns()
    rmtab.add_missing_columns()
    # Readd complex test
    rmtab["complex_test"]  = "sigma_add OR Second moment"


    if outfile is None:
        log.info(pformat(rmtab))

    if outfile is not None:
        _ , ext = os.path.splitext(outfile)
        if ext == ".xml" or ext == ".vot":
            write_votable(rmtab, outfile)
        else:
            rmtab.table.write(outfile, overwrite=True)
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

    args = parser.parse_args()

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
    )


if __name__ == "__main__":
    cli()
