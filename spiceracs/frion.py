#!/usr/bin/env python3
import os
import numpy as np
import astropy.units as u
from astropy.time import Time, TimeDelta
from spiceracs.utils import get_field_db, getfreq, get_db, tqdm_dask, test_db, try_mkdir
from FRion import predict, correct
from dask import delayed
from dask.distributed import Client, progress, LocalCluster, wait
import pymongo
import time
from pprint import pprint
from shutil import copyfile
from glob import glob


@delayed
def correct_worker(beam, outdir, field, predict_file, island_id):
    """FR correction

    Args:
        beam (dict): Beam database
        outdir (str): Directory containing cutouts
        field (str): RACS field
        predict_file (str): Name of predict file to read
        island_id (str): Source ID

    Returns:
        UpdateOne: Mongo update
    """    
    qfile = os.path.join(outdir, beam['beams'][field]['q_file'])
    ufile = os.path.join(outdir, beam['beams'][field]['u_file'])

    qout = beam['beams'][field]['q_file'].replace('.fits', '.ion.fits')
    uout = beam['beams'][field]['u_file'].replace('.fits', '.ion.fits')

    qout_f = os.path.join(outdir, qout)
    uout_f = os.path.join(outdir, uout)

    correct.apply_correction_to_files(
        qfile, ufile, predict_file, qout_f, uout_f, overwrite=True
    )

    myquery = {"Source_ID": island_id}

    newvalues = {
        '$set': {
            f'beams.{field}.q_file_ion': qout,
            f'beams.{field}.u_file_ion': uout,
        }
    }
    return pymongo.UpdateOne(myquery, newvalues)


@delayed
def predict_worker(island, field, beam, start_time, end_time, freq, cutdir, plotdir):
    """FR prediction

    Args:
        island (dict): Island databse entry
        field (str): RACS field
        beam (dict): Beam database entry
        start_time (Time): Obs start time
        end_time (Time): Obs end time
        freq (array): Frequency array
        cutdir (str): Directory containing cutouts

    Returns:
        str: Location of prediction file
    """    
    ifile = os.path.join(cutdir, beam["beams"][field]["i_file"])
    i_dir = os.path.dirname(ifile)
    iname = island["Source_ID"]
    ra = island['RA']
    dec = island['Dec']

    times, RMs, theta = predict.calculate_modulation(
        start_time=start_time.fits,
        end_time=end_time.fits,
        freq_array=freq,
        telescope_location=predict.get_telescope_coordinates('ASKAP'),
        ra=ra,
        dec=dec,
        timestep=300.0,
        ionexPath=os.path.join(os.path.dirname(cutdir), 'IONEXdata')
    )
    predict_file = os.path.join(i_dir, f"{iname}_ion.txt")
    predict.write_modulation(
        freq_array=freq,
        theta=theta,
        filename=predict_file
    )

    plot_file = os.path.join(i_dir, f'{iname}_ion.pdf')
    predict.generate_plots(
        times,
        RMs,
        theta,
        freq,
        position=[ra, dec],
        savename=plot_file
    )
    plot_files = glob(os.path.join(i_dir, '*ion.pdf'))
    for src in plot_files:
        base = os.path.basename(src)
        dst = os.path.join(plotdir, base)
        copyfile(src, dst)
    return predict_file


def main(
    field,
    outdir,
    host,
    client,
    username=None,
    password=None,
    database=False,
    verbose=True
):
    """Main script

    Args:
        field (str): RACS field
        outdir (str): Dir containing 'cutouts'
        host (str): MongoDB host
        client (client): Dask client
        username (str, optional): MongoDB username. Defaults to None.
        password (str, optional): Mongodb password. Defaults to None.
        verbose (bool, optional): Verbose output. Defaults to True.
    """
    # Query database for data
    outdir = os.path.abspath(outdir)
    cutdir = os.path.join(outdir, "cutouts")

    plotdir = os.path.join(cutdir, 'plots')
    try_mkdir(plotdir)

    beams_col, island_col, comp_col = get_db(
        host=host, username=username, password=password
    )

    query = {
        "$and": [{f"beams.{field}": {"$exists": True}}, {f"beams.{field}.DR1": True}]
    }

    beams = list(beams_col.find(query).sort("Source_ID"))
    island_ids = sorted(beams_col.distinct("Source_ID", query))

    # Get FRion arguments
    query = {"Source_ID": {"$in": island_ids}}
    islands = list(island_col.find(query).sort("Source_ID"))

    field_col = get_field_db(
        host,
        username=username,
        password=password
    )
    query = {"FIELD_NAME": f"RACS_{field}"}
    # Get most recent SBID
    if field_col.count_documents(query) > 1:
        field_datas = list(field_col.find({"FIELD_NAME": f"RACS_{field}"}))
        sbids = [f['CAL_SBID'] for f in field_datas]
        max_idx = np.argmax(sbids)
        if verbose:
            print(f"Using CAL_SBID {sbids[max_idx]}")
        field_data = field_datas[max_idx]
    else:
        field_data = field_col.find_one({"FIELD_NAME": f"RACS_{field}"})

    start_time = Time(field_data['SCAN_START']*u.second, format='mjd')
    end_time = start_time + TimeDelta(field_data['SCAN_TINT']*u.second)

    freq = getfreq(
        os.path.join(cutdir, f"{beams[0]['beams'][f'{field}']['q_file']}"),
        verbose=verbose,
    )

    
    # Loop over islands in parallel
    outputs = []
    for island in islands:
        island_id = island['Source_ID']
        beam_idx = [i for i, b in enumerate(
            beams) if b['Source_ID'] == island_id][0]
        beam = beams[beam_idx]
        # Get FRion predictions
        predict_file = predict_worker(
            island=island,
            field=field,
            beam=beam,
            start_time=start_time,
            end_time=end_time,
            freq=freq.to(u.Hz).value,
            cutdir=cutdir,
            plotdir=plotdir,
        )
        # Apply FRion predictions
        output = correct_worker(
            beam=beam,
            outdir=cutdir,
            field=field,
            predict_file=predict_file,
            island_id=island_id
        )
        outputs.append(output)

    # Wait for IONEX data I guess...
    _ = outputs[0].compute()
    time.sleep(10)
    # Execute
    futures = client.persist(outputs)
    # dumb solution for https://github.com/dask/distributed/issues/4831
    time.sleep(10)
    tqdm_dask(futures, desc="Running FRion", disable=(not verbose), total=len(islands)*2)
    if database:
        if verbose:
            print("Updating database...")
        updates = [f.compute() for f in futures]
        db_res = beams_col.bulk_write(
            updates,
            ordered=False
            )
        if verbose:
            pprint(db_res.bulk_api_result)


def cli():
    """Command-line interface
    """
    import argparse
    import warnings
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
    SPICE-RACS Stage:
    Correct for ionospheric Faraday rotation

    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "field", metavar="field", type=str, help="RACS field to mosaic - e.g. 2132-50A."
    )
    parser.add_argument(
        "outdir",
        metavar="outdir",
        type=str,
        help="Directory containing cutouts (in subdir outdir/cutouts).",
    )

    parser.add_argument(
        "host",
        metavar="host",
        type=str,
        help="Host of mongodb (probably $hostname -i).",
    )

    parser.add_argument(
        "--username", type=str, default='admin', help="Username of mongodb."
    )

    parser.add_argument(
        "--password", type=str, default=None, help="Password of mongodb."
    )
    
    parser.add_argument(
        "-m", "--database", action="store_true", help="Add data to MongoDB [False]."
    )

    parser.add_argument(
        "-v", dest="verbose", action="store_true", help="verbose output [False]."
    )


    args = parser.parse_args()

    verbose = args.verbose

    cluster = LocalCluster(
        n_workers=12, processes=True, threads_per_worker=1, local_directory="/dev/shm"
    )
    client = Client(cluster)
    print(client)

    test_db(host=args.host, username=args.username,
            password=args.password, verbose=verbose)

    main(
        field=args.field,
        outdir=args.outdir,
        host=args.host,
        client=client,
        username=args.username,
        password=args.password,
        database=args.database,
        verbose=verbose
    )


if __name__ == "__main__":
    cli()
