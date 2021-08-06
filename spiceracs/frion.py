#!/usr/bin/env python3
import os
import numpy as np
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord
from spiceracs.utils import get_field_db, getfreq, get_db
from FRion import predict, correct
from dask import delayed
import pymongo
import time
from pprint import pprint


@delayed
def correct_worker(beam, outdir, field, predict_file, island_id):
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
            'beams.field.q_file_ion': qout,
            'beams.field.u_file_ion': uout,
        }
    }
    return pymongo.UpdateOne(myquery, newvalues)

def main(
    field,
    outdir,
    host,
    client,
    username=None,
    password=None,
    verbose=True
):

    # Query database for data
    outdir = os.path.abspath(outdir)
    cutdir = os.path.join(outdir, "cutouts")

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

    field_col = get_field_db(
        host,
        username=username,
        password=password
    )
    query = {"FIELD_NAME": f"RACS_{field}"}
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

    # Get FRion predictions
    times, RMs, theta = predict.calculate_modulation(
        start_time=start_time.fits,
        end_time=end_time.fits,
        freq_array=freq.to(u.Hz).value,
        telescope_location=predict.get_telescope_coordinates('ASKAP'),
        ra=field_data['RA_DEG'],
        dec=field_data['DEC_DEG'],
        timestep=300.0,
        ionexPath=os.path.join(outdir, 'IONEXdata')
    )

    predict_file = os.path.join(outdir, f"{field}_ion.txt")
    predict.write_modulation(
        freq_array=freq, theta=theta, filename=predict_file)
    plot_file = os.path.join(cutdir, 'plots', f'{field}_ion.pdf')
    predict.generate_plots(times, RMs, theta, freq, position=[
                           field_data['RA_DEG'], field_data['DEC_DEG']], savename=plot_file)

    # Apply FRion predictions
    outputs = []
    for i, island_id in enumerate(island_ids):
        beam_idx = [i for i, b in enumerate(
            beams) if b['Source_ID'] == island_id][0]
        beam = beams[beam_idx]
        output = correct_worker(
            beam=beam,
            outdir=cutdir,
            field=field,
            predict_file=predict_file,
            island_id=island_id
        )
        outputs.append(output)

    futures = client.persist(outputs)
    # dumb solution for https://github.com/dask/distributed/issues/4831
    time.sleep(5)
    tqdm_dask(futures, desc="Running FRion", disable=(not verbose))
    if verbose:
        print("Updating database...")
    updates = [f.compute() for f in futures]
    db_res = beams_col.bulk_write(updates)
        if verbose:
            pprint(db_res.bulk_api_result)