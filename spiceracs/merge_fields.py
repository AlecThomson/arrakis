import os
from shutil import copyfile
import pymongo
from spiceracs.utils import getfreq, MyEncoder, test_db, tqdm_dask, try_mkdir, get_db
from tqdm import tqdm
from typing import List, Tuple, Dict
from dask import delayed
from dask.distributed import Client, LocalCluster
from spiceracs.linmos import linmos, get_yanda


@delayed
def copy_singleton(
        beam: dict,
        vals: dict,
        merge_name: str,
        field_dir: str,
        data_dir: str
) -> pymongo.UpdateOne:
    try:
        i_file_old = os.path.join(field_dir, vals["i_file"])
        q_file_old = os.path.join(field_dir, vals["q_file_ion"])
        u_file_old = os.path.join(field_dir, vals["u_file_ion"])
    except KeyError:
        raise KeyError("Ion files not found. Have you run FRion?")
    new_dir = os.path.join(
        data_dir,
        beam["Source_ID"]
    )

    try_mkdir(
        new_dir
    )

    i_file_new = os.path.join(new_dir, os.path.basename(
        i_file_old)).replace('.fits', '.edge.linmos.fits')
    q_file_new = os.path.join(new_dir, os.path.basename(
        q_file_old)).replace('.fits', '.edge.linmos.fits')
    u_file_new = os.path.join(new_dir, os.path.basename(
        u_file_old)).replace('.fits', '.edge.linmos.fits')

    for src, dst in zip([i_file_old, q_file_old, u_file_old], [i_file_new, q_file_new, u_file_new]):
        copyfile(i_file_old, dst)

    query = {"Source_ID": beam["Source_ID"]}
    newvalues = {
        "$set": {
            f"beams.{merge_name}.i_file": i_file_new,
            f"beams.{merge_name}.q_file": q_file_new,
            f"beams.{merge_name}.u_file": u_file_new,
        }
    }

    return pymongo.UpdateOne(query, newvalues)


def copy_singletons(
        field_dict: Dict[str,str],
        data_dir: str,
        beams_col: pymongo.collection.Collection,
        merge_name: str
) -> list:
    # Find all islands with the given fields that DON'T overlap another field
    query = {
        "$or": [
            {
                "$and": [
                    {f"beams.{field}": {"$exists": True}},
                    {f"beams.{field}.DR1": True},
                    {'n_fields_DR1': 1}
                ]
            } for field in field_dict.keys()
        ]
    }

    island_ids = sorted(beams_col.distinct("Source_ID", query))
    big_beams = list(
        beams_col.find({"Source_ID": {"$in": island_ids}}).sort("Source_ID")
    )
    updates = []
    for beam in tqdm(big_beams):
        for field, vals in beam["beams"].items():
            if field not in field_dict.keys():
                continue
            field_dir = field_dict[field]
            update = copy_singleton(
                beam, vals, merge_name, field_dir, data_dir)
            updates.append(update)
    return updates


@delayed
def genparset(
    old_ims: list,
    stokes: str,
    new_dir: str,

) -> str:


    imlist = "[" + ','.join([im.replace('.fits', '') for im in old_ims]) + "]"
    weightlist = f"[{','.join([im.replace('.fits', '').replace('.image.restored.','.weights.') for im in old_ims])}]"

    im_outname = os.path.join(new_dir, os.path.basename(
        old_ims[0])).replace('.fits', '.egde.linmos')
    wt_outname = os.path.join(new_dir, os.path.basename(
        old_ims[0])).replace('.fits', '.egde.linmos').replace('.image.restored.','.weights.') 


    parset_file = os.path.join(new_dir, f"edge_linmos_{stokes}.in")
    parset = f"""# LINMOS parset
linmos.names            = {imlist}
linmos.weights          = {weightlist}
linmos.imagetype        = fits
linmos.outname          = {im_outname}
linmos.outweight        = {wt_outname}
# For ASKAPsoft>1.3.0
linmos.weighttype       = FromWeightImages
linmos.weightstate      = Corrected
"""

    with open(parset_file, "w") as f:
        f.write(parset)

    return parset_file


def merge_multiple_fields(
        field_dict: Dict[str,str],
        data_dir: str,
        beams_col: pymongo.collection.Collection,
        merge_name: str,
        image: str,
) -> list:

    # Find all islands with the given fields that overlap another field
    query = {
        "$or": [
            {
                "$and": [
                    {f"beams.{field}": {"$exists": True}},
                    {f"beams.{field}.DR1": True},
                    {'n_fields_DR1': {"$gt": 1}}
                ]
            } for field in field_dict.keys()
        ]
    }

    island_ids = sorted(beams_col.distinct("Source_ID", query))
    big_beams = list(
        beams_col.find({"Source_ID": {"$in": island_ids}}).sort("Source_ID")
    )

    # Sanity check that the FRion has been applied
    # "beams.{field}.q_file_ion"
    updates = []
    for beam in tqdm(big_beams):
        i_files_old = []
        q_files_old = []
        u_files_old = []
        for field, vals in beam["beams"].items():
            if field not in field_dict.keys():
                continue
            field_dir = field_dict[field]
            try:
                i_file_old = os.path.join(field_dir, vals["i_file"])
                q_file_old = os.path.join(field_dir, vals["q_file_ion"])
                u_file_old = os.path.join(field_dir, vals["u_file_ion"])
            except KeyError:
                raise KeyError("Ion files not found. Have you run FRion?")
            i_files_old.append(i_file_old)
            q_files_old.append(q_file_old)
            u_files_old.append(u_file_old)

        new_dir = os.path.join(
            data_dir,
            beam["Source_ID"]
        )

        try_mkdir(
            new_dir
        )

        for stokes, imlist in zip(['I', 'Q', 'U'], [i_files_old, q_files_old, u_files_old]):
            parset_file = genparset(imlist, stokes, new_dir)
            update = linmos(parset_file, merge_name, image)
            updates.append(update)

    return updates


def main(
    fields: List[str],
    field_dirs: List[str],
    merge_name: str,
    output_dir: str,
    host: str = None,
    username: str = None,
    password: str = None,
    yanda="1.3.0",
    verbose: bool = True,
) -> None:

    assert len(fields) == len(
        field_dirs), f"List of fields must be the same length as length of field dirs."

    field_dict ={field: field_dir for field, field_dir in zip(fields, field_dirs)}

    image = get_yanda(version=yanda)

    beams_col, island_col, comp_col = get_db(
        host=host, username=username, password=password
    )

    output_dir = os.path.abspath(output_dir)
    data_dir = os.path.join(output_dir, merge_name)
    try_mkdir(data_dir)


    singleton_updates = copy_singletons(
        field_dict=field_dict,
        data_dir=data_dir,
        beams_col=beams_col,
        merge_name=merge_name
    )

    mutilple_updates = merge_multiple_fields(
        field_dict=field_dict,
        data_dir=data_dir,
        beams_col=beams_col,
        merge_name=merge_name,
        image=image,
    )

    from IPython import embed
    embed()


def cli():
    """Command-line interface"""
    import argparse

    # Help string to be shown using the -h option
    descStr = f"""
    Mosaic RACS beam fields with linmos.

    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "--fields",
        type=str,
        nargs='+',
        help="RACS fields to mosaic - e.g. 2132-50A."
    )

    parser.add_argument(
        "--datadirs",
        type=str,
        nargs='+',
        help="Directories containing cutouts (in subdir outdir/cutouts)..",
    )

    parser.add_argument(
        "--merge_name",
        type=str,
        help="Name of the merged region",
    )

    parser.add_argument(
        "--output_dir",
        type=str,
        help="Path to save merged data (in output_dir/merge_name)",
    )

    parser.add_argument(
        "--yanda",
        type=str,
        default="1.3.0",
        help="Yandasoft version to pull from DockerHub [1.3.0].",
    )

    parser.add_argument(
        "--host",
        type=str,
        help="Host of mongodb (probably $hostname -i).",
    )

    parser.add_argument(
        "-v", dest="verbose", action="store_true", help="Verbose output [False]."
    )

    parser.add_argument(
        "--username", type=str, default=None, help="Username of mongodb."
    )

    parser.add_argument(
        "--password", type=str, default=None, help="Password of mongodb."
    )

    args = parser.parse_args()

    cluster = LocalCluster(n_workers=1)
    client = Client(cluster)

    verbose = args.verbose
    test_db(
        host=args.host, username=args.username, password=args.password, verbose=verbose
    )

    main(
        fields=args.fields,
        field_dirs=args.datadirs,
        merge_name=args.merge_name,
        output_dir=args.output_dir,
        host=args.host,
        username=args.username,
        password=args.password,
        yanda=args.yanda,
        verbose=verbose,
    )


if __name__ == "__main__":
    cli()
