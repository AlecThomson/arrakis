import numpy as np
import os
import pymongo

def rmsythoncut(args):
    i, collection, doc, verbose = args

    iname = doc['island_name']
    qfile = doc['q_file']
    ufile = doc['u_file']



    freqfile = '/avatar/athomson/cubes/8585/RACS_test4_1.05_1049-31A/cutouts_python/frequecies.txt'
    NSAMPLES = 10
    command = f"rmsynth3d {qfile} {ufile} {freqfile} -s {NSAMPLES} -o {iname}."
    if verbose: command += " -v"
    os.system(command)

    myquery = { "island_name": iname }
    newvalues = { "$set": { "rmsynth": True } }
    collection.update_one(myquery, newvalues)

def main(pool, verbose=False):

    client = pymongo.MongoClient()  # default connection (ie, local)
    mydb = client['racs']  # Create/open database
    mycol = mydb['spice']  # Create/open collection

    # Basic querey
    myquery = { "resolved": False }


    mydoc = mycol.find(myquery).sort("flux_int", -1)
    count = mydoc.count()

    if verbose: print(f'Running RMsynth on {count} sources')

    inputs = [[i, mycol, mydoc[i], verbose] for i in range(count)]
    list(pool.map(rmsythoncut, inputs))
    pool.close()

    if verbose: print('Done!')


if __name__ == "__main__":
    import argparse
    import schwimmbad

    # Help string to be shown using the -h option
    descStr = """
    Run RMsynthesis on all cubelets with unresolved sources.
    """

    # Parse the command line options
    parser = argparse.ArgumentParser(description=descStr,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-v", dest="verbose", action="store_true",
                        help="verbose output [False].")

    group = parser.add_mutually_exclusive_group()

    group.add_argument("--ncores", dest="n_cores", default=1,
                       type=int, help="Number of processes (uses multiprocessing).")
    group.add_argument("--mpi", dest="mpi", default=False,
                       action="store_true", help="Run with MPI.")


    args = parser.parse_args()
    pool = schwimmbad.choose_pool(mpi=args.mpi, processes=args.n_cores)

    verbose=args.verbose

    if args.mpi:
        if not pool.is_master():
            pool.wait()
            sys.exit(0)

    datadir = args.datadir[0]

    main(pool,
        verbose=verbose
        )