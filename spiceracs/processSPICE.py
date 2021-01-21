#!/usr/bin/env python
from spiceracs import cutout_rolling
#from spiceracs import linmos
from spiceracs import rmsynth_oncuts
from spiceracs import rmclean_oncuts
from spiceracs import makecat
import subprocess
import shlex
import schwimmbad
from dask_jobqueue import SLURMCluster
from distributed import Client, progress
from dask import delayed


def start_mongo(dbpath):
    proc = subprocess.run(shlex.split('hostname -i'),
                          stderr=subprocess.PIPE,
                          encoding='utf-8',
                          stdout=subprocess.PIPE
                          )
    host = proc.stdout.split()[0]

    command = f"numactl --interleave=all mongod --dbpath={dbpath} --bind_ip {host}"
    print(command)
    proc = subprocess.Popen(shlex.split(command),
                            stderr=subprocess.PIPE,
                            encoding='utf-8',
                            stdout=subprocess.PIPE
                            )
    return proc, host


def main(args, pool):
    proc, host = start_mongo(args.dbpath)
    cluster = SLURMCluster(cores=20,
                           memory="60GB",
                           project='askap',
                           queue='workq',
                           walltime='12:00:00',
                           job_extra=['-M galaxy'],
                           interface="ipogif0",
                           )
    cluster.scale(10)

    cutout_islands = delayed(cutout_rolling.cutout_islands)
    cutout_islands(args.field,
                   args.datadir,
                   pool,
                   host,
                   verbose=args.verbose,
                   pad=args.pad,
                   verbose_worker=args.verbose_worker,
                   dryrun=args.dryrun
                   )

    # cutout_rolling.cutout_islands(args.field,
    #                               args.datadir,
    #                               pool,
    #                               host,
    #                               verbose=verbose,
    #                               pad=args.pad,
    #                               verbose_worker=args.verbose_worker,
    #                               dryrun=args.dryrun
    #                               )


def cli():
    """Command-line interface
    """
    import argparse

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

    descStr = f"""
    {logostr}
    SPICE-RACS
    """

    # Parse the command line options
    parser = argparse.ArgumentParser(
        description=descStr, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'datadir',
        metavar='datadir',
        type=str,
        help='Directory containing data cubes in FITS format.')

    parser.add_argument(
        'field',
        metavar='field',
        type=str,
        help='Name of field (e.g. 2132-50A).')

    parser.add_argument(
        'dbpath',
        metavar='dbpath',
        type=str,
        help='Location of SPICE-RACS mongodb.')

    parser.add_argument(
        '-i',
        dest='islandcat',
        type=str,
        help='Master island RACS catalogue.')
    parser.add_argument(
        '-c',
        dest='compcat',
        type=str,
        help='Master component RACS catalogue.')
    parser.add_argument(
        "-v",
        dest="verbose",
        action="store_true",
        help="Verbose output [False]."
    )
    parser.add_argument(
        "-l",
        dest="load",
        action="store_true",
        help="Load catalogue into database [False]."
    )
    parser.add_argument(
        '-p',
        '--pad',
        dest='pad',
        type=float,
        default=3,
        help='Number of beamwidths to pad around source [3].')

    parser.add_argument(
        "-vw",
        dest="verbose_worker",
        action="store_true",
        help="Verbose worker output [False]."
    )
    parser.add_argument(
        "-d",
        dest="dryrun",
        action="store_true",
        help="Do a dry-run [False]."
    )

    group = parser.add_mutually_exclusive_group()

    group.add_argument(
        "--ncores",
        dest="n_cores",
        default=1,
        type=int, help="Number of processes (uses multiprocessing)."
    )
    group.add_argument(
        "--mpi",
        dest="mpi",
        default=False,
        action="store_true",
        help="Run with MPI."
    )

    pool = schwimmbad.choose_pool(mpi=args.mpi, processes=args.n_cores)
    if args.mpi:
        if not pool.is_master():
            pool.wait()
            sys.exit(0)
    # make it so we can use imap in serial and mpi mode
    if not isinstance(pool, schwimmbad.MultiPool):
        pool.imap = pool.map

    if verbose:
        print(f"Using pool: {pool.__class__.__name__}")

    main(args, pool)


if __name__ == "__main__":
    cli()
