#!/usr/bin/env python
from prefect import task, Task, Flow
from prefect.engine.executors import DaskExecutor
from spiceracs import cutout_rolling
from spiceracs import linmos
import pymongo
# from spiceracs import rmsynth_oncuts
# from spiceracs import rmclean_oncuts
# from spiceracs import makecat
import subprocess
import shlex
from dask_jobqueue import SLURMCluster
from distributed import Client, progress, performance_report
from dask.diagnostics import ProgressBar
from dask import delayed
from IPython import embed
from time import sleep


def main(args):
    # proc, host = start_mongo(args.dbpath)
    host = args.host
    cut_task = task(cutout_rolling.cutout_islands, name='cutout')
    linmos_task = task(linmos.main, name='LINMOS')
    cluster = SLURMCluster(cores=20,
                           memory="60GB",
                           project='askap',
                           queue='workq',
                           walltime='12:00:00',
                           job_extra=['-M galaxy'],
                           #    interface="eth2"
                           interface="ipogif0",
                           log_directory='logs',
                           env_extra=['module load askapsoft']
                           )
    # cluster.scale(nworkers)
    cluster.adapt(minimum=1, maximum=50)
    client = Client(cluster)
    # while ((client.status == "running") and (len(client.scheduler_info()["workers"]) < nworkers)):
    #     sleep(1.0)
    print(client.scheduler_info()['services'])
    with Flow(f'SPICE-RACS: {args.field}') as flow:
        cuts = cut_task(args.field,
                        args.datadir,
                        host,
                        client,
                        verbose=args.verbose,
                        pad=args.pad,
                        verbose_worker=args.verbose_worker,
                        dryrun=args.dryrun
                        )
        mosaic = linmos_task(args.field,
                             args.datadir,
                             client,
                             host,
                             dryrun=False,
                             prefix="",
                             stokeslist=None,
                             verbose=True,
                             upstream_tasks=[cuts]
                             )
    # executor = DaskExecutor(address=client.scheduler.address)
    embed()
    with performance_report(f'{args.field}-report.html'):
        flow.run()



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
        'host',
        metavar='host',
        type=str,
        help='Host of mongodb (probably $hostname -i).')

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

    args = parser.parse_args()

    verbose = args.verbose

    if verbose:
        print('Testing MongoDB connection...')
    # default connection (ie, local)
    with pymongo.MongoClient(host=args.host, connect=False) as client:
        try:
            client.list_database_names()
        except pymongo.errors.ServerSelectionTimeoutError:
            raise Exception("Please ensure 'mongod' is running")
        else:
            if verbose:
                print('MongoDB connection succesful!')

    main(args)


if __name__ == "__main__":
    cli()
