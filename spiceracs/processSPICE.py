#!/usr/bin/env python
import pymongo
from prefect import task, Task, Flow
from prefect.engine.executors import DaskExecutor
from spiceracs import cutout_rolling
from spiceracs import linmos
from spiceracs import rmsynth_oncuts
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
    rmsynth_task = task(rmsynth_oncuts.main, name='RM Synthesis')

    # Set up for Galaxy
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

    cluster.adapt(minimum=1, maximum=50)
    client = Client(cluster)

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
        mosaics = linmos_task(args.field,
                              args.datadir,
                              client,
                              host,
                              dryrun=False,
                              prefix="",
                              stokeslist=None,
                              verbose=True,
                              upstream_tasks=[cuts]
                              )
        dirty_spec = rmsynth_task(field=args.field,
                                  outdir=args.datadir,
                                  host=host,
                                  client=client,
                                  dimension=args.dimension,
                                  verbose=args.verbose,
                                  database=args.database,
                                  validate=args.validate,
                                  limit=args.limit,
                                  savePlots=args.savePlots,
                                  weightType=args.weightType,
                                  fitRMSF=args.fitRMSF,
                                  phiMax_radm2=args.phiMax_radm2,
                                  dPhi_radm2=args.dPhi_radm2,
                                  nSamples=args.nSamples,
                                  polyOrd=args.polyOrd,
                                  noStokesI=args.noStokesI,
                                  showPlots=args.showPlots,
                                  not_RMSF=args.not_RMSF,
                                  rm_verbose=args.verbose_worker,
                                  debug=args.debug,
                                  upstream_tasks=[mosaics]
                                  )

    # executor = DaskExecutor(address=client.scheduler.address)
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
        'field',
        metavar='field',
        type=str,
        help='Name of field (e.g. 2132-50A).')

    parser.add_argument(
        'datadir',
        metavar='datadir',
        type=str,
        help='Directory containing data cubes in FITS format.')

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
        '--islandcat',
        dest='islandcat',
        type=str,
        help='Master island RACS catalogue.')
    parser.add_argument(
        '-compcat',
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
        "--load",
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
        "--dryrun",
        dest="dryrun",
        action="store_true",
        help="Do a dry-run [False]."
    )

    parser.add_argument("--dimension", dest="dimension", default="1d",
                        help="How many dimensions for RMsynth [1d] or '3d'.")

    parser.add_argument(
        "-m",
        dest="database",
        action="store_true",
        help="Add RMsynth data to MongoDB [False]."
    )

    parser.add_argument("--validate", dest="validate", action="store_true",
                        help="Run on RMsynth Stokes I [False].")

    parser.add_argument("--limit", dest="limit", default=None,
                        type=int, help="Limit number of sources [All].")
    
    # RM-tools args
    parser.add_argument("-sp", dest="savePlots", action="store_true",
                        help="save the plots [False].")
    parser.add_argument("-w", dest="weightType", default="variance",
                        help="weighting [variance] (all 1s) or 'variance'.")
    parser.add_argument("-t", dest="fitRMSF", action="store_true",
                        help="Fit a Gaussian to the RMSF [False]")
    parser.add_argument("-l", dest="phiMax_radm2", type=float, default=None,
                        help="Absolute max Faraday depth sampled (overrides NSAMPLES) [Auto].")
    parser.add_argument("-d", dest="dPhi_radm2", type=float, default=None,
                        help="Width of Faraday depth channel [Auto].")
    parser.add_argument("-s", dest="nSamples", type=float, default=5,
                        help="Number of samples across the FWHM RMSF.")
    parser.add_argument("-o", dest="polyOrd", type=int, default=3,
                        help="polynomial order to fit to I spectrum [3].")
    parser.add_argument("-i", dest="noStokesI", action="store_true",
                        help="ignore the Stokes I spectrum [False].")
    parser.add_argument("--plots", dest="showPlots", action="store_true",
                        help="show the plots [False].")
    parser.add_argument("-R", dest="not_RMSF", action="store_true",
                        help="Skip calculation of RMSF? [False]")
    parser.add_argument("-rmv", dest="rm_verbose", action="store_true",
                        help="Verbose RMsynth [False].")
    parser.add_argument("-D", dest="debug", action="store_true",
                        help="turn on debugging messages & plots [False].")

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
