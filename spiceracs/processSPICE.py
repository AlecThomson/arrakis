#!/usr/bin/env python3
import pymongo
from prefect import task, Task, Flow
from prefect.engine.executors import DaskExecutor
from spiceracs import cutout
from spiceracs import linmos
from spiceracs import cleanup
from spiceracs import rmsynth_oncuts
from spiceracs import rmclean_oncuts
from spiceracs import makecat
from dask_jobqueue import SLURMCluster
from distributed import Client, progress, performance_report, LocalCluster
from dask.diagnostics import ProgressBar
from dask import delayed
from IPython import embed
from time import sleep


def main(args):
    # proc, host = start_mongo(args.dbpath)
    host = args.host
    cut_task = task(cutout.cutout_islands, name='cutout')
    linmos_task = task(linmos.main, name='LINMOS')
    cleanup_task = task(cleanup.main, name='Clean up')
    rmsynth_task = task(rmsynth_oncuts.main, name='RM Synthesis')
    rmclean_task = task(rmclean_oncuts.main, name='RM-CLEAN')
    cat_task = task(makecat.main, name='Catalogue')

    # Set up for Galaxy
    cluster = SLURMCluster(cores=20,
                           memory="60GB",
                           project='askap',
                           queue='workq',
                           walltime='12:00:00',
                           job_extra=['-M galaxy'],
                           interface="ipogif0",
                           log_directory='logs',
                           env_extra=[
                               'module load askapsoft'
                               'unset PYTHONPATH',
                               'source /home/$(whoami)/.bashrc',
                               'conda activate base'
                           ],
                           dashboard_address=":9999",
                           python='srun -n 1 -c 20 python',
                           )
    # Set up for Zeus
    # cluster = SLURMCluster(cores=28,
    #                     memory="120GB",
    #                     project='askaprt',
    #                     queue='workq',
    #                     walltime='12:00:00',
    #                     job_extra=['-M zeus'],
    #                     interface='ib0',
    #                     log_directory='logs',
    #                     env_extra=['module load askapsoft'],
    #                     dashboard_address=":9999"
    #                     )

    # Request up to 20 nodes
    cluster.adapt(minimum=0, maximum=20)
    # cluster.scale(20)
    # cluster = LocalCluster(n_workers=20)
    client = Client(cluster)

    # Prin out Dask client info
    print(client.scheduler_info()['services'])

    # Define flow
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
                              stokeslist=["I", "Q", "U"],
                              verbose=True,
                              upstream_tasks=[cuts]
                              )
        tidy = cleanup_task(datadir=args.datadir,
                            client=client,
                            stokeslist=["I", "Q", "U"],
                            verbose=True,
                            upstream_tasks=[mosaics]
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
                                  rm_verbose=args.rm_verbose,
                                  debug=args.debug,
                                  upstream_tasks=[tidy]
                                  # upstream_tasks=[mosaics]
                                  )
        clean_spec = rmclean_task(field=args.field,
                                  outdir=args.datadir,
                                  host=host,
                                  client=client,
                                  dimension=args.dimension,
                                  verbose=args.verbose,
                                  database=args.database,
                                  validate=args.validate,
                                  limit=args.limit,
                                  cutoff=args.cutoff,
                                  maxIter=args.maxIter,
                                  gain=args.gain,
                                  showPlots=args.showPlots,
                                  rm_verbose=args.rm_verbose,
                                  upstream_tasks=[dirty_spec]
                                  )
        catalogue = cat_task(args.field,
                             host,
                             verbose=args.verbose,
                             limit=args.limit,
                             outfile=f'{args.field}.pipe.test.fits',
                             cat_format='fits',
                             upstream_tasks=[clean_spec]
                             )

    with performance_report(f'{args.field}-report.html'):
        flow.run()


def cli():
    """Command-line interface
    """
    import configargparse

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
    parser = configargparse.ArgParser(default_config_files=['.default_config.txt'],
                                      description=descStr,
                                      formatter_class=configargparse.RawTextHelpFormatter
                                      )
    parser.add('--config', required=False,
               is_config_file=True, help='Config file path')
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
    options = parser.add_argument_group("output options")
    options.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Verbose output [True]."
    )
    options.add_argument(
        "-vw",
        "--verbose_worker",
        action="store_true",
        help="Verbose worker output [False]."
    )
    cutout = parser.add_argument_group("cutout arguments")
    cutout.add_argument(
        '-p',
        '--pad',
        type=float,
        default=5,
        help='Number of beamwidths to pad around source [5].')
    synth = parser.add_argument_group("RM-synth/CLEAN arguments")
    cutout.add_argument(
        "--dryrun",
        action="store_true",
        help="Do a dry-run [False]."
    )

    synth.add_argument("--dimension", default="1d",
                       help="How many dimensions for RMsynth [1d] or '3d'.")

    synth.add_argument(
        "-m",
        "--database",
        action="store_true",
        help="Add RMsynth data to MongoDB [False]."
    )

    synth.add_argument("--validate",
                       action="store_true",
                       help="Run on RMsynth Stokes I [False].")

    synth.add_argument("--limit",
                       default=None,
                       help="Limit number of sources [All].")
    tools = parser.add_argument_group("RM-tools arguments")
    # RM-tools args
    tools.add_argument("-sp",
                       "--savePlots",
                       action="store_true",
                       help="save the plots [False].")
    tools.add_argument("-w",
                       "--weightType",
                       default="variance",
                       help="weighting [variance] (all 1s) or 'variance'.")
    tools.add_argument("-t", "--fitRMSF", action="store_true",
                       help="Fit a Gaussian to the RMSF [False]")
    tools.add_argument("-l", "--phiMax_radm2", type=float, default=None,
                       help="Absolute max Faraday depth sampled (overrides NSAMPLES) [Auto].")
    tools.add_argument("-d", "--dPhi_radm2", type=float, default=None,
                       help="Width of Faraday depth channel [Auto].")
    tools.add_argument("-s", "--nSamples", type=float, default=5,
                       help="Number of samples across the FWHM RMSF.")
    tools.add_argument("-o", "--polyOrd", type=int, default=3,
                       help="polynomial order to fit to I spectrum [3].")
    tools.add_argument("-i", "--noStokesI", action="store_true",
                       help="ignore the Stokes I spectrum [False].")
    tools.add_argument("--showPlots", action="store_true",
                       help="show the plots [False].")
    tools.add_argument("-R", "--not_RMSF", action="store_true",
                       help="Skip calculation of RMSF? [False]")
    tools.add_argument("-rmv", "--rm_verbose", action="store_true",
                       help="Verbose RMsynth/CLEAN [False].")
    tools.add_argument("-D", "--debug", action="store_true",
                       help="turn on debugging messages & plots [False].")
    # RM-tools args
    tools.add_argument("-c", "--cutoff", type=float, default=-3,
                       help="CLEAN cutoff (+ve = absolute, -ve = sigma) [-3].")
    tools.add_argument("-n", "--maxIter", type=int, default=10000,
                       help="maximum number of CLEAN iterations [10000].")
    tools.add_argument("-g", "--gain", type=float, default=0.1,
                       help="CLEAN loop gain [0.1].")

    cat = parser.add_argument_group("catalogue arguments")
    # Cat args
    cat.add_argument("--outfile", default=None,
                     type=str, help="File to save table to [None].")

    cat.add_argument("-f", "--format", default=None,
                     type=str, help="Format for output file [None].")
    args = parser.parse_args()

    verbose = args.verbose

    if verbose:
        print('Testing MongoDB connection...')
    # default connection (ie, local)
    with pymongo.MongoClient(host=args.host, connect=False) as dbclient:
        try:
            dbclient.list_database_names()
        except pymongo.errors.ServerSelectionTimeoutError:
            raise Exception("Please ensure 'mongod' is running")
        else:
            if verbose:
                print('MongoDB connection succesful!')

    main(args)


if __name__ == "__main__":
    cli()
