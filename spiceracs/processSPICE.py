#!/usr/bin/env python3
"""SPICE-RACS pipeline script"""
import os
import pymongo
from prefect import task, Task, Flow
from prefect.engine.executors import DaskExecutor
from prefect.engine import signals
from spiceracs import cutout
from spiceracs import linmos
from spiceracs import cleanup
from spiceracs import rmsynth_oncuts
from spiceracs import rmclean_oncuts
from spiceracs import makecat
from spiceracs import frion
from spiceracs.utils import port_forward, test_db
from dask_jobqueue import SLURMCluster
from dask_mpi import initialize
from dask import distributed
from dask.distributed import Client, progress, performance_report, LocalCluster
from dask.diagnostics import ProgressBar
from dask import delayed
from IPython import embed
from time import sleep
from astropy.time import Time
import yaml
import socket
import configargparse


@task(name="Cutout")
def cut_task(skip: bool, **kwargs) -> Task:
    """Cutout task

    Kwargs passed to cutout.cutout_islands

    Args:
        skip (bool): Whether to skip this task

    Raises:
        signals.SUCCESS: If task is skipped

    Returns:
        Task: Runs cutout.cutout_islands
    """
    if skip:
        check_cond = True
    else:
        check_cond = False
    if check_cond:
        raise signals.SUCCESS
    else:
        return cutout.cutout_islands(**kwargs)


@task(name="LINMOS")
def linmos_task(skip: bool, **kwargs) -> Task:
    """LINOS task

    Kwargs passed to linmos.main

    Args:
        skip (bool): Whether to skip this task

    Raises:
        signals.SUCCESS: If task is skipped

    Returns:
        Task: Runs linmos.main
    """
    if skip:
        check_cond = True
    else:
        check_cond = False
    if check_cond:
        raise signals.SUCCESS
    else:
        return linmos.main(**kwargs)


@task(name="FRion")
def frion_task(skip: bool, **kwargs) -> Task:
    """FRion task

    Kwargs passed to frion.main

    Args:
        skip (bool): Whether to skip this task

    Raises:
        signals.SUCCESS: If task is skipped

    Returns:
        Task: Runs frion.main
    """
    if skip:
        check_cond = True
    else:
        check_cond = False
    if check_cond:
        raise signals.SUCCESS
    else:
        return frion.main(**kwargs)


@task(name="Clean up")
def cleanup_task(skip: bool, **kwargs) -> Task:
    """Cleanup task

    Kwargs passed to cleanup.main

    Args:
        skip (bool): Whether to skip this task

    Raises:
        signals.SUCCESS: If task is skipped

    Returns:
        Task: Runs cleanup.main
    """
    if skip:
        check_cond = True
    else:
        check_cond = False
    if check_cond:
        raise signals.SUCCESS
    else:
        return cleanup.main(**kwargs)


@task(name="RM Synthesis")
def rmsynth_task(skip: bool, **kwargs) -> Task:
    """RM synth task

    Kwargs passed to rmsynth_oncuts.main

    Args:
        skip (bool): Whether to skip this task

    Raises:
        signals.SUCCESS: If task is skipped

    Returns:
        Task: Runs rmsynth_oncuts.main
    """
    if skip:
        check_cond = True
    else:
        check_cond = False
    if check_cond:
        raise signals.SUCCESS
    else:
        return rmsynth_oncuts.main(**kwargs)


@task(name="RM-CLEAN")
def rmclean_task(skip: bool, **kwargs) -> Task:
    """RM-CLEAN task

    Kwargs passed to rmclean_oncuts.main

    Args:
        skip (bool): Whether to skip this task

    Raises:
        signals.SUCCESS: If task is skipped

    Returns:
        Task: Runs rmclean_oncuts.main
    """
    if skip:
        check_cond = True
    else:
        check_cond = False
    if check_cond:
        raise signals.SUCCESS
    else:
        return rmclean_oncuts.main(**kwargs)


@task(name="Catalogue")
def cat_task(skip: bool, **kwargs) -> Task:
    """Catalogue task

    Kwargs passed to makecat.main

    Args:
        skip (bool): Whether to skip this task

    Raises:
        signals.SUCCESS: If task is skipped

    Returns:
        Task: Runs makecat.main
    """
    if skip:
        check_cond = True
    else:
        check_cond = False
    if check_cond:
        raise signals.SUCCESS
    else:
        return makecat.main(**kwargs)


def main(args: configargparse.Namespace) -> None:
    """Main script

    Args:
        args (configargparse.Namespace): Command line arguments.
    """
    host = args.host

    if args.dask_config is None:
        scriptdir = os.path.dirname(os.path.realpath(__file__))
        config_dir = f"{scriptdir}/../configs"
        args.dask_config = f"{config_dir}/default.yaml"

    if args.outfile is None:
        args.outfile = f"{args.field}.pipe.test.fits"

    # Following https://github.com/dask/dask-jobqueue/issues/499
    with open(args.dask_config) as f:
        config = yaml.safe_load(f)

    config.update(
        {
            # 'scheduler_options': {
            # "dashboard_address": f":{args.port}"
            # },
            "log_directory": f"{args.field}_{Time.now().fits}_spice_logs/"
        }
    )
    if args.use_mpi:
        initialize(
            interface=config["interface"],
            local_directory=config["local_directory"],
            # dashboard_address=f":{args.port}",
        )
        client = Client()
    else:
        cluster = SLURMCluster(
            **config,
        )
        print("Submitted scripts will look like: \n", cluster.job_script())

        # Request 15 nodes
        cluster.scale(jobs=15)
        # cluster = LocalCluster(n_workers=10, processes=True, threads_per_worker=1, local_directory="/dev/shm",dashboard_address=f":{args.port}")
        client = Client(cluster)

    test_db(
        host=args.host,
        username=args.username,
        password=args.password,
        verbose=args.verbose,
    )

    args_yaml = yaml.dump(vars(args))
    args_yaml_f = os.path.abspath(f"{args.field}-config-{Time.now().fits}.yaml")
    if args.verbose:
        print(f"Saving config to '{args_yaml_f}'")
    with open(args_yaml_f, "w") as f:
        f.write(args_yaml)

    port = client.scheduler_info()["services"]["dashboard"]

    # Forward ports
    if args.port_forward is not None:
        for p in args.port_forward:
            port_forward(port, p)

    # Prin out Dask client info
    print(client.scheduler_info()["services"])

    # Define flow
    with Flow(f"SPICE-RACS: {args.field}") as flow:
        cuts = cut_task(
            args.skip_cutout,
            field=args.field,
            directory=args.datadir,
            host=host,
            username=args.username,
            password=args.password,
            client=client,
            verbose=args.verbose,
            pad=args.pad,
            stokeslist=["I", "Q", "U"],
            verbose_worker=args.verbose_worker,
            dryrun=args.dryrun,
        )
        mosaics = linmos_task(
            args.skip_linmos,
            field=args.field,
            datadir=args.datadir,
            client=client,
            host=host,
            holofile=args.holofile,
            username=args.username,
            password=args.password,
            yanda=args.yanda,
            dryrun=False,
            prefix="",
            stokeslist=["I", "Q", "U"],
            verbose=True,
            upstream_tasks=[cuts],
        )
        tidy = cleanup_task(
            args.skip_cleanup,
            datadir=args.datadir,
            client=client,
            stokeslist=["I", "Q", "U"],
            verbose=True,
            upstream_tasks=[mosaics],
        )
        frion_run = frion_task(
            args.skip_frion,
            field=args.field,
            outdir=args.datadir,
            host=host,
            client=client,
            username=args.username,
            password=args.password,
            database=args.database,
            verbose=args.verbose,
            upstream_tasks=[mosaics],
        )
        dirty_spec = rmsynth_task(
            args.skip_rmsynth,
            field=args.field,
            outdir=args.datadir,
            host=host,
            username=args.username,
            password=args.password,
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
            fit_function=args.fit_function,
            tt0=args.tt0,
            tt1=args.tt1,
            ion=True,
            upstream_tasks=[frion_run],
        )
        clean_spec = rmclean_task(
            args.skip_rmclean,
            field=args.field,
            outdir=args.datadir,
            host=host,
            username=args.username,
            password=args.password,
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
            upstream_tasks=[dirty_spec],
        )
        catalogue = cat_task(
            args.skip_cat,
            field=args.field,
            host=host,
            username=args.username,
            password=args.password,
            verbose=args.verbose,
            outfile=args.outfile,
            cat_format=args.format,
            upstream_tasks=[clean_spec],
        )

    with performance_report(f"{args.field}-report-{Time.now().fits}.html"):
        flow.run()


def cli():
    """Command-line interface"""
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
    SPICE-RACS pipeline.

    Before running make sure to start a session of mongodb e.g.
        $ mongod --dbpath=/path/to/database --bind_ip $(hostname -i)

    """

    # Parse the command line options
    parser = configargparse.ArgParser(
        default_config_files=[".default_config.txt"],
        description=descStr,
        formatter_class=configargparse.RawTextHelpFormatter,
    )
    parser.add("--config", required=False, is_config_file=True, help="Config file path")
    parser.add_argument(
        "field", metavar="field", type=str, help="Name of field (e.g. 2132-50A)."
    )

    parser.add_argument(
        "datadir",
        metavar="datadir",
        type=str,
        help="Directory containing data cubes in FITS format.",
    )

    parser.add_argument(
        "--host",
        default=None,
        type=str,
        help="Host of mongodb (probably $hostname -i).",
    )

    parser.add_argument(
        "--username", type=str, default=None, help="Username of mongodb."
    )

    parser.add_argument(
        "--password", type=str, default=None, help="Password of mongodb."
    )

    # parser.add_argument(
    #     '--port',
    #     type=int,
    #     default=9999,
    #     help="Port to run Dask dashboard on."
    # )
    parser.add_argument(
        "--use_mpi",
        action="store_true",
        help="Use Dask-mpi to parallelise -- must use srun/mpirun to assign resources.",
    )
    parser.add_argument(
        "--port_forward",
        default=None,
        help="Platform to fowards dask port [None].",
        nargs="+",
    )

    parser.add_argument(
        "--dask_config",
        type=str,
        default=None,
        help="Config file for Dask SlurmCLUSTER.",
    )
    parser.add_argument("--holofile", type=str, help="Path to holography image")

    parser.add_argument(
        "--yanda",
        type=str,
        default="1.3.0",
        help="Yandasoft version to pull from DockerHub [1.3.0].",
    )

    flowargs = parser.add_argument_group("pipeline flow options")
    flowargs.add_argument(
        "--skip_cutout", action="store_true", help="Skip cutout stage [False]."
    )
    flowargs.add_argument(
        "--skip_linmos", action="store_true", help="Skip LINMOS stage [False]."
    )
    flowargs.add_argument(
        "--skip_cleanup", action="store_true", help="Skip cleanup stage [False]."
    )
    flowargs.add_argument(
        "--skip_frion", action="store_true", help="Skip cleanup stage [False]."
    )
    flowargs.add_argument(
        "--skip_rmsynth", action="store_true", help="Skip RM Synthesis stage [False]."
    )
    flowargs.add_argument(
        "--skip_rmclean", action="store_true", help="Skip RM-CLEAN stage [False]."
    )
    flowargs.add_argument(
        "--skip_cat", action="store_true", help="Skip catalogue stage [False]."
    )

    options = parser.add_argument_group("output options")
    options.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose output [False]."
    )
    options.add_argument(
        "-vw",
        "--verbose_worker",
        action="store_true",
        help="Verbose worker output [False].",
    )
    cutargs = parser.add_argument_group("cutout arguments")
    cutargs.add_argument(
        "-p",
        "--pad",
        type=float,
        default=5,
        help="Number of beamwidths to pad around source [5].",
    )

    cutargs.add_argument("--dryrun", action="store_true", help="Do a dry-run [False].")

    synth = parser.add_argument_group("RM-synth/CLEAN arguments")

    synth.add_argument(
        "--dimension",
        default="1d",
        help="How many dimensions for RMsynth [1d] or '3d'.",
    )

    synth.add_argument(
        "-m",
        "--database",
        action="store_true",
        help="Add RMsynth data to MongoDB [False].",
    )

    synth.add_argument(
        "--tt0",
        default=None,
        type=str,
        help="TT0 MFS image -- will be used for model of Stokes I -- also needs --tt1.",
    )

    synth.add_argument(
        "--tt1",
        default=None,
        type=str,
        help="TT1 MFS image -- will be used for model of Stokes I -- also needs --tt0.",
    )

    synth.add_argument(
        "--validate", action="store_true", help="Run on RMsynth Stokes I [False]."
    )

    synth.add_argument("--limit", default=None, help="Limit number of sources [All].")
    tools = parser.add_argument_group("RM-tools arguments")
    # RM-tools args
    tools.add_argument(
        "-sp", "--savePlots", action="store_true", help="save the plots [False]."
    )
    tools.add_argument(
        "-w",
        "--weightType",
        default="variance",
        help="weighting [variance] (all 1s) or 'uniform'.",
    )
    tools.add_argument(
        "--fit_function",
        type=str,
        default="log",
        help="Stokes I fitting function: 'linear' or ['log'] polynomials.",
    )
    tools.add_argument(
        "-t",
        "--fitRMSF",
        action="store_true",
        help="Fit a Gaussian to the RMSF [False]",
    )
    tools.add_argument(
        "-l",
        "--phiMax_radm2",
        type=float,
        default=None,
        help="Absolute max Faraday depth sampled (overrides NSAMPLES) [Auto].",
    )
    tools.add_argument(
        "-d",
        "--dPhi_radm2",
        type=float,
        default=None,
        help="Width of Faraday depth channel [Auto].",
    )
    tools.add_argument(
        "-s",
        "--nSamples",
        type=float,
        default=5,
        help="Number of samples across the FWHM RMSF.",
    )
    tools.add_argument(
        "-o",
        "--polyOrd",
        type=int,
        default=3,
        help="polynomial order to fit to I spectrum [3].",
    )
    tools.add_argument(
        "-i",
        "--noStokesI",
        action="store_true",
        help="ignore the Stokes I spectrum [False].",
    )
    tools.add_argument(
        "--showPlots", action="store_true", help="show the plots [False]."
    )
    tools.add_argument(
        "-R",
        "--not_RMSF",
        action="store_true",
        help="Skip calculation of RMSF? [False]",
    )
    tools.add_argument(
        "-rmv",
        "--rm_verbose",
        action="store_true",
        help="Verbose RMsynth/CLEAN [False].",
    )
    tools.add_argument(
        "-D",
        "--debug",
        action="store_true",
        help="turn on debugging messages & plots [False].",
    )
    # RM-tools args
    tools.add_argument(
        "-c",
        "--cutoff",
        type=float,
        default=-3,
        help="CLEAN cutoff (+ve = absolute, -ve = sigma) [-3].",
    )
    tools.add_argument(
        "-n",
        "--maxIter",
        type=int,
        default=10000,
        help="maximum number of CLEAN iterations [10000].",
    )
    tools.add_argument(
        "-g", "--gain", type=float, default=0.1, help="CLEAN loop gain [0.1]."
    )

    cat = parser.add_argument_group("catalogue arguments")
    # Cat args
    cat.add_argument(
        "--outfile", default=None, type=str, help="File to save table to [None]."
    )

    cat.add_argument(
        "-f", "--format", default=None, type=str, help="Format for output file [None]."
    )
    args = parser.parse_args()
    if not args.use_mpi:
        parser.print_values()

    main(args)


if __name__ == "__main__":
    cli()
