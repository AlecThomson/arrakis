#!/usr/bin/env python3
"""SPICE-RACS multi-field pipeline"""
import logging as log
import os
from time import sleep

import configargparse
import yaml
from astropy.time import Time
from dask import delayed, distributed
from dask.diagnostics import ProgressBar
from dask.distributed import Client, LocalCluster, performance_report, progress
from dask_jobqueue import SLURMCluster
from dask_mpi import initialize
from IPython import embed
from prefect import Flow, Task, task
from prefect.engine import signals
from prefect.engine.executors import DaskExecutor

from spiceracs import merge_fields, process_spice
from spiceracs.utils import port_forward, test_db


@task(name="Merge fields", skip_on_upstream_skip=False)
def merge_task(skip: bool, **kwargs) -> Task:
    """Cutout task

    Kwargs passed to merge_fields.main

    Args:
        skip (bool): Whether to skip this task

    Raises:
        signals.SKIP: If task is skipped

    Returns:
        Task: Runs merge_fields.main
    """
    if skip:
        raise signals.SKIP("Skipping merge task")
    return merge_fields.main(**kwargs)


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
        args.outfile = f"{args.merge_name}.pipe.test.fits"

    # Following https://github.com/dask/dask-jobqueue/issues/499
    with open(args.dask_config) as f:
        config = yaml.safe_load(f)

    config.update(
        {
            # 'scheduler_options': {
            # "dashboard_address": f":{args.port}"
            # },
            "log_directory": f"{args.merge_name}_{Time.now().fits}_spice_logs/"
        }
    )
    if args.use_mpi:
        initialize(
            interface=config["interface"],
            local_directory=config["local_directory"],
            nthreads=config["cores"] / config["processes"],
        )
        client = Client()
    else:
        cluster = SLURMCluster(**config,)
        log.debug(f"Submitted scripts will look like: \n {cluster.job_script()}")

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
    args_yaml_f = os.path.abspath(f"{args.merge_name}-config-{Time.now().fits}.yaml")
    log.info(f"Saving config to '{args_yaml_f}'")
    with open(args_yaml_f, "w") as f:
        f.write(args_yaml)

    port = client.scheduler_info()["services"]["dashboard"]

    # Forward ports
    if args.port_forward is not None:
        for p in args.port_forward:
            port_forward(port, p)

    # Prin out Dask client info
    log.info(client.scheduler_info()["services"])
    # Define flow
    inter_dir = os.path.join(os.path.abspath(args.output_dir), args.merge_name)
    with Flow(f"SPICE-RACS: {args.merge_name}") as flow:
        merge = merge_task(
            args.skip_merge,
            fields=args.fields,
            field_dirs=args.datadirs,
            merge_name=args.merge_name,
            output_dir=args.output_dir,
            client=client,
            host=host,
            username=args.username,
            password=args.password,
            yanda=args.yanda,
            verbose=args.verbose,
        )
        dirty_spec = process_spice.rmsynth_task(
            args.skip_rmsynth,
            field=args.merge_name,
            outdir=inter_dir,
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
            ion=False,
            do_own_fit=args.do_own_fit,
            upstream_tasks=[merge],
        )
        clean_spec = process_spice.rmclean_task(
            args.skip_rmclean,
            field=args.merge_name,
            outdir=inter_dir,
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
            window=args.window,
            showPlots=args.showPlots,
            rm_verbose=args.rm_verbose,
            upstream_tasks=[dirty_spec],
        )
        catalogue = process_spice.cat_task(
            args.skip_cat,
            field=args.merge_name,
            host=host,
            username=args.username,
            password=args.password,
            verbose=args.verbose,
            outfile=args.outfile,
            upstream_tasks=[clean_spec],
        )

    with performance_report(f"{args.merge_name}-report-{Time.now().fits}.html"):
        flow.run()

    client.close()


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
    SPICE-RACS regional pipeline.

    Before running make sure to start a session of mongodb e.g.
        $ mongod --dbpath=/path/to/database --bind_ip $(hostname -i)

    """

    # Parse the command line options
    parser = configargparse.ArgParser(
        default_config_files=[".default_field_config.txt"],
        description=descStr,
        formatter_class=configargparse.RawTextHelpFormatter,
    )
    parser.add("--config", required=False, is_config_file=True, help="Config file path")

    parser.add_argument(
        "--merge_name", type=str, help="Name of the merged region",
    )

    parser.add_argument(
        "--fields", type=str, nargs="+", help="RACS fields to mosaic - e.g. 2132-50A."
    )

    parser.add_argument(
        "--datadirs",
        type=str,
        nargs="+",
        help="Directories containing cutouts (in subdir outdir/cutouts)..",
    )

    parser.add_argument(
        "--output_dir",
        type=str,
        help="Path to save merged data (in output_dir/merge_name/cutouts)",
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

    parser.add_argument(
        "--yanda",
        type=str,
        default="1.3.0",
        help="Yandasoft version to pull from DockerHub [1.3.0].",
    )

    flowargs = parser.add_argument_group("pipeline flow options")
    flowargs.add_argument(
        "--skip_merge", action="store_true", help="Skip merge stage [False]."
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
        "--debugger", action="store_true", help="Debug output [False]."
    )
    options.add_argument(
        "-vw",
        "--verbose_worker",
        action="store_true",
        help="Verbose worker output [False].",
    )

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

    synth.add_argument(
        "--limit", default=None, type=int, help="Limit number of sources [All]."
    )
    synth.add_argument(
        "--own_fit",
        dest="do_own_fit",
        action="store_true",
        help="Use own Stokes I fit function [False].",
    )

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
    tools.add_argument(
        "--window",
        type=float,
        default=None,
        help="Further CLEAN in mask to this threshold [False].",
    )

    cat = parser.add_argument_group("catalogue arguments")
    # Cat args
    cat.add_argument(
        "--outfile", default=None, type=str, help="File to save table to [None]."
    )

    args = parser.parse_args()
    if not args.use_mpi:
        parser.print_values()

    verbose = args.verbose
    if verbose:
        log.basicConfig(
            level=log.INFO,
            format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            force=True,
        )
    if args.debugger:
        log.basicConfig(
            level=log.DEBUG,
            format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            force=True,
        )
    else:
        log.basicConfig(
            format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            force=True,
        )

    main(args)


if __name__ == "__main__":
    cli()
