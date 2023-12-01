#!/usr/bin/env python3
"""Arrakis multi-field pipeline"""
import os

import configargparse
import pkg_resources
import yaml
from astropy.time import Time
from dask.distributed import Client
from dask_jobqueue import SLURMCluster
from dask_mpi import initialize
from prefect import flow, task
from prefect_dask import DaskTaskRunner

from arrakis import merge_fields, process_spice
from arrakis.logger import logger
from arrakis.utils.database import test_db
from arrakis.utils.pipeline import logo_str, port_forward

merge_task = task(merge_fields.main, name="Merge fields")


@flow
def process_merge(args, host: str, inter_dir: str) -> None:
    """Workflow to merge spectra from overlapping fields together

    Args:
        args (Namespace): Parameters to use for this process
        host (str): Address of the mongoDB servicing the processing
        inter_dir (str): Location to store data from merged fields
    """
    previous_future = None
    previous_future = (
        merge_task.submit(
            fields=args.fields,
            field_dirs=args.datadirs,
            merge_name=args.merge_name,
            output_dir=args.output_dir,
            host=host,
            epoch=args.epoch,
            username=args.username,
            password=args.password,
            yanda=args.yanda,
            verbose=args.verbose,
        )
        if not args.skip_merge
        else previous_future
    )

    previous_future = (
        process_spice.rmsynth_task.submit(
            field=args.merge_name,
            outdir=inter_dir,
            host=host,
            epoch=args.epoch,
            username=args.username,
            password=args.password,
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
            wait_for=[previous_future],
        )
        if not args.skip_rmsynth
        else previous_future
    )

    previous_future = (
        process_spice.rmclean_task.submit(
            field=args.merge_name,
            outdir=inter_dir,
            host=host,
            epoch=args.epoch,
            username=args.username,
            password=args.password,
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
            wait_for=[previous_future],
        )
        if not args.skip_rmclean
        else previous_future
    )

    previous_future = (
        process_spice.cat_task.submit(
            field=args.merge_name,
            host=host,
            epoch=args.epoch,
            username=args.username,
            password=args.password,
            verbose=args.verbose,
            outfile=args.outfile,
            wait_for=[previous_future],
        )
        if not args.skip_cat
        else previous_future
    )


def main(args: configargparse.Namespace) -> None:
    """Main script

    Args:
        args (configargparse.Namespace): Command line arguments.
    """
    host = args.host

    if args.dask_config is None:
        config_dir = pkg_resources.resource_filename("arrakis", "configs")
        args.dask_config = os.path.join(config_dir, "default.yaml")

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
        cluster = SLURMCluster(
            **config,
        )
        logger.debug(f"Submitted scripts will look like: \n {cluster.job_script()}")

        client = Client(cluster)

    test_db(
        host=args.host,
        username=args.username,
        password=args.password,
    )

    args_yaml = yaml.dump(vars(args))
    args_yaml_f = os.path.abspath(f"{args.merge_name}-config-{Time.now().fits}.yaml")
    logger.info(f"Saving config to '{args_yaml_f}'")
    with open(args_yaml_f, "w") as f:
        f.write(args_yaml)

    port = client.scheduler_info()["services"]["dashboard"]

    # Forward ports
    if args.port_forward is not None:
        for p in args.port_forward:
            port_forward(port, p)

    # Prin out Dask client info
    logger.info(client.scheduler_info()["services"])

    dask_runner = DaskTaskRunner(address=client.scheduler.address)

    inter_dir = os.path.join(os.path.abspath(args.output_dir), args.merge_name)

    process_merge.with_options(
        name=f"SPICE-RACS: {args.merge_name}", task_runner=dask_runner
    )(args, host, inter_dir)

    client.close()


def cli():
    """Command-line interface"""
    # Help string to be shown using the -h option
    descStr = f"""
    {logo_str}
    Arrakis regional pipeline.

    Before running make sure to start a session of mongodb e.g.
        $ mongod --dbpath=/path/to/database --bind_ip $(hostname -i)

    """

    # Parse the command line options
    parser = configargparse.ArgParser(
        default_config_files=[".default_field_config.txt"],
        description=descStr,
        formatter_class=configargparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add("--config", required=False, is_config_file=True, help="Config file path")

    parser.add_argument(
        "--merge_name",
        type=str,
        help="Name of the merged region",
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
        "--epoch",
        type=int,
        default=0,
        help="Epoch to read field data from",
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
        logger.setLevel(logger.INFO)
    if args.debugger:
        logger.setLevel(logger.DEBUG)

    main(args)


if __name__ == "__main__":
    cli()
