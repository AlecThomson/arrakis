#!/usr/bin/env python3
"""Arrakis single-field pipeline"""
import logging
import os
from pathlib import Path

import configargparse
import pkg_resources
import yaml
from astropy.time import Time
from prefect import flow
from prefect.task_runners import BaseTaskRunner
from prefect_dask import DaskTaskRunner

from arrakis import (
    cleanup,
    cutout,
    frion,
    imager,
    linmos,
    makecat,
    rmclean_oncuts,
    rmsynth_oncuts,
)
from arrakis.logger import logger
from arrakis.utils.database import test_db
from arrakis.utils.pipeline import logo_str


@flow(name="Combining+Synthesis on Arrakis")
def process_spice(args, host: str, task_runner: BaseTaskRunner) -> None:
    """Workflow to process the SPIRCE-RACS data

    Args:
        args (configargparse.Namespace): Configuration parameters for this run
        host (str): Host address of the mongoDB.
    """
    outfile = f"{args.field}.pipe.test.fits" if args.outfile is None else args.outfile

    previous_future = None
    previous_future = (
        cutout.cutout_islands.with_options(
            task_runner=task_runner,
        )(
            field=args.field,
            directory=str(args.outdir),
            host=host,
            epoch=args.epoch,
            username=args.username,
            password=args.password,
            pad=args.pad,
            stokeslist=["I", "Q", "U"],
            dryrun=args.dryrun,
            limit=args.limit,
        )
        if not args.skip_cutout
        else previous_future
    )

    previous_future = (
        linmos.main.with_options(
            task_runner=task_runner,
        )(
            field=args.field,
            datadir=Path(args.outdir),
            host=host,
            epoch=args.epoch,
            holofile=Path(args.holofile),
            username=args.username,
            password=args.password,
            yanda=args.yanda,
            yanda_img=args.yanda_image,
            stokeslist=["I", "Q", "U"],
            limit=args.limit,
        )
        if not args.skip_linmos
        else previous_future
    )

    previous_future = (
        frion.main.with_options(task_runner=task_runner)(
            field=args.field,
            outdir=args.outdir,
            host=host,
            epoch=args.epoch,
            username=args.username,
            password=args.password,
            database=args.database,
            ionex_server=args.ionex_server,
            ionex_prefix=args.ionex_prefix,
            ionex_proxy_server=args.ionex_proxy_server,
            ionex_formatter=args.ionex_formatter,
            ionex_predownload=args.ionex_predownload,
            limit=args.limit,
        )
        if not args.skip_frion
        else previous_future
    )

    previous_future = (
        rmsynth_oncuts.main.with_options(task_runner=task_runner)(
            field=args.field,
            outdir=args.outdir,
            host=host,
            epoch=args.epoch,
            username=args.username,
            password=args.password,
            dimension=args.dimension,
            verbose=args.verbose,
            database=args.database,
            do_validate=args.validate,
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
            do_own_fit=args.do_own_fit,
        )
        if not args.skip_rmsynth
        else previous_future
    )

    previous_future = (
        rmclean_oncuts.main.with_options(task_runner=task_runner)(
            field=args.field,
            outdir=args.outdir,
            host=host,
            epoch=args.epoch,
            username=args.username,
            password=args.password,
            dimension=args.dimension,
            database=args.database,
            limit=args.limit,
            cutoff=args.cutoff,
            maxIter=args.maxIter,
            gain=args.gain,
            window=args.window,
            showPlots=args.showPlots,
            rm_verbose=args.rm_verbose,
        )
        if not args.skip_rmclean
        else previous_future
    )

    previous_future = (
        makecat.main.with_options(task_runner=task_runner)(
            field=args.field,
            host=host,
            epoch=args.epoch,
            username=args.username,
            password=args.password,
            verbose=args.verbose,
            outfile=outfile,
        )
        if not args.skip_cat
        else previous_future
    )

    previous_future = (
        cleanup.main.with_options(task_runner=task_runner)(
            datadir=args.outdir,
        )
        if not args.skip_cleanup
        else previous_future
    )


def save_args(args: configargparse.Namespace) -> Path:
    """Helper function to create a record of the input configuration arguments that
    govern the pipeline instance

    Args:
        args (configargparse.Namespace): Supplied arguments for the Arrakis pipeline instance

    Returns:
        Path: Output path of the saved file
    """
    args_yaml = yaml.dump(vars(args))
    args_yaml_f = os.path.abspath(f"{args.field}-config-{Time.now().fits}.yaml")
    logger.info(f"Saving config to '{args_yaml_f}'")
    with open(args_yaml_f, "w") as f:
        f.write(args_yaml)

    return Path(args_yaml_f)


def create_dask_runner(
    dask_config: str,
    overload: bool = False,
) -> DaskTaskRunner:
    """Create a DaskTaskRunner

    Args:
        dask_config (str): Configuraiton file for the DaskTaskRunner
        overload (bool, optional): Overload the options for threadded work. Defaults to False.

    Returns:
        DaskTaskRunner: The prefect DaskTaskRunner instance
    """
    logger.setLevel(logging.INFO)
    logger.info("Creating a Dask Task Runner.")
    if dask_config is None:
        config_dir = pkg_resources.resource_filename("arrakis", "configs")
        dask_config = f"{config_dir}/default.yaml"

    with open(dask_config) as f:
        logger.info(f"Loading {dask_config}")
        yaml_config: dict = yaml.safe_load(f)

    cluster_class_str = yaml_config.get("cluster_class", "distributed.LocalCluster")
    cluster_kwargs = yaml_config.get("cluster_kwargs", {})
    adapt_kwargs = yaml_config.get("adapt_kwargs", {})

    if overload:
        logger.info("Overwriting config attributes.")
        cluster_kwargs["job_cpu"] = cluster_kwargs["cores"]
        cluster_kwargs["cores"] = 1
        cluster_kwargs["processes"] = 1

    config = {
        "cluster_class": cluster_class_str,
        "cluster_kwargs": cluster_kwargs,
        "adapt_kwargs": adapt_kwargs,
    }

    return DaskTaskRunner(**config)


def main(args: configargparse.Namespace) -> None:
    """Main script

    Args:
        args (configargparse.Namespace): Command line arguments.
    """
    host = args.host

    # Lets save the args as a record for the ages
    output_args_path = save_args(args)
    logger.info(f"Saved arguments to {output_args_path}.")

    # Test the mongoDB
    test_db(
        host=host,
        username=args.username,
        password=args.password,
    )

    if not args.skip_imager:
        # This is the client for the imager component of the arrakis
        # pipeline.
        dask_runner = create_dask_runner(
            dask_config=args.imager_dask_config,
            overload=True,
        )

        logger.info("Obtained DaskTaskRunner, executing the imager workflow. ")
        imager.main.with_options(
            name=f"Arrakis Imaging -- {args.field}", task_runner=dask_runner
        )(
            msdir=args.msdir,
            out_dir=args.outdir,
            temp_dir=args.temp_dir,
            cutoff=args.psf_cutoff,
            robust=args.robust,
            pols=args.pols,
            nchan=args.nchan,
            local_rms=args.local_rms,
            local_rms_window=args.local_rms_window,
            size=args.size,
            scale=args.scale,
            mgain=args.mgain,
            niter=args.niter,
            nmiter=args.nmiter,
            auto_mask=args.auto_mask,
            force_mask_rounds=args.force_mask_rounds,
            auto_threshold=args.auto_threshold,
            minuv=args.minuv,
            purge=args.purge,
            taper=args.taper,
            parallel_deconvolution=args.parallel,
            gridder=args.gridder,
            wsclean_path=(
                Path(args.local_wsclean) if args.local_wsclean else args.hosted_wsclean
            ),
            multiscale=args.multiscale,
            multiscale_scale_bias=args.multiscale_scale_bias,
            absmem=args.absmem,
            ms_glob_pattern=args.ms_glob_pattern,
            data_column=args.data_column,
            skip_fix_ms=args.skip_fix_ms,
            no_mf_weighting=args.no_mf_weighting,
        )
        client = dask_runner._client
        if client is not None:
            client.close()
        del dask_runner
    else:
        logger.warn("Skipping the image creation step. ")

    if args.imager_only:
        logger.info("Not running any stages after the imager. ")
        return

    # This is the client and pipeline for the RM extraction
    dask_runner_2 = create_dask_runner(
        dask_config=args.dask_config,
    )

    # Define flow
    process_spice.with_options(
        name=f"Arrakis Synthesis -- {args.field}",
    )(args, host, dask_runner_2)


def cli():
    """Command-line interface"""
    # Help string to be shown using the -h option

    descStr = f"""
    {logo_str}

    Arrakis pipeline.

    Before running make sure to start a session of mongodb e.g.
        $ mongod --dbpath=/path/to/database --bind_ip $(hostname -i)

    """

    imager_parser = imager.imager_parser(parent_parser=True)

    # Parse the command line options
    parser = configargparse.ArgParser(
        default_config_files=[".default_config.txt"],
        description=descStr,
        formatter_class=configargparse.ArgumentDefaultsHelpFormatter,
        parents=[imager_parser],
    )
    parser.add("--config", required=False, is_config_file=True, help="Config file path")

    parser.add_argument(
        "field", metavar="field", type=str, help="Name of field (e.g. 2132-50A)."
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
        "--dask_config",
        type=str,
        default=None,
        help="Config file for Dask SlurmCLUSTER.",
    )
    parser.add_argument(
        "--imager_dask_config",
        type=str,
        default=None,
        help="Config file for Dask SlurmCLUSTER.",
    )
    parser.add_argument(
        "--holofile", type=str, default=None, help="Path to holography image"
    )

    parser.add_argument(
        "--yanda",
        type=str,
        default="1.3.0",
        help="Yandasoft version to pull from DockerHub [1.3.0].",
    )

    parser.add_argument(
        "--yanda_image",
        default=None,
        type=Path,
        help="Path to an existing yandasoft singularity container image. ",
    )

    flowargs = parser.add_argument_group("pipeline flow options")
    flowargs.add_argument(
        "--imager_only",
        action="store_true",
        help="Only run the imager component of the pipeline. ",
    )
    flowargs.add_argument(
        "--skip_imager", action="store_true", help="Skip imaging stage [False]."
    )
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
    tools.add_argument(
        "--ionex_server",
        type=str,
        default="ftp://ftp.aiub.unibe.ch/CODE/",
        help="IONEX server [ftp://ftp.aiub.unibe.ch/CODE/].",
    )
    tools.add_argument(
        "--ionex_prefix",
        type=str,
        default="codg",
        help="IONEX prefix.",
    )
    tools.add_argument(
        "--ionex_proxy_server",
        type=str,
        default=None,
        help="Proxy server [None].",
    )
    tools.add_argument(
        "--ionex_formatter",
        type=str,
        default=None,
        help="IONEX formatter [None].",
    )
    tools.add_argument(
        "--ionex_predownload",
        action="store_true",
        help="Pre-download IONEX files [False].",
    )
    cat = parser.add_argument_group("catalogue arguments")
    # Cat args
    cat.add_argument(
        "--outfile", default=None, type=str, help="File to save table to [None]."
    )
    args = parser.parse_args()
    parser.print_values()

    verbose = args.verbose
    if verbose:
        logger.setLevel(logging.INFO)

    logger.info(logo_str)
    logger.info("\n\nArguments: ")
    logger.info(args)

    main(args)


if __name__ == "__main__":
    cli()
