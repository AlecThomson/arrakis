#!/usr/bin/env python3
"""Arrakis single-field pipeline"""
import argparse
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
from arrakis.logger import UltimateHelpFormatter, logger
from arrakis.utils.database import test_db
from arrakis.utils.pipeline import generic_parser, logo_str


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

    if not args.imager_only:
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


def pipeline_parser(parent_parser: bool = False) -> argparse.ArgumentParser:
    descStr = f"""
    {logo_str}

    Arrakis pipeline.

    Before running make sure to start a session of mongodb e.g.
        $ mongod --dbpath=/path/to/database --bind_ip $(hostname -i)

    """
    # Parse the command line options
    pipeline_parser = argparse.ArgumentParser(
        add_help=not parent_parser,
        description=descStr,
        formatter_class=UltimateHelpFormatter,
    )
    parser = pipeline_parser.add_argument_group("pipeline arguments")
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
        "--imager_only",
        action="store_true",
        help="Only run the imager component of the pipeline. ",
    )
    parser.add_argument(
        "--skip_imager", action="store_true", help="Skip imaging stage [False]."
    )
    parser.add_argument(
        "--skip_cutout", action="store_true", help="Skip cutout stage [False]."
    )
    parser.add_argument(
        "--skip_linmos", action="store_true", help="Skip LINMOS stage [False]."
    )
    parser.add_argument(
        "--skip_frion", action="store_true", help="Skip cleanup stage [False]."
    )
    parser.add_argument(
        "--skip_rmsynth", action="store_true", help="Skip RM Synthesis stage [False]."
    )
    parser.add_argument(
        "--skip_rmclean", action="store_true", help="Skip RM-CLEAN stage [False]."
    )
    parser.add_argument(
        "--skip_cat", action="store_true", help="Skip catalogue stage [False]."
    )
    parser.add_argument(
        "--skip_cleanup", action="store_true", help="Skip cleanup stage [False]."
    )

    return pipeline_parser


def cli():
    """Command-line interface"""
    # Help string to be shown using the -h option

    pipe_parser = pipeline_parser(parent_parser=True)
    gen_parser = generic_parser(parent_parser=True)
    imager_parser = imager.imager_parser(parent_parser=True)
    cutout_parser = cutout.cutout_parser(parent_parser=True)
    linmos_parser = linmos.linmos_parser(parent_parser=True)
    common_parser = rmsynth_oncuts.rm_common_parser(parent_parser=True)
    synth_parser = rmsynth_oncuts.rmsynth_parser(parent_parser=True)
    rmclean_parser = rmclean_oncuts.clean_parser(parent_parser=True)
    catalogue_parser = makecat.cat_parser(parent_parser=True)
    clean_parser = cleanup.cleanup_parser(parent_parser=True)
    # Parse the command line options
    parser = configargparse.ArgParser(
        default_config_files=[".default_config.cfg"],
        description=pipe_parser.description,
        formatter_class=UltimateHelpFormatter,
        parents=[
            pipe_parser,
            gen_parser,
            imager_parser,
            cutout_parser,
            linmos_parser,
            common_parser,
            synth_parser,
            rmclean_parser,
            catalogue_parser,
            clean_parser,
        ],
    )

    parser.add("--config", required=False, is_config_file=True, help="Config file path")

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
