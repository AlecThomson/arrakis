#!/usr/bin/env python3
"""Arrakis multi-field pipeline."""

from __future__ import annotations

import argparse
import logging
from importlib import resources
from pathlib import Path

import configargparse
import yaml
from astropy.time import Time
from prefect import flow

from arrakis import (
    cleanup,
    linmos,
    makecat,
    merge_fields,
    process_spice,
    rmclean_oncuts,
    rmsynth_oncuts,
    validate,
)
from arrakis.logger import UltimateHelpFormatter, logger
from arrakis.utils.database import test_db
from arrakis.utils.pipeline import logo_str
from arrakis.validate import validation_parser


@flow
def process_merge(
    args: argparse.Namespace, host: str, inter_dir: Path, task_runner
) -> None:
    """Workflow to merge spectra from overlapping fields together.

    Args:
        args (Namespace): Parameters to use for this process
        host (str): Address of the mongoDB servicing the processing
        inter_dir (Path): Location to store data from merged fields
        task_runner (TaskRunner): Task runner to use for this process
    """
    previous_future = None
    previous_future = (
        merge_fields.with_options(task_runner=task_runner)(
            fields=args.fields,
            field_dirs=args.datadirs,
            merge_name=args.merge_name,
            output_dir=args.output_dir,
            host=host,
            epoch=args.epoch,
            username=args.username,
            password=args.password,
            yanda=args.yanda,
        )
        if not args.skip_merge
        else previous_future
    )

    previous_future = (
        rmsynth_oncuts.main.with_options(task_runner=task_runner)(
            field=args.field,
            outdir=Path(args.datadir),
            host=args.host,
            epoch=args.epoch,
            sbid=args.sbid,
            username=args.username,
            password=args.password,
            dimension=args.dimension,
            verbose=args.verbose,
            database=args.database,
            limit=args.limit,
            savePlots=args.save_plots,
            weightType=args.weight_type,
            fitRMSF=args.fit_rmsf,
            phiMax_radm2=args.phi_max,
            dPhi_radm2=args.dphi,
            nSamples=args.n_samples,
            polyOrd=args.poly_ord,
            noStokesI=args.no_stokes_i,
            showPlots=args.show_plots,
            not_RMSF=args.not_rmsf,
            rm_verbose=args.rm_verbose,
            debug=args.debug,
            fit_function=args.fit_function,
            tt0=args.tt0,
            tt1=args.tt1,
            ion=False,
            do_own_fit=args.do_own_fit,
        )
        if not args.skip_rmsynth
        else previous_future
    )

    previous_future = (
        rmclean_oncuts.main.with_options(task_runner=task_runner)(
            field=args.merge_name,
            outdir=inter_dir,
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
            showPlots=args.show_plots,
            rm_verbose=args.rm_verbose,
        )
        if not args.skip_rmclean
        else previous_future
    )

    previous_future = (
        makecat.main.with_options(task_runner=task_runner)(
            field=args.merge_name,
            host=host,
            epoch=args.epoch,
            username=args.username,
            password=args.password,
            verbose=args.verbose,
            outfile=args.outfile,
        )
        if not args.skip_cat
        else previous_future
    )
    previous_future = (
        validate.main.with_options(task_runner=task_runner)(
            catalogue_path=Path(args.outfile),
            npix=args.npix,
            map_size=args.map_size,
            snr_cut=args.leakage_snr,
            bins=args.leakage_bins * 2,
        )
        if not args.skip_validate
        else previous_future
    )


def main(args: configargparse.Namespace) -> None:
    """Main script.

    Args:
        args (configargparse.Namespace): Command line arguments.
    """
    if args.dask_config is None:
        with resources.path("arrakis", "configs") as config_dir:
            args.dask_config = config_dir / "default.yaml"

    if args.outfile is None:
        args.outfile = f"{args.merge_name}.pipe.test.fits"

    test_db(
        host=args.host,
        username=args.username,
        password=args.password,
    )

    args_yaml = yaml.dump(vars(args))
    args_yaml_f = Path(f"{args.merge_name}-config-{Time.now().fits}.yaml").absolute()
    logger.info(f"Saving config to '{args_yaml_f}'")
    with args_yaml_f.open("w") as f:
        f.write(args_yaml)

    dask_runner = process_spice.create_dask_runner(
        dask_config=Path(args.dask_config),
    )

    inter_dir = Path(args.output_dir).absolute() / args.merge_name

    process_merge.with_options(
        name=f"Arrakis Merge: {args.merge_name}", task_runner=dask_runner
    )(args, args.host, inter_dir, dask_runner)


def pipeline_parser(parent_parser: bool = False) -> argparse.ArgumentParser:
    """Pipeline parser.

    Args:
        parent_parser (bool, optional): Parent parser. Defaults to False.

    Returns:
        argparse.ArgumentParser: Pipeline parser
    """
    descStr = f"""
    {logo_str}
    Arrakis regional pipeline.

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

    parser.add_argument("--skip_frion", action="store_true", help="Skip cleanup stage.")
    parser.add_argument(
        "--skip_rmsynth", action="store_true", help="Skip RM Synthesis stage."
    )
    parser.add_argument(
        "--skip_rmclean", action="store_true", help="Skip RM-CLEAN stage."
    )
    parser.add_argument("--skip_cat", action="store_true", help="Skip catalogue stage.")
    parser.add_argument(
        "--skip_validate", action="store_true", help="Skip validation stage."
    )
    parser.add_argument(
        "--skip_cleanup", action="store_true", help="Skip cleanup stage."
    )
    return pipeline_parser


def cli():
    """Command-line interface."""
    # Help string to be shown using the -h option

    # Parse the command line options
    pipe_parser = pipeline_parser(parent_parser=True)
    merge_parser = merge_fields.merge_parser(parent_parser=True)
    linmos_parser = linmos.linmos_parser(parent_parser=True)
    common_parser = rmsynth_oncuts.rm_common_parser(parent_parser=True)
    synth_parser = rmsynth_oncuts.rmsynth_parser(parent_parser=True)
    rmclean_parser = rmclean_oncuts.clean_parser(parent_parser=True)
    catalogue_parser = makecat.cat_parser(parent_parser=True)
    val_parser = validation_parser(parent_parser=True)
    clean_parser = cleanup.cleanup_parser(parent_parser=True)
    # Parse the command line options
    parser = configargparse.ArgParser(
        description=pipe_parser.description,
        formatter_class=UltimateHelpFormatter,
        parents=[
            pipe_parser,
            merge_parser,
            linmos_parser,
            common_parser,
            synth_parser,
            rmclean_parser,
            catalogue_parser,
            val_parser,
            clean_parser,
        ],
    )
    parser.add("--config", required=False, is_config_file=True, help="Config file path")
    args = parser.parse_args()

    verbose = args.verbose
    if verbose:
        logger.setLevel(logging.INFO)
    if args.debugger:
        logger.setLevel(logging.DEBUG)

    main(args)


if __name__ == "__main__":
    cli()
