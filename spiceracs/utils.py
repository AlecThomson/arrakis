#!/usr/bin/env python
"""Utility functions"""
import copy
import dataclasses
import functools
import json
import logging
import os
import shlex
import stat
import subprocess
import time
import warnings
from dataclasses import asdict, dataclass, make_dataclass
from functools import partial
from glob import glob
from itertools import zip_longest
from os import name
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import astropy.units as u
import dask
import dask.array as da
import dask.distributed as distributed
import numpy as np
import pymongo
from astropy.coordinates import SkyCoord
from astropy.coordinates.angles import dms_tuple, hms_tuple
from astropy.io import fits
from astropy.stats import akaike_info_criterion_lsq
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning
from astropy.wcs import WCS
from casacore.tables import table
from casatasks import listobs
from dask import delayed
from dask.delayed import Delayed
from dask.distributed import Client, get_client
from distributed.client import futures_of
from distributed.diagnostics.progressbar import ProgressBar
from distributed.utils import LoopRunner, is_kernel
from FRion.correct import find_freq_axis
from pymongo.collection import Collection
from scipy.optimize import curve_fit
from scipy.stats import normaltest
from spectral_cube import SpectralCube
from spectral_cube.utils import SpectralCubeWarning
from tornado.ioloop import IOLoop
from tqdm.auto import tqdm, trange

from spiceracs.logger import logger

warnings.filterwarnings(action="ignore", category=SpectralCubeWarning, append=True)
warnings.simplefilter("ignore", category=AstropyWarning)

print = functools.partial(print, flush=True)


def chi_squared(model: np.ndarray, data: np.ndarray, error: np.ndarray) -> float:
    """Calculate chi squared.

    Args:
        model (np.ndarray): Model flux.
        data (np.ndarray): Data flux.
        error (np.ndarray): Data error.

    Returns:
        np.ndarray: Chi squared.
    """
    return np.sum(((model - data) / error) ** 2)


def inspect_client(
    client: Union[distributed.Client, None] = None
) -> Tuple[str, int, int, u.Quantity, int, u.Quantity]:
    """_summary_

    Args:
        client (Union[distributed.Client,None]): Dask client to inspect.
            if None, will use the default client.

    Returns:
        Tuple[ str, int, int, u.Quantity, float, u.Quantity ]: addr, nworkers,
            nthreads, memory, threads_per_worker, memory_per_worker
    """
    """Inspect a client"""
    if client is None:
        client = get_client()
    logger.debug(f"Client: {client}")
    info = client._scheduler_identity
    addr = info.get("address")
    workers = info.get("workers", {})
    nworkers = len(workers)
    nthreads = sum(w["nthreads"] for w in workers.values())
    memory = sum([w["memory_limit"] for w in workers.values()]) * u.byte
    threads_per_worker = nthreads // nworkers
    memory_per_worker = memory / nworkers
    return addr, nworkers, nthreads, memory, threads_per_worker, memory_per_worker


def beam_from_ms(ms: str) -> int:
    """Work out which beam is in this MS"""
    t = table(ms, readonly=True, ack=False)
    vis_feed = t.getcol("FEED1", 0, 1)
    beam = vis_feed[0]
    t.close()
    return beam


def field_idx_from_ms(ms: str) -> int:
    """Get the field from MS metadata"""
    obs = listobs(vis=ms, verbose=True)
    fields = [k for k in obs.keys() if "field_" in k]
    assert len(fields) == 1
    field = fields[0]
    idx = int(field.replace("field_", ""))
    return idx


def wsclean(
    mslist: list,
    use_mpi: bool,
    version: bool = False,
    j: int = None,
    parallel_gridding: int = None,
    parallel_reordering: int = None,
    no_work_on_master: bool = False,
    mem: float = None,
    abs_mem: float = None,
    verbose: bool = False,
    log_time: bool = False,
    quiet: bool = False,
    reorder: bool = False,
    no_reorder: bool = False,
    temp_dir: str = None,
    update_model_required: bool = False,
    no_update_model_required: bool = False,
    no_dirty: bool = False,
    save_first_residual: bool = False,
    save_weights: bool = False,
    save_uv: bool = False,
    reuse_psf: str = None,
    reuse_dirty: str = None,
    apply_primary_beam: bool = False,
    reuse_primary_beam: bool = False,
    use_differential_lofar_beam: bool = False,
    primary_beam_limit: float = None,
    mwa_path: str = None,
    save_psf_pb: bool = False,
    pb_grid_size: int = None,
    beam_model: str = None,
    beam_mode: str = None,
    beam_normalisation_mode: str = None,
    dry_run: bool = False,
    weight: str = None,
    super_weight: float = None,
    mf_weighting: bool = False,
    no_mf_weighting: bool = False,
    weighting_rank_filter: float = None,
    weighting_rank_filter_size: float = None,
    taper_gaussian: str = None,
    taper_tukey: float = None,
    taper_inner_tukey: float = None,
    taper_edge: float = None,
    taper_edge_tukey: float = None,
    use_weights_as_taper: bool = False,
    store_imaging_weights: bool = False,
    name: str = None,
    size: str = None,
    padding: float = None,
    scale: str = None,
    predict: bool = False,
    ws_continue: bool = False,
    subtract_model: bool = False,
    gridder: str = None,
    channels_out: int = None,
    shift: str = None,
    gap_channel_division: bool = False,
    channel_division_frequencies: str = None,
    nwlayers: int = None,
    nwlayers_factor: float = None,
    nwlayers_for_size: str = None,
    no_small_inversion: bool = False,
    small_inversion: bool = False,
    grid_mode: str = None,
    kernel_size: int = None,
    oversampling: int = None,
    make_psf: bool = False,
    make_psf_only: bool = False,
    visibility_weighting_mode: str = None,
    no_normalize_for_weighting: bool = False,
    baseline_averaging: float = None,
    simulate_noise: float = None,
    simulate_baseline_noise: str = None,
    idg_mode: str = None,
    wgridder_accuracy: float = None,
    aterm_config: str = None,
    grid_with_beam: bool = False,
    beam_aterm_update: int = False,
    aterm_kernel_size: float = None,
    apply_facet_solutions: str = None,
    apply_facet_beam: bool = False,
    facet_beam_update: int = False,
    save_aterms: bool = False,
    pol: str = None,
    interval: str = None,
    intervals_out: int = None,
    even_timesteps: bool = False,
    odd_timesteps: bool = False,
    channel_range: str = None,
    field: str = None,
    spws: str = None,
    data_column: str = None,
    maxuvw_m: float = None,
    minuvw_m: float = None,
    maxuv_l: float = None,
    minuv_l: float = None,
    maxw: float = None,
    niter: int = None,
    nmiter: int = None,
    threshold: float = None,
    auto_threshold: float = None,
    auto_mask: float = None,
    force_mask_rounds: int = None,
    local_rms: bool = False,
    local_rms_window: bool = False,
    local_rms_method: bool = False,
    gain: float = None,
    mgain: float = None,
    join_polarizations: bool = False,
    link_polarizations: str = None,
    facet_regions: str = None,
    join_channels: bool = False,
    spectral_correction: str = None,
    no_fast_subminor: bool = False,
    multiscale: bool = False,
    multiscale_scale_bias: bool = False,
    multiscale_max_scales: int = None,
    multiscale_scales: str = None,
    multiscale_shape: str = None,
    multiscale_gain: float = None,
    multiscale_convolution_padding: float = None,
    no_multiscale_fast_subminor: bool = False,
    python_deconvolution: str = None,
    iuwt: bool = False,
    iuwt_snr_test: bool = False,
    no_iuwt_snr_test: bool = False,
    moresane_ext: str = None,
    moresane_arg: str = None,
    moresane_sl: str = None,
    save_source_list: bool = False,
    clean_border: float = None,
    fits_mask: str = None,
    casa_mask: str = None,
    horizon_mask: str = None,
    no_negative: bool = False,
    negative: bool = False,
    stop_negative: bool = False,
    fit_spectral_pol: int = None,
    fit_spectral_log_pol: int = None,
    force_spectrum: str = None,
    deconvolution_channels: int = None,
    squared_channel_joining: bool = False,
    parallel_deconvolution: int = None,
    deconvolution_threads: int = None,
    restore: str = None,
    restore_list: str = None,
    beam_size: float = None,
    beam_shape: str = None,
    fit_beam: bool = False,
    no_fit_beam: bool = False,
    beam_fitting_size: float = None,
    theoretic_beam: bool = False,
    circular_beam: bool = False,
    elliptical_beam: bool = False,
) -> str:
    """Construct a wsclean command.
    If False or None is passed as a parameter, the parameter is not included
    in the command (i.e. wsclean will assume a default value).
    Args:
        mslist (list): List of MSs to be processed.
        use_mpi (bool): Use wsclean-mp for parallel processing.
        version (bool, optional): Print WSClean's version and exit.
            Defaults to False.
        j (int, optional): Specify number of computing threads to use,
            i.e., number of cpu cores that will be used.
            Default: use all cpu cores. to None.
        parallel_gridding (int, optional): Will execute multiple gridders
            simultaneously. This can make things faster in certain cases,
            but will increase memory usage. Defaults to None.
        parallel_reordering (int, optional): Process the reordering with
            multipliple threads. Defaults to None.
        no_work_on_master (bool, optional): In MPI runs, do not use the master
            for gridding. This may be useful if the resources such as memory
            of the master are limited. Defaults to False.
        mem (float, optional): Limit memory usage to the given fraction of the
            total system memory. This is an approximate value.
            Default: 100. Defaults to None.
        abs_mem (float, optional): Like -mem, but this specifies a fixed amount
            of memory in gigabytes. Defaults to None.
        verbose (bool, optional): Increase verbosity of output.
            Defaults to False.
        log_time (bool, optional): Add date and time to each line in the
            output. Defaults to False.
        quiet (bool, optional): Do not output anything but errors.
            Defaults to False.
        reorder (bool, optional): Force reordering of Measurement Set.
            This can be faster when the measurement set needs to be iterated
            several times, such as with many major iterations or in channel
            imaging mode. Default: only reorder when in channel imaging mode.
            Defaults to False.
        no_reorder (bool, optional): Disable reordering of Measurement Set.
            This can be faster when the measurement set needs to be iterated
            several times, such as with many major iterations or in channel
            imaging mode. Default: only reorder when in channel imaging mode.
            Defaults to False.
        temp_dir (str, optional): Set the temporary directory used when
            reordering files. Default: same directory as input measurement set.
            Defaults to None.
        update_model_required (bool, optional): Default. Defaults to False.
        no_update_model_required (bool, optional): These two options specify
            whether the model data column is required to contain valid model
            data after imaging. It can save time to not update the model data
            column. Defaults to False.
        no_dirty (bool, optional): Do not save the dirty image.
            Defaults to False.
        save_first_residual (bool, optional): Save the residual after the
            first iteration. Defaults to False.
        save_weights (bool, optional): Save the gridded weights in the a fits
            file named <image-prefix>-weights.fits. Defaults to False.
        save_uv (bool, optional): Save the gridded uv plane, i.e., the FFT of
            the residual image. The UV plane is complex, hence two images will
            be output: <prefix>-uv-real.fits and <prefix>-uv-imag.fits.
            Defaults to False.
        reuse_psf (str, optional): Load the psf(s) from the given prefix and
            skip the inversion for the psf image. Defaults to None.
        reuse_dirty (str, optional): Load the dirty from the given prefix and
            skip the inversion for the dirty image. Defaults to None.
        apply_primary_beam (bool, optional): Calculate and apply the primary
            beam and save images for the Jones components, with weighting
            identical to the weighting as used by the imager. Only available
            for instruments supported by EveryBeam. Defaults to False.
        reuse_primary_beam (bool, optional): If a primary beam image exists
            on disk, reuse those images. Defaults to False.
        use_differential_lofar_beam (bool, optional): Assume the visibilities
            have already been beam-corrected for the reference direction.
            By default, WSClean will use the information in the measurement
            set to determine if the differential beam should be applied for
            obtaining proper flux levels. Defaults to False.
        primary_beam_limit (float, optional): Level at which to trim the beam
            when performing image-based beam correction,. Default: 0.005.
            Defaults to None.
        mwa_path (str, optional): Set path where to find the MWA beam file(s).
            Defaults to None.
        save_psf_pb (bool, optional): When applying beam correction,
            also save the primary-beam corrected PSF image. Defaults to False.
        pb_grid_size (int, optional): Specify the grid size in number of
            pixels at which to evaluate the primary beam.
            Typically, the primary beam is calculated at a coarse resolution
            grid and interpolated, to reduce the time spent in evaluating the
            beam. This parameter controls the resolution of the grid at which
            to evaluate the primary beam. For rectangular images, pb-grid-size
            indicates the number of pixels along the shortest dimension.
            The total number of pixels in the primary beam grid thus amounts
            to:
                max(width, height) / min(width, height) * pb-grid-size**2.
            Default: 32. Defaults to None.
        beam_model (str, optional): Specify the beam model, only relevant for
            SKA and LOFAR. Available models are Hamaker, Lobes, OskarDipole,
            OskarSphericalWave. Input is case insensitive. Default is Hamaker
            for LOFAR and OskarSphericalWave for SKA. Defaults to None.
        beam_mode (str, optional): [DEBUGGING ONLY] Manually specify the
            beam mode. Only relevant for simulated SKA measurement sets.
            Available modes are array_factor, element and full. Input is case
            insensitive. Default is full. Defaults to None.
        beam_normalisation_mode (str, optional): [DEBUGGING ONLY]
            Manually specify the normalisation of the beam. Only relevant
            for simulated SKA measurement sets. Available modes are none,
            preapplied, full, and amplitude. Default is preapplied.
            Defaults to None.
        dry_run (bool, optional): Parses the command line and quits afterwards.
            No imaging is done. Defaults to False.
        weight (str, optional): Weightmode can be: natural, uniform, briggs.
            Default: uniform. When using Briggs' weighting, add the robustness
            parameter, like: "-weight briggs 0.5". Defaults to None.
        super_weight (float, optional): Increase the weight gridding box size,
            similar to Casa's superuniform weighting scheme. Default: 1.0
            The factor can be rational and can be less than one for subpixel
            weighting. Defaults to None.
        mf_weighting (bool, optional): In spectral mode, calculate the weights
            as if the image was made using MF. This makes sure that the sum of
            channel images equals the MF weights. Otherwise, the channel image
            will become a bit more naturally weighted. This is only relevant
            for weighting modes that require gridding (i.e., Uniform, Briggs').
            Default: off, unless -join-channels is specified.
            Defaults to False.
        no_mf_weighting (bool, optional): Opposite of -ms-weighting;
            can be used to turn off MF weighting in -join-channels mode.
            Defaults to False.
        weighting_rank_filter (float, optional): Filter the weights and set
            high weights to the local mean. The level parameter specifies the
            filter level; any value larger than level*localmean will be set to
            level*localmean. Defaults to None.
        weighting_rank_filter_size (float, optional): Set size of weighting
            rank filter. Default: 16. Defaults to None.
        taper_gaussian (str, optional): Taper the weights with a Gaussian
            function. This will reduce the contribution of long baselines.
            The beamsize is by default in asec, but a unit can be specified
            ("2amin"). Defaults to None.
        taper_tukey (float, optional): Taper the outer weights with a Tukey
            transition. Lambda specifies the size of the transition; use in
            combination with -maxuv-l. Defaults to None.
        taper_inner_tukey (float, optional): Taper the weights with a Tukey
            transition. Lambda specifies the size of the transition; use in
            combination with -minuv-l. Defaults to None.
        taper_edge (float, optional): Taper the weights with a rectangle,
            to keep a space of lambda between the edge and gridded
            visibilities. Defaults to None.
        taper_edge_tukey (float, optional): Taper the edge weights with a Tukey
            window. Lambda is the size of the Tukey transition. When
            -taper-edge is also specified, the Tukey transition starts inside
            the inner rectangle. Defaults to None.
        use_weights_as_taper (bool, optional): Will not use visibility weights
            when determining the imaging weights. This has the effect that e.g.
            uniform weighting can be modified by increasing the visibility
            weight of certain baselines. Without this option, uniform imaging
            weights absorb the visibility weight to make the weighting truly
            uniform. Defaults to False.
        store_imaging_weights (bool, optional): Will store the imaging weights
            in a column named 'IMAGING_WEIGHT_SPECTRUM'. Defaults to False.
        name (str, optional): Use image-prefix as prefix for output files.
            Default is 'wsclean'. Defaults to None.
        size (str, optional): Set the output image size in number of pixels
            (without padding). Defaults to None.
        padding (float, optional): Pad images by the given factor during
            inversion to avoid aliasing. Default: 1.2 (=20%). Defaults to None.
        scale (str, optional): Scale of a pixel. Default unit is degrees, but
            can be specificied, e.g. -scale 20asec. Default: 0.01deg.
            Defaults to None.
        predict (bool, optional): Only perform a single prediction for an
            existing image. Doesn't do any imaging or cleaning. The input
            images should have the same name as the model output images would
            have in normal imaging mode. Defaults to False.
        ws_continue (bool, optional): Will continue an earlier WSClean run.
            Earlier model images will be read and model visibilities will be
            subtracted to create the first dirty residual. CS should have been
            used in the earlier run, and model datashould have been written
            to the measurement set for this to work. Default: off.
            Defaults to False.
        subtract_model (bool, optional): Subtract the model from the
            data column in the first iteration. This can be used to reimage
            an already cleaned image, e.g. at a different resolution.
            Defaults to False.
        gridder (str): Set gridder type: direct-ft, idg, wgridder,
            tuned-wgridder, or wstacking.
        channels_out (int, optional): Splits the bandwidth and makes count
            nr. of images. Default: 1. Defaults to None.
        shift (str, optional): Shift the phase centre to the given location.
            The shift is along the tangential plane. Defaults to None.
        gap_channel_division (bool, optional): In case of irregular frequency
            spacing, this option can be used to not try and split channels to
            make the output channel bandwidth similar, but instead to split
            largest gaps first. Defaults to False.
        channel_division_frequencies (str, optional): Split the bandwidth at
            the specified frequencies (in Hz) before the normal bandwidth
            division is performed. This can e.g. be useful for imaging multiple
            bands with irregular number of channels. Defaults to None.
        nwlayers (int, optional): Number of w-layers to use. Default: minimum
            suggested #w-layers for first MS. Defaults to None.
        nwlayers_factor (float, optional): Use automatic calculation of the
            number of w-layers, but multiple that number by the given factor.
            This can e.g. be useful for increasing w-accuracy. Defaults to None.
        nwlayers_for_size (str, optional): Use the minimum suggested w-layers
            for an image of the given size. Can e.g. be used to increase
            accuracy when predicting small part of full image.
            Defaults to None.
        no_small_inversion (bool, optional): Perform inversion at the Nyquist
            resolution and upscale the image to the requested image size
            afterwards. This speeds up inversion considerably, but makes
            aliasing slightly worse. This effect is in most cases <1%.
            Default: on. Defaults to False.
        small_inversion (bool, optional): Perform inversion at the
            Nyquist resolution and upscale the image to the requested
            image size afterwards. This speeds up inversion considerably,
            but makes aliasing slightly worse. This effect is in most cases
            <1%. Default: on. Defaults to False.
        grid_mode (str, optional): Kernel and mode used for gridding:
            kb = Kaiser-Bessel (default with 7 pixels), nn = nearest neighbour
            (no kernel), more options: rect, kb-no-sinc, gaus, bn. Default: kb.
            Defaults to None.
        kernel_size (int, optional): Gridding antialiasing kernel size.
            Default: 7. Defaults to None.
        oversampling (int, optional): Oversampling factor used during gridding.
            Default: 63. Defaults to None.
        make_psf (bool, optional): Always make the psf, even when no cleaning
            is performed. Defaults to False.
        make_psf_only (bool, optional): Only make the psf, no images are made.
            Defaults to False.
        visibility_weighting_mode (str, optional): Specify visibility weighting
            modi. Affects how the weights (normally) stored in WEIGHT_SPECTRUM
            column are applied. Useful for estimating e.g. EoR power spectra
            errors. Normally one would use this in combination with
            -no-normalize-for-weighting. Defaults to None.
        no_normalize_for_weighting (bool, optional): Disable the normalization
            for the weights, which makes the PSF's peak one.
            See -visibility-weighting-mode. Only useful with natural weighting.
            Defaults to False.
        baseline_averaging (float, optional): Enable baseline-dependent
            averaging. The specified size is in number of wavelengths
            (i.e., uvw-units). One way to calculate this is with
                <baseline in nr. of lambdas> * 2pi *
                <acceptable integration in s> / (24*60*60).
            Defaults to None.
        simulate_noise (float, optional): Will replace every visibility by a
            Gaussian distributed value with given standard deviation before
            imaging. Defaults to None.
        simulate_baseline_noise (str, optional): Like -simulate-noise, but the
            stddevs are provided per baseline, in a text file with antenna1 and
            antenna2 indices and the stddev per line, separated by spaces,
            e.g. "0 1 3.14". Defaults to None.
        idg_mode (str, optional): Sets the IDG mode. Default: cpu. Hybrid is
            recommended when a GPU is available. Defaults to None.
        wgridder_accuracy (float, optional): Set the w-gridding accuracy.
            Default: 1e-4 Useful range: 1e-2 to 1e-. Defaults to None.
        aterm_config (str, optional): Specify a parameter set describing how
            a-terms should be applied. Please refer to the documentation for
            details of the configuration file format. Applying a-terms is only
            possible when IDG is enabled. Defaults to None.
        grid_with_beam (bool, optional): Apply a-terms to correct for the
            primary beam. This is only possible when IDG is enabled.
            Defaults to False.
        beam_aterm_update (int, optional): Set the ATerm update time in
            seconds. The default is every 300 seconds. It also sets the
            interval over which to calculate the primary beam when using
            -apply-primary-beam when not gridding with the beam.
            Defaults to False.
        aterm_kernel_size (float, optional): Kernel size reserved for aterms
            by IDG. Defaults to None.
        apply_facet_solutions (str, optional): Apply solutions from the
            provided (h5) file per facet when gridding facet based images.
            Provided file is assumed to be in H5Parm format. Filename is
            followed by a comma separated list of strings specifying which sol
            tabs from the provided H5Parm file are used. Defaults to None.
        apply_facet_beam (bool, optional): Apply beam gains to facet center
            when gridding facet based image. Defaults to False.
        facet_beam_update (int, optional): Set the facet beam update time in
            seconds. The default is every 120 seconds. Defaults to False.
        save_aterms (bool, optional): Output a fits file for every aterm
            update, containing the applied image for every station.
            Defaults to False.
        pol (str, optional): Default: 'I'.
            Possible values: XX, XY, YX, YY, I, Q, U, V, RR, RL, LR or LL
            (case insensitive). It is allowed but not necessary to separate
            with commas, e.g.: 'xx,xy,yx,yy'.Two or four polarizations can be
            joinedly cleaned (see '-joinpolarizations'), but this is not the
            default. I, Q, U and V polarizations will be directly calculated
            from the visibilities, which might require correction to get to
            real IQUV values. The 'xy' polarization will output both a real
            and an imaginary image, which allows calculating true Stokes
            polarizations for those telescopes. Defaults to None.
        interval (str, optional): Only image the given time interval. Indices
            specify the timesteps, end index is exclusive. Default: image all
            time steps. Defaults to None.
        intervals_out (int, optional): Number of intervals to image inside the
            selected global interval. Default: . Defaults to None.
        even_timesteps (bool, optional): Only select even timesteps. Can be
            used together with -odd-timesteps to determine noise values.
            Defaults to False.
        odd_timesteps (bool, optional): Only select odd timesteps.
            Defaults to False.
        channel_range (str, optional): Only image the given channel range.
            Indices specify channel indices, end index is exclusive.
            Default: image all channels. Defaults to None.
        field (str, optional): Image the given field id(s). A comma-separated
            list of field ids can be provided. When multiple fields are given,
            all fields should have the same phase centre.
            Specifying '-field all' will image all fields in the
            measurement set. Default: first field (id 0). Defaults to None.
        spws (str, optional): Selects only the spws given in the list. list
            should be a comma-separated list of integers. Default: all spws.
            Defaults to None.
        data_column (str, optional): Default: CORRECTED_DATA if it exists,
            otherwise DATA will be used. Defaults to None.
        maxuvw_m (float, optional): Set the min max uv distance in lambda.
            Defaults to None.
        minuvw_m (float, optional): Set the max baseline distance in meters.
            Defaults to None.
        maxuv_l (float, optional): Set the min max uv distance in lambda.
            Defaults to None.
        minuv_l (float, optional): Set the max uv distance in lambda.
            Defaults to None.
        maxw (float, optional): Do not grid visibilities with a w-value
            higher than the given percentage of the max w, to save speed.
            Default: grid everythin. Defaults to None.
        niter (int, optional): Maximum number of clean iterations to perform.
            Default: 0 (=no cleaning). Defaults to None.
        nmiter (int, optional): Maximum number of major clean
            (inversion/prediction) iterations. Default: 20.A value of 0 means
            no limit. Defaults to None.
        threshold (float, optional): Stopping clean thresholding in Jy.
            Default: 0.0. Defaults to None.
        auto_threshold (float, optional): Estimate noise level using a robust
            estimator and stop at sigma x stddev. Defaults to None.
        auto_mask (float, optional): Construct a mask from found components
            and when a threshold of sigma is reached, continue cleaning with
            the mask down to the normal threshold. Defaults to None.
        force_mask_rounds (int, optional): Will force the derivation of the
            mask to be carried out across a set number of major cleaning rounds.
        local_rms (bool, optional): Instead of using a single RMS for auto
            thresholding/masking, use a spatially varying RMS image.
            Defaults to False.
        local_rms_window (bool, optional): Size of window for creating the
            RMS background map, in number of PSFs. Default: 25 psfs.
            Defaults to False.
        local_rms_method (bool, optional): Either 'rms'
            (default, uses sliding window RMS) or 'rms-with-min'
            (use max(window rms, 0.3 x window min)). Defaults to False.
        gain (float, optional): Cleaning gain: Ratio of peak that will be
            subtracted in each iteration. Default: 0.1. Defaults to None.
        mgain (float, optional): Cleaning gain for major iterations: Ratio of
            peak that will be subtracted in each major iteration.
            To use major iterations, 0.85 is a good value.
            Default: 1.0. Defaults to None.
        join_polarizations (bool, optional): Perform deconvolution by
            searching for peaks in the sum of squares of the polarizations,
            but subtract components from the individual images. Only possible
            when imaging two or four Stokes or linear parameters.
            Default: off. Defaults to False.
        link_polarizations (str, optional): Links all polarizations to be
            cleaned from the given list: components are found in the given
            list, but cleaned from all polarizations.  Defaults to None.
        facet_regions (str, optional): Split the image into facets using the
            facet regions defined in  the facets.reg file. Default: off.
            Defaults to None.
        join_channels (bool, optional): Perform deconvolution by searching for
            peaks in the MF image, but subtract components from individual
            channels. This will turn on mf-weighting by default.
            Default: off. Defaults to False.
        spectral_correction (str, optional): Enable correction of the given
            spectral function inside deconvolution. This can e.g. avoid
            downweighting higher frequencies because of reduced flux density.
            1st term is total flux, 2nd is si, 3rd curvature, etc.
            Example: -spectral-correction 150e6 83.084,-0.699,-0.110
            Defaults to None.
        no_fast_subminor (bool, optional): Do not use the subminor loop
            optimization during (non-multiscale) cleaning.
            Default: use the optimization. Defaults to False.
        multiscale (bool, optional): Clean on different scales.
            This is a new algorithm. Default: off. This parameter invokes the
            optimized multiscale algorithm published by Offringa & Smirnov
            (2017). Defaults to False.
        multiscale_scale_bias (bool, optional): Parameter to prevent cleaning
            small scales in the large-scale iterations. A lower bias will give
            more focus to larger scales. Default: 0.6 Defaults to False.
        multiscale_max_scales (int, optional): Set the maximum number of scales
            that WSClean should use in multiscale cleaning. Only relevant when
            -multiscale-scales is not set. Default: unlimited.
            Defaults to None.
        multiscale_scales (str, optional): Sets a list of scales to use in
            multi-scale cleaning. If unset, WSClean will select the delta
            (zero) scale, scales starting at four times the synthesized PSF,
            and increase by a factor of two until the maximum scale is reached
            or the maximum number of scales is reached.
            Example: -multiscale-scales 0,5,12.5 Defaults to None.
        multiscale_shape (str, optional): Sets the shape function used during
            multi-scale clean. Either 'tapered-quadratic' (default) or
            'gaussian'. Defaults to None.
        multiscale_gain (float, optional): Size of step made in the subminor
            loop of multi-scale. Default currently 0.2, but shows sign of
            instability. A value of 0.1 might be more stable. Defaults to None.
        multiscale_convolution_padding (float, optional): Size of zero-padding
            for convolutions during the multi-scale cleaning.
            Default: 1.1 Defaults to None.
        no_multiscale_fast_subminor (bool, optional): Disable the
            'fast subminor loop' optimization, that will only search a part of
            the image during the multi-scale subminor loop. The optimization
            is on by default. Defaults to False.
        python_deconvolution (str, optional): Run a custom deconvolution
            algorithm written in Python. See manual for the interface.
            Defaults to None.
        iuwt (bool, optional): Use the IUWT deconvolution algorithm.
            Defaults to False.
        iuwt_snr_test (bool, optional): Stop IUWT when the SNR decreases.
            This might help limitting divergence, but can occasionally also
            stop the algorithm too early. Default: no SNR test.
            Defaults to False.
        no_iuwt_snr_test (bool, optional): Do not stop IUWT when the SNR
            decreases. This might help limitting divergence, but can
            occasionally also stop the algorithm too early.
            Default: no SNR test. Defaults to False.
        moresane_ext (str, optional): Use the MoreSane deconvolution algorithm,
            installed at the specified location. Defaults to None.
        moresane_arg (str, optional): Pass the specified arguments to moresane.
            Note that multiple parameters have to be enclosed in quotes.
            Defaults to None.
        moresane_sl (str, optional): MoreSane --sigmalevel setting for each
            major loop iteration. Useful to start at high levels and go down
            with subsequent loops, e.g. 20,10,5 Defaults to None.
        save_source_list (bool, optional): Saves the found clean components
            as a BBS/DP3 text sky model. This parameter enables Gaussian shapes
            during multi-scale cleaning (-multiscale-shape gaussian).
            Defaults to False.
        clean_border (float, optional): Set the border size in which no
            cleaning is performed, in percentage of the width/height of the
            image. With an image size of 1000 and clean border of 1%,
            each border is 10 pixels. Default: 0% Defaults to None.
        fits_mask (str, optional): Use the specified fits-file as mask during
            cleaning. Defaults to None.
        casa_mask (str, optional): Use the specified CASA mask as mask
        during cleaning. Defaults to None.
        horizon_mask (str, optional): Use a mask that avoids cleaning emission
            beyond the horizon. Distance is an angle (e.g. "5deg") that
            (when positive) decreases the size of the mask to stay further away
            from the horizon. Defaults to None.
        no_negative (bool, optional): Do not allow negative components during
            cleaning. Not the default. Defaults to False.
        negative (bool, optional): Default on: opposite of -nonegative.
            Defaults to False.
        stop_negative (bool, optional): Stop on negative components.
            Not the default. Defaults to False.
        fit_spectral_pol (int, optional): Fit a polynomial over frequency to
            each clean component. This has only effect when the channels are
            joined with -join-channels. Defaults to None.
        fit_spectral_log_pol (int, optional): Like fit-spectral-pol, but fits
            a logarithmic polynomial over frequency instead. Defaults to None.
        force_spectrum (str, optional): Uses the fits file to force spectral
            indices (or other/more terms)during the deconvolution.
            Defaults to None.
        deconvolution_channels (int, optional): Decrease the number of channels
            as specified by -channels-out to the given number for
            deconvolution. Only possible in combination with one of the
            -fit-spectral options. Proper residuals/restored images will
            only be returned when mgain < 1. Defaults to None.
        squared_channel_joining (bool, optional): Use with -join-channels to
            perform peak finding in the sum of squared values over channels,
            instead of the normal sum. This is useful for imaging QU
            polarizations with non-zero rotation measures, for which the normal
             sum is insensitive. Defaults to False.
        parallel_deconvolution (int, optional): Deconvolve subimages in
            parallel. Subimages will be at most of the given size.
            Defaults to None.
        deconvolution_threads (int, optional): Number of threads to use during
            deconvolution. On machines with a large nr of cores, this may be
             used to decrease the memory usage. Defaults to None.
        restore (str, optional): Restore the model image onto the residual
            image and save it in output image. By default, the beam parameters
             are read from the residual image. If this parameter is given,
              wsclean will do the restoring and then exit:
            no cleaning is performed. Defaults to None.
        restore_list (str, optional): Restore a source list onto the residual
            image and save it in output image. Except for the model input
            format, this parameter behaves equal to -restore. Defaults to None.
        beam_size (float, optional): Set a circular beam size (FWHM) in arcsec
            for restoring the clean components. This is the same as
            -beam-shape <size> <size> 0. Defaults to None.
        beam_shape (str, optional): Set the FWHM beam shape for restoring the
            clean components. Defaults units for maj and min are arcsec, and
            degrees for PA. Can be overriden,
            e.g. '-beam-shape 1amin 1amin 3deg'.
            Default: shape of PSF. Defaults to None.
        fit_beam (bool, optional): Determine beam shape by fitting the PSF
            (default if PSF is made). Defaults to False.
        no_fit_beam (bool, optional): Do not determine beam shape from the PSF.
            Defaults to False.
        beam_fitting_size (float, optional): Use a fitting box the size of
            <factor> times the theoretical beam size for fitting a Gaussian
            to the PSF. Defaults to None.
        theoretic_beam (bool, optional): Write the beam in output fits files as
            calculated from the longest projected baseline. This method results
            in slightly less accurate beam size/integrated fluxes, but provides
             a beam size without making the PSF for quick imaging.
             Default: off. Defaults to False.
        circular_beam (bool, optional): Force the beam to be circular:
            bmin will be set to bmaj. Defaults to False.
        elliptical_beam (bool, optional): Allow the beam to be elliptical.
            Default. Defaults to False.
    Returns:
        str: WSClean command
    """

    arguments = copy.deepcopy(locals())
    mslist = arguments.pop("mslist")
    use_mpi = arguments.pop("use_mpi")
    # Check for MPI
    if use_mpi:
        command = "mpirun wsclean-mp"
    else:
        command = "wsclean "

    # Check for square channels and multiscale
    if arguments["squared_channel_joining"] and arguments["multiscale"]:
        logger.info("CAUTION - square channel joining and multiscale is unstable!")

    for key, value in arguments.items():
        if type(value) is bool:
            if value:
                command += f" -{key.replace('_', '-')}"
        elif value:
            if "ws_" in key:  # Catch for ws_continue command
                key.replace("ws_", "")
            command += f" -{key.replace('_','-')} {value}"
    command += f" {' '.join(mslist)}"
    return command


def best_aic_func(aics: np.ndarray, n_param: np.ndarray) -> Tuple[float, int, int]:
    """Find the best AIC for a set of AICs using Occam's razor."""
    # Find the best AIC
    best_aic_idx = int(np.nanargmin(aics))
    best_aic = float(aics[best_aic_idx])
    best_n = int(n_param[best_aic_idx])
    logger.debug(f"Lowest AIC is {best_aic}, with {best_n} params.")
    # Check if lower have diff < 2 in AIC
    aic_abs_diff = np.abs(aics - best_aic)
    bool_min_idx = np.zeros_like(aics).astype(bool)
    bool_min_idx[best_aic_idx] = True
    potential_idx = (aic_abs_diff[~bool_min_idx] < 2) & (
        n_param[~bool_min_idx] < best_n
    )
    if not any(potential_idx):
        return best_aic, best_n, best_aic_idx

    bestest_n = int(np.min(n_param[~bool_min_idx][potential_idx]))
    bestest_aic_idx = int(np.where(n_param == bestest_n)[0][0])
    bestest_aic = float(aics[bestest_aic_idx])
    logger.debug(
        f"Model within 2 of lowest AIC found. Occam says to take AIC of {bestest_aic}, with {bestest_n} params."
    )
    return bestest_aic, bestest_n, bestest_aic_idx


# Stolen from GLEAM-X - thanks Uncle Timmy!
def power_law(nu: np.ndarray, norm: float, alpha: float, ref_nu: float) -> np.ndarray:
    """A power law model.

    Args:
        nu (np.ndarray): Frequency array.
        norm (float): Reference flux.
        alpha (float): Spectral index.
        ref_nu (float): Reference frequency.

    Returns:
        np.ndarray: Model flux.
    """
    return norm * (nu / ref_nu) ** alpha


def flat_power_law(nu: np.ndarray, norm: float, ref_nu: float) -> np.ndarray:
    """A flat power law model.

    Args:
        nu (np.ndarray): Frequency array.
        norm (float): Reference flux.
        ref_nu (float): Reference frequency.

    Returns:
        np.ndarray: Model flux.
    """
    x = ref_nu * np.ones_like(nu)
    return norm * x


def curved_power_law(
    nu: np.ndarray, norm: float, alpha: float, beta: float, ref_nu: float
) -> np.ndarray:
    """A curved power law model.

    Args:
        nu (np.ndarray): Frequency array.
        norm (float): Reference flux.
        alpha (float): Spectral index.
        beta (float): Spectral curvature.
        ref_nu (float): Reference frequency.

    Returns:
        np.ndarray: Model flux.
    """
    x = nu / ref_nu
    power = alpha + beta * np.log10(x)
    return norm * x**power


def fit_pl(
    freq: np.ndarray, flux: np.ndarray, fluxerr: np.ndarray, nterms: int
) -> dict:
    """Perform a power law fit to a spectrum.

    Args:
        freq (np.ndarray): Frequency array.
        flux (np.ndarray): Flux array.
        fluxerr (np.ndarray): Error array.
        nterms (int): Number of terms to use in the fit.

    Returns:
        dict: Best fit parameters.
    """
    try:
        goodchan = np.logical_and(
            np.isfinite(flux), np.isfinite(fluxerr)
        )  # Ignore NaN channels!
        ref_nu = np.nanmean(freq[goodchan])
        p0_long = (np.median(flux[goodchan]), -0.8, 0.0)
        model_func_dict = {
            0: partial(flat_power_law, ref_nu=ref_nu),
            1: partial(power_law, ref_nu=ref_nu),
            2: partial(curved_power_law, ref_nu=ref_nu),
        }

        # Initialise the save dict
        save_dict = {
            n: {} for n in range(nterms + 1)
        }  # type: Dict[int, Dict[str, Any]]
        for n in range(nterms + 1):
            p0 = p0_long[: n + 1]
            save_dict[n]["aics"] = np.nan
            save_dict[n]["params"] = np.ones_like(p0) * np.nan
            save_dict[n]["errors"] = np.ones_like(p0) * np.nan
            save_dict[n]["models"] = np.ones_like(freq)
            save_dict[n]["highs"] = np.ones_like(freq)
            save_dict[n]["lows"] = np.ones_like(freq)
            # 4 possible flags
            save_dict[n]["fit_flags"] = {
                "is_negative": True,
                "is_not_finite": True,
                "is_not_normal": True,
                "is_close_to_zero": True,
            }

        # Now iterate over the number of terms
        for n in range(nterms + 1):
            p0 = p0_long[: n + 1]
            model_func = model_func_dict[n]
            try:
                fit_res = curve_fit(
                    model_func,
                    freq[goodchan],
                    flux[goodchan],
                    p0=p0,
                    sigma=fluxerr[goodchan],
                    absolute_sigma=True,
                )
            except RuntimeError:
                logger.critical(f"Failed to fit {n}-term power law")
                continue

            best, covar = fit_res
            model_arr = model_func(freq, *best)
            model_high = model_func(freq, *(best + np.sqrt(np.diag(covar))))
            model_low = model_func(freq, *(best - np.sqrt(np.diag(covar))))
            model_err = model_high - model_low
            ssr = np.sum((flux[goodchan] - model_arr[goodchan]) ** 2)
            aic = akaike_info_criterion_lsq(ssr, len(p0), goodchan.sum())

            # Save the results
            save_dict[n]["aics"] = aic
            save_dict[n]["params"] = best
            save_dict[n]["errors"] = np.sqrt(np.diag(covar))
            save_dict[n]["models"] = model_arr
            save_dict[n]["highs"] = model_high
            save_dict[n]["lows"] = model_low

            # Calculate the flags
            # Flag if model is negative
            is_negative = (model_arr < 0).any()
            if is_negative:
                logger.warning(f"Stokes I flag: Model {n} is negative")
            # Flag if model is NaN or Inf
            is_not_finite = ~np.isfinite(model_arr).all()
            if is_not_finite:
                logger.warning(f"Stokes I flag: Model {n} is not finite")
            # # Flag if model and data are statistically different
            residuals = flux[goodchan] - model_arr[goodchan]
            # Assume errors on resdiuals are the same as the data
            # i.e. assume the model has no error
            residuals_err = fluxerr[goodchan]
            residuals_norm = residuals / residuals_err
            # Test if the residuals are normally distributed
            ks, pval = normaltest(residuals_norm)
            is_not_normal = pval < 1e-6  # 1 in a million chance of being unlucky
            if is_not_normal:
                logger.warning(
                    f"Stokes I flag: Model {n} is not normally distributed - {pval=}, {ks=}"
                )

            # Test if model is close to 0 within 1 sigma
            is_close_to_zero = (model_arr[goodchan] / fluxerr[goodchan] < 1).any()
            if is_close_to_zero:
                logger.warning(f"Stokes I flag: Model {n} is close (1sigma) to 0")
            fit_flag = {
                "is_negative": is_negative,
                "is_not_finite": is_not_finite,
                "is_not_normal": is_not_normal,
                "is_close_to_zero": is_close_to_zero,
            }
            save_dict[n]["fit_flags"] = fit_flag
            logger.debug(f"{n}: {aic}")

        # Now find the best model
        best_aic, best_n, best_aic_idx = best_aic_func(
            np.array([save_dict[n]["aics"] for n in range(nterms + 1)]),
            np.array([n for n in range(nterms + 1)]),
        )
        logger.debug(f"Best fit: {best_n}, {best_aic}")
        best_p = save_dict[best_n]["params"]
        best_e = save_dict[best_n]["errors"]
        best_m = save_dict[best_n]["models"]
        best_f = model_func_dict[best_n]
        best_flag = save_dict[best_n]["fit_flags"]
        best_h = save_dict[best_n]["highs"]
        best_l = save_dict[best_n]["lows"]
        chi_sq = chi_squared(
            model=best_m[goodchan],
            data=flux[goodchan],
            error=fluxerr[goodchan],
        )
        chi_sq_red = chi_sq / (goodchan.sum() - len(best_p))
        return dict(
            best_n=best_n,
            best_p=best_p,
            best_e=best_e,
            best_m=best_m,
            best_h=best_h,
            best_l=best_l,
            best_f=best_f,
            fit_flag=best_flag,
            ref_nu=ref_nu,
            chi_sq=chi_sq,
            chi_sq_red=chi_sq_red,
        )
    except Exception as e:
        logger.critical(f"Failed to fit power law: {e}")
        return dict(
            best_n=np.nan,
            best_p=[np.nan],
            best_e=[np.nan],
            best_m=np.ones_like(freq),
            best_h=np.ones_like(freq),
            best_l=np.ones_like(freq),
            best_f=None,
            fit_flag={
                "is_negative": True,
                "is_not_finite": True,
                "is_not_normal": True,
                "is_close_to_zero": True,
            },
            ref_nu=np.nan,
            chi_sq=np.nan,
            chi_sq_red=np.nan,
        )


# stolen from https://stackoverflow.com/questions/32954486/zip-iterators-asserting-for-equal-length-in-python
def zip_equal(*iterables):
    sentinel = object()
    for combo in zip_longest(*iterables, fillvalue=sentinel):
        if sentinel in combo:
            raise ValueError("Iterables have different lengths")
        yield combo


def chunk_dask(
    outputs: list,
    batch_size: int = 10_000,
    task_name="",
    progress_text="",
    verbose=True,
) -> list:
    client = get_client()
    chunk_outputs = []
    for i in trange(
        0, len(outputs), batch_size, desc=f"Chunking {task_name}", disable=(not verbose)
    ):
        outputs_chunk = outputs[i : i + batch_size]
        futures = client.persist(outputs_chunk)
        # dumb solution for https://github.com/dask/distributed/issues/4831
        if i == 0:
            logger.debug("I sleep!")
            time.sleep(10)
            logger.debug("I awake!")
        tqdm_dask(futures, desc=progress_text, disable=(not verbose))
        chunk_outputs.extend(futures)
    return chunk_outputs


def latexify(fig_width=None, fig_height=None, columns=1):
    """Set up matplotlib's RC params for LaTeX plotting.
    Call this before plotting a figure.

    Parameters
    ----------
    fig_width : float, optional, inches
    fig_height : float,  optional, inches
    columns : {1, 2}
    """
    from math import sqrt

    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    # code adapted from http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples
    # Width and max height in inches for IEEE journals taken from
    # computer.org/cms/Computer.org/Journal%20templates/transactions_art_guide.pdf

    assert columns in [1, 2]

    if fig_width is None:
        fig_width = 3.39 if columns == 1 else 6.9  # width in inches

    if fig_height is None:
        golden_mean = (sqrt(5) - 1.0) / 2.0  # Aesthetic ratio
        fig_height = fig_width * golden_mean  # height in inches

    MAX_HEIGHT_INCHES = 8.0
    if fig_height > MAX_HEIGHT_INCHES:
        logger.waning(
            "fig_height too large:"
            + fig_height
            + "so will reduce to"
            + MAX_HEIGHT_INCHES
            + "inches."
        )
        fig_height = MAX_HEIGHT_INCHES

    params = {
        "backend": "pdf",
        "axes.labelsize": 8,  # fontsize for x and y labels (was 10)
        "axes.titlesize": 8,
        "font.size": 8,  # was 10
        "legend.fontsize": 8,  # was 10
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "text.usetex": False,
        "figure.figsize": [fig_width, fig_height],
        "font.family": "serif",
    }

    matplotlib.rcParams.update(params)


def delayed_to_da(
    list_of_delayed: List[Delayed], chunk: Union[int, None] = None
) -> da.Array:
    """Convert list of delayed arrays to a dask array

    Args:
        list_of_delayed (List[delayed]): List of delayed objects
        chunk (int, optional): Chunksize to use. Defaults to None.

    Returns:
        da.Array: Dask array
    """
    sample = list_of_delayed[0].compute()
    dim = (len(list_of_delayed),) + sample.shape
    if chunk is None:
        c_dim = dim
    else:
        c_dim = (chunk,) + sample.shape
    darray_list = [
        da.from_delayed(lazy, dtype=sample.dtype, shape=sample.shape)
        for lazy in list_of_delayed
    ]
    darray = da.stack(darray_list, axis=0).reshape(dim).rechunk(c_dim)

    return darray


def yes_or_no(question: str) -> bool:
    """Ask a yes or no question via input()

    Args:
        question (str): Question to ask

    Returns:
        bool: True for yes, False for no
    """
    while "Please answer 'y' or 'n'":
        reply = str(input(question + " (y/n): ")).lower().strip()
        if reply[:1] == "y":
            ret = True
        if reply[:1] == "n":
            ret = False
    return ret


def fix_header(cutout_header: fits.Header, original_header: fits.Header) -> fits.Header:
    """Make cutout header the same as original header

    Args:
        cutout_header (fits.Header): Cutout header
        original_header (fits.Header): Original header

    Returns:
        fits.Header: Fixed header
    """
    axis_cut = find_freq_axis(cutout_header)
    axis_orig = find_freq_axis(original_header)
    fixed_header = cutout_header.copy()
    if axis_cut != axis_orig:
        for key, val in cutout_header.items():
            if key[-1] == str(axis_cut):
                fixed_header[f"{key[:-1]}{axis_orig}"] = val
                fixed_header[key] = original_header[key]

    return fixed_header


def deg_to_hms(deg: float) -> hms_tuple:
    """Convert degree to hms without astropy.

    Args:
        deg (float): Decimal degrees

    Returns:
        hms_tuple: HMS, like coord.ra.hms
    """
    h_per_d = 24 / 360
    hours = deg * h_per_d
    hour = float(int(hours))
    minutes = (hours - hour) * 60
    minute = float(int(minutes))
    seconds = (minutes - minute) * 60
    return hms_tuple(hour, minute, seconds)


def deg_to_dms(deg: float) -> dms_tuple:
    """Convert degree to hms without astropy.

    Args:
        deg (float): Decimal degrees

    Returns:
        hms_tuple: DMS, like coord.dec.dms
    """
    degree = float(int(deg))
    minutes = (deg - degree) * 60
    minute = float(int(minutes))
    seconds = (minutes - minute) * 60
    return dms_tuple(degree, minute, seconds)


def coord_to_string(coord: SkyCoord) -> Tuple[str, str]:
    """Convert coordinate to string without astropy

    Args:
        coord (SkyCoord): Coordinate

    Returns:
        Tuple[str,str]: Tuple of RA string, Dec string
    """
    ra = coord.ra
    dec = coord.dec

    ra_hms = deg_to_hms(ra.value)
    dec_dms = deg_to_dms(dec.value)

    ra_str = f"{ra_hms.h:02.0f}:{ra_hms.m:02.0f}:{ra_hms.s:06.3f}"
    dec_str = f"{dec_dms.d:02.0f}:{abs(dec_dms.m):02.0f}:{abs(dec_dms.s):05.2f}"
    return ra_str, dec_str


def test_db(
    host: str, username: Union[str, None] = None, password: Union[str, None] = None
) -> bool:
    """Test connection to MongoDB

    Args:
        host (str): Mongo host IP.
        username (str, optional): Mongo username. Defaults to None.
        password (str, optional): Mongo password. Defaults to None.
        verbose (bool, optional): Verbose output. Defaults to True.


    Returns:
        bool: True if connection succesful

    Raises:
        Exception: If connection fails.
    """
    logger.info("Testing MongoDB connection...")
    # default connection (ie, local)
    with pymongo.MongoClient(
        host=host,
        connect=False,
        username=username,
        password=password,
        authMechanism="SCRAM-SHA-256",
    ) as dbclient:  # type: pymongo.MongoClient
        try:
            dbclient.list_database_names()
        except pymongo.errors.ServerSelectionTimeoutError:
            raise Exception("Please ensure 'mongod' is running")

        logger.info("MongoDB connection succesful!")

    return True


def get_db(
    host: str, username: Union[str, None] = None, password: Union[str, None] = None
) -> Tuple[Collection, Collection, Collection,]:
    """Get MongoDBs

    Args:
        host (str): Mongo host IP.
        username (str, optional): Username. Defaults to None.
        password (str, optional): Password. Defaults to None.

    Returns:
        Tuple[Collection, Collection, Collection]: beams_col, island_col, comp_col
    """
    dbclient = pymongo.MongoClient(
        host=host,
        connect=False,
        username=username,
        password=password,
        authMechanism="SCRAM-SHA-256",
    )  # type: pymongo.MongoClient
    mydb = dbclient["spiceracs"]  # Create/open database
    comp_col = mydb["components"]  # Create/open collection
    island_col = mydb["islands"]  # Create/open collection
    beams_col = mydb["beams"]  # Create/open collection
    return beams_col, island_col, comp_col


def get_field_db(host: str, username=None, password=None) -> Collection:
    """Get MongoDBs

    Args:
        host (str): Mongo host IP.
        username (str, optional): Username. Defaults to None.
        password (str, optional): Password. Defaults to None.

    Returns:
        pymongo.Collection: beams_col, island_col, comp_col
    """
    dbclient = pymongo.MongoClient(
        host=host,
        connect=False,
        username=username,
        password=password,
        authMechanism="SCRAM-SHA-256",
    )  # type: pymongo.MongoClient
    mydb = dbclient["spiceracs"]  # Create/open database
    field_col = mydb["fields"]  # Create/open collection
    return field_col


# stolen from https://github.com/tqdm/tqdm/issues/278
class TqdmProgressBar(ProgressBar):
    """Tqdm for Dask"""

    def __init__(
        self,
        keys,
        scheduler=None,
        interval="100ms",
        loop=None,
        complete=True,
        start=True,
        **tqdm_kwargs,
    ):
        super(TqdmProgressBar, self).__init__(keys, scheduler, interval, complete)
        self.tqdm = tqdm(keys, **tqdm_kwargs)
        self.loop = loop or IOLoop()

        if start:
            loop_runner = LoopRunner(self.loop)
            loop_runner.run_sync(self.listen)

    def _draw_bar(self, remaining, all, **kwargs):
        update_ct = (all - remaining) - self.tqdm.n
        self.tqdm.update(update_ct)

    def _draw_stop(self, **kwargs):
        self.tqdm.close()


def tqdm_dask(futures: distributed.Future, **kwargs) -> None:
    """Tqdm for Dask futures"""
    futures = futures_of(futures)
    if not isinstance(futures, (set, list)):
        futures = [futures]
    TqdmProgressBar(futures, **kwargs)


def port_forward(port: int, target: str) -> None:
    """Forward ports to local host

    Args:
        port (int): port to forward
        target (str): Target host
    """
    logger.info(f"Forwarding {port} from localhost to {target}")
    cmd = f"ssh -N -f -R {port}:localhost:{port} {target}"
    command = shlex.split(cmd)
    output = subprocess.Popen(command)


def try_mkdir(dir_path: str, verbose=True):
    """Create directory if it doesn't exist

    Args:
        dir_path (str): Path to directory
        verbose (bool, optional): Verbose output. Defaults to True.
    """
    # Create output dir if it doesn't exist
    try:
        os.mkdir(dir_path)
        logger.info(f"Made directory '{dir_path}'.")
    except FileExistsError:
        logger.info(f"Directory '{dir_path}' exists.")


def try_symlink(src: str, dst: str, verbose=True):
    """Create symlink if it doesn't exist

    Args:
        src (str): Source path
        dst (str): Destination path
        verbose (bool, optional): Verbose output. Defaults to True.
    """
    # Create output dir if it doesn't exist
    try:
        os.symlink(src, dst)
        logger.info(f"Made symlink '{dst}'.")
    except FileExistsError:
        logger.info(f"Symlink '{dst}' exists.")


def head2dict(h: fits.Header) -> Dict[str, Any]:
    """Convert FITS header to a dict.

    Writes a cutout, as stored in source_dict, to disk. The file location
    should already be specified in source_dict. This format is intended
    for parallel use with pool.map syntax.

    Args:
        h: An astropy FITS header.

    Returns:
        data (dict): The FITS head converted to a dict.

    """
    data = {}
    for c in h.__dict__["_cards"]:
        if c[0] == "":
            continue
        data[c[0]] = c[1]
    return data


class MyEncoder(json.JSONEncoder):
    """Cutom JSON encorder.

    Parses the data stored in source_dict to JSON without
    errors.

    """

    def default(self, obj):  # pylint: disable=E0202
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.complex):
            return (obj.real, obj.imag)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, fits.Header):
            return head2dict(obj)
        elif dataclasses.is_dataclass(obj):
            return dataclasses.asdict(obj)
        else:
            return super(MyEncoder, self).default(obj)


def cpu_to_use(max_cpu: int, count: int) -> int:
    """Find number of cpus to use.

    Find the right number of cpus to use when dividing up a task, such
    that there are no remainders.

    Args:
        max_cpu (int): Maximum number of cores to use for a process.
        count (int): Number of tasks.

    Returns:
        Maximum number of cores to be used that divides into the number

    """
    factors = []
    for i in range(1, count + 1):
        if count % i == 0:
            factors.append(i)
    factors_arr = np.array(factors)
    return np.max(factors_arr[factors_arr <= max_cpu])


def getfreq(
    cube: str, outdir: Union[str, None] = None, filename: Union[str, None] = None
):
    """Get list of frequencies from FITS data.

    Gets the frequency list from a given cube. Can optionally save
    frequency list to disk.

    Args:
        cube (str): File to get spectral axis from.

    Kwargs:
        outdir (str): Where to save the output file. If not given, data
            will not be saved to disk.

        filename (str): Name of frequency list file. Requires 'outdir'
            to also be specified.

        verbose (bool): Whether to print messages.

    Returns:
        freq (list): Frequencies of each channel in the input cube.

    """
    with fits.open(cube, memmap=True, mode="denywrite") as hdulist:
        hdu = hdulist[0]
        data = hdu.data
    wcs = WCS(hdu)
    freq = wcs.spectral.pixel_to_world(np.arange(data.shape[0]))  # Type: u.Quantity

    # Write to file if outdir is specified
    if outdir is None:
        return freq  # Type: u.Quantity
    else:
        if outdir[-1] == "/":
            outdir = outdir[:-1]
        if filename is None:
            outfile = f"{outdir}/frequencies.txt"
        else:
            outfile = f"{outdir}/{filename}"
        logger.info(f"Saving to {outfile}")
        np.savetxt(outfile, np.array(freq))
        return freq, outfile  # Type: Tuple[u.Quantity, str]


def gettable(tabledir: str, keyword: str, verbose=True) -> Tuple[Table, str]:
    """Get a table from a directory given a keyword to glob.

    Args:
        tabledir (str): Directory.
        keyword (str): Keyword to glob for.
        verbose (bool, optional): Verbose output. Defaults to True.

    Returns:
        Tuple[Table, str]: Table and it's file location.
    """
    if tabledir[-1] == "/":
        tabledir = tabledir[:-1]
    # Glob out the necessary files
    files = glob(f"{tabledir}/*.{keyword}*.xml")  # Selvay VOTab
    filename = files[0]
    logger.info(f"Getting table data from {filename}...")

    # Get selvay data from VOTab
    table = Table.read(filename, format="votable")
    table = table.to_pandas()
    str_df = table.select_dtypes([object])
    str_df = str_df.stack().str.decode("utf-8").unstack()
    for col in str_df:
        table[col] = str_df[col]
    return table, filename


def getdata(cubedir="./", tabledir="./", mapdata=None, verbose=True):
    """Get the spectral and source-finding data.

    Args:
        cubedir: Directory containing data cubes in FITS format.
        tabledir: Directory containing Selavy results.
        mapdata: 2D FITS image which corresponds to Selavy table.

    Kwargs:
        verbose (bool): Whether to print messages.

    Returns:
        datadict (dict): Dictionary of necessary astropy tables and
            Spectral cubes.

    """
    if cubedir[-1] == "/":
        cubedir = cubedir[:-1]

    if tabledir[-1] == "/":
        tabledir = tabledir[:-1]
    # Glob out the necessary files
    # Data cubes
    icubes = glob(f"{cubedir}/image.restored.i.*contcube*linmos.fits")
    qcubes = glob(f"{cubedir}/image.restored.q.*contcube*linmos.fits")
    ucubes = glob(f"{cubedir}/image.restored.u.*contcube*linmos.fits")
    vcubes = glob(f"{cubedir}/image.restored.v.*contcube*linmos.fits")

    cubes = [icubes, qcubes, ucubes, vcubes]
    # Selavy images
    selavyfits = mapdata
    # Get selvay data from VOTab
    i_tab, voisle = gettable(tabledir, "islands", verbose=verbose)  # Selvay VOTab
    components, tablename = gettable(tabledir, "components", verbose=verbose)

    logger.info(f"Getting spectral data from: {cubes}\n")
    logger.info(f"Getting source location data from: {selavyfits}\n")

    # Read data using Spectral cube
    i_taylor = SpectralCube.read(selavyfits, mode="denywrite")
    wcs_taylor = WCS(i_taylor.header)
    i_cube = SpectralCube.read(icubes[0], mode="denywrite")
    wcs_cube = WCS(i_cube.header)
    q_cube = SpectralCube.read(qcubes[0], mode="denywrite")
    u_cube = SpectralCube.read(ucubes[0], mode="denywrite")
    if len(vcubes) != 0:
        v_cube = SpectralCube.read(vcubes[0], mode="denywrite")
    else:
        v_cube = None
    # Mask out using Stokes I == 0 -- seems to be the current fill value
    mask = ~(i_cube == 0 * u.jansky / u.beam)
    i_cube = i_cube.with_mask(mask)
    mask = ~(q_cube == 0 * u.jansky / u.beam)
    q_cube = q_cube.with_mask(mask)
    mask = ~(u_cube == 0 * u.jansky / u.beam)
    u_cube = u_cube.with_mask(mask)

    datadict = {
        "i_tab": i_tab,
        "i_tab_comp": components,
        "i_taylor": i_taylor,
        "wcs_taylor": wcs_taylor,
        "wcs_cube": wcs_cube,
        "i_cube": i_cube,
        "q_cube": q_cube,
        "u_cube": u_cube,
        "v_cube": v_cube,
        "i_file": icubes[0],
        "q_file": qcubes[0],
        "u_file": ucubes[0],
        "v_file": vcubes[0],
    }

    return datadict


class Error(OSError):
    pass


class SameFileError(Error):
    """Raised when source and destination are the same file."""


class SpecialFileError(OSError):
    """Raised when trying to do a kind of operation (e.g. copying) which is
    not supported on a special file (e.g. a named pipe)"""


class ExecError(OSError):
    """Raised when a command could not be executed"""


class ReadError(OSError):
    """Raised when an archive cannot be read"""


class RegistryError(Exception):
    """Raised when a registry operation with the archiving
    and unpacking registeries fails"""


def _samefile(src, dst):
    # Macintosh, Unix.
    if hasattr(os.path, "samefile"):
        try:
            return os.path.samefile(src, dst)
        except OSError:
            return False


def copyfile(src, dst, *, follow_symlinks=True, verbose=True):
    """Copy data from src to dst.

    If follow_symlinks is not set and src is a symbolic link, a new
    symlink will be created instead of copying the file it points to.

    """
    if _samefile(src, dst):
        raise SameFileError("{!r} and {!r} are the same file".format(src, dst))

    for fn in [src, dst]:
        try:
            st = os.stat(fn)
        except OSError:
            # File most likely does not exist
            pass
        else:
            # XXX What about other special files? (sockets, devices...)
            if stat.S_ISFIFO(st.st_mode):
                raise SpecialFileError("`%s` is a named pipe" % fn)

    if not follow_symlinks and os.path.islink(src):
        os.symlink(os.readlink(src), dst)
    else:
        with open(src, "rb") as fsrc:
            with open(dst, "wb") as fdst:
                copyfileobj(fsrc, fdst, verbose=verbose)
    return dst


def copyfileobj(fsrc, fdst, length=16 * 1024, verbose=True):
    # copied = 0
    total = os.fstat(fsrc.fileno()).st_size
    with tqdm(
        total=total, disable=(not verbose), unit_scale=True, desc="Copying file"
    ) as pbar:
        while True:
            buf = fsrc.read(length)
            if not buf:
                break
            fdst.write(buf)
            copied = len(buf)
            pbar.update(copied)
