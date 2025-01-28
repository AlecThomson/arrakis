#!/usr/bin/env python
"""MeasurementSet utilities"""

import copy
import warnings
from pathlib import Path
from typing import Optional

import astropy.units as u
from astropy.utils.exceptions import AstropyWarning
from casacore.tables import table
from spectral_cube.utils import SpectralCubeWarning

from arrakis.logger import logger

warnings.filterwarnings(action="ignore", category=SpectralCubeWarning, append=True)
warnings.simplefilter("ignore", category=AstropyWarning)


def get_pol_axis(
    ms: Path, feed_idx: Optional[int] = None, col: str = "RECEPTOR_ANGLE"
) -> u.Quantity:
    """Get the polarization axis from the ASKAP MS. Checks are performed
    to ensure this polarisation axis angle is constant throughout the observation.


    Args:
        ms (Path): The path to the measurement set that will be inspected
        feed_idx (Optional[int], optional): Specify the entery in the FEED
        table of `ms` to return. This might be required when a subset of a
        measurement set has been extracted from an observation with a varying
        orientation.
        col (str, optional): The column to extract the polarization angle from.

    Returns:
        astropy.units.Quantity: The rotation of the PAF throughout the observing.
    """
    _known_cols = ("RECEPTOR_ANGLE", "INSTRUMENT_RECEPTOR_ANGLE")
    if col not in _known_cols:
        raise ValueError(f"Unknown column {col=}, please use one of {_known_cols}")
    with table((ms / "FEED").as_posix(), readonly=True, ack=False) as tf:
        ms_feed = tf.getcol(col) * u.rad
        # PAF is at 45deg to feeds
        # 45 - feed_angle = pol_angle
        pol_axes = -(ms_feed - 45.0 * u.deg)

    if feed_idx is None:
        assert (ms_feed[:, 0] == ms_feed[0, 0]).all() & (
            ms_feed[:, 1] == ms_feed[0, 1]
        ).all(), f"The {col} changes with time, please check the MS"

        pol_ang = pol_axes[0, 0].to(u.deg)
    else:
        logger.debug(f"Extracting the third-axis orientation for {feed_idx=}")
        pol_ang = pol_axes[feed_idx, 0].to(u.deg)

    return pol_ang


def beam_from_ms(ms: str) -> int:
    """Work out which beam is in this MS"""
    with table(ms, readonly=True, ack=False) as t:
        vis_feed = t.getcol("FEED1", 0, 1)
        beam = vis_feed[0]
    return beam


def field_idx_from_ms(ms: str) -> int:
    """Get the field from MS metadata"""
    with table(f"{ms}/FIELD", readonly=True, ack=False) as field:
        idxs = list(field.SOURCE_ID)
        assert len(idxs) == 1 or all([idx == idxs[0] for idx in idxs]), (
            "More than one field in MS"
        )
        idx = idxs[0]
    return idx


def field_name_from_ms(ms: str) -> str:
    """Get the field name from MS metadata"""
    with table(f"{ms}/FIELD", readonly=True, ack=False) as field:
        names = list(field.NAME)
        assert len(names) == 1, "More than one field in MS"
        name = names[0]
    return name


def wsclean(
    mslist: list,
    use_mpi: bool,
    version: bool = False,
    j: Optional[int] = None,
    parallel_gridding: Optional[int] = None,
    parallel_reordering: Optional[int] = None,
    no_work_on_master: bool = False,
    mem: Optional[float] = None,
    abs_mem: Optional[float] = None,
    verbose: bool = False,
    log_time: bool = False,
    quiet: bool = False,
    reorder: bool = False,
    no_reorder: bool = False,
    temp_dir: Optional[str] = None,
    update_model_required: bool = False,
    no_update_model_required: bool = False,
    no_dirty: bool = False,
    save_first_residual: bool = False,
    save_weights: bool = False,
    save_uv: bool = False,
    reuse_psf: Optional[str] = None,
    reuse_dirty: Optional[str] = None,
    apply_primary_beam: bool = False,
    reuse_primary_beam: bool = False,
    use_differential_lofar_beam: bool = False,
    primary_beam_limit: Optional[float] = None,
    mwa_path: Optional[str] = None,
    save_psf_pb: bool = False,
    pb_grid_size: Optional[int] = None,
    beam_model: Optional[str] = None,
    beam_mode: Optional[str] = None,
    beam_normalisation_mode: Optional[str] = None,
    dry_run: bool = False,
    weight: Optional[str] = None,
    super_weight: Optional[float] = None,
    mf_weighting: bool = False,
    no_mf_weighting: bool = False,
    weighting_rank_filter: Optional[float] = None,
    weighting_rank_filter_size: Optional[float] = None,
    taper_gaussian: Optional[str] = None,
    taper_tukey: Optional[float] = None,
    taper_inner_tukey: Optional[float] = None,
    taper_edge: Optional[float] = None,
    taper_edge_tukey: Optional[float] = None,
    use_weights_as_taper: bool = False,
    store_imaging_weights: bool = False,
    name: Optional[str] = None,
    size: Optional[str] = None,
    padding: Optional[float] = None,
    scale: Optional[str] = None,
    predict: bool = False,
    ws_continue: bool = False,
    subtract_model: bool = False,
    gridder: Optional[str] = None,
    channels_out: Optional[int] = None,
    shift: Optional[str] = None,
    gap_channel_division: bool = False,
    channel_division_frequencies: Optional[str] = None,
    nwlayers: Optional[int] = None,
    nwlayers_factor: Optional[float] = None,
    nwlayers_for_size: Optional[str] = None,
    no_small_inversion: bool = False,
    small_inversion: bool = False,
    grid_mode: Optional[str] = None,
    kernel_size: Optional[int] = None,
    oversampling: Optional[int] = None,
    make_psf: bool = False,
    make_psf_only: bool = False,
    visibility_weighting_mode: Optional[str] = None,
    no_normalize_for_weighting: bool = False,
    baseline_averaging: Optional[float] = None,
    simulate_noise: Optional[float] = None,
    simulate_baseline_noise: Optional[str] = None,
    idg_mode: Optional[str] = None,
    wgridder_accuracy: Optional[float] = None,
    aterm_config: Optional[str] = None,
    grid_with_beam: bool = False,
    beam_aterm_update: Optional[int] = False,
    aterm_kernel_size: Optional[float] = None,
    apply_facet_solutions: Optional[str] = None,
    apply_facet_beam: bool = False,
    facet_beam_update: Optional[int] = False,
    save_aterms: bool = False,
    pol: Optional[str] = None,
    interval: Optional[str] = None,
    intervals_out: Optional[int] = None,
    even_timesteps: bool = False,
    odd_timesteps: bool = False,
    channel_range: Optional[str] = None,
    field: Optional[int] = None,
    spws: Optional[str] = None,
    data_column: Optional[str] = None,
    maxuvw_m: Optional[float] = None,
    minuvw_m: Optional[float] = None,
    maxuv_l: Optional[float] = None,
    minuv_l: Optional[float] = None,
    maxw: Optional[float] = None,
    niter: Optional[int] = None,
    nmiter: Optional[int] = None,
    threshold: Optional[float] = None,
    auto_threshold: Optional[float] = None,
    auto_mask: Optional[float] = None,
    force_mask_rounds: Optional[int] = None,
    local_rms: bool = False,
    local_rms_window: Optional[float] = False,
    local_rms_method: Optional[str] = None,
    gain: Optional[float] = None,
    mgain: Optional[float] = None,
    join_polarizations: bool = False,
    link_polarizations: Optional[str] = None,
    facet_regions: Optional[str] = None,
    join_channels: bool = False,
    spectral_correction: Optional[str] = None,
    no_fast_subminor: bool = False,
    multiscale: bool = False,
    multiscale_scale_bias: Optional[float] = None,
    multiscale_max_scales: Optional[int] = None,
    multiscale_scales: Optional[str] = None,
    multiscale_shape: Optional[str] = None,
    multiscale_gain: Optional[float] = None,
    multiscale_convolution_padding: Optional[float] = None,
    no_multiscale_fast_subminor: bool = False,
    python_deconvolution: Optional[str] = None,
    iuwt: bool = False,
    iuwt_snr_test: bool = False,
    no_iuwt_snr_test: bool = False,
    moresane_ext: Optional[str] = None,
    moresane_arg: Optional[str] = None,
    moresane_sl: Optional[str] = None,
    save_source_list: bool = False,
    clean_border: Optional[float] = None,
    fits_mask: Optional[str] = None,
    casa_mask: Optional[str] = None,
    horizon_mask: Optional[str] = None,
    no_negative: bool = False,
    negative: bool = False,
    stop_negative: bool = False,
    fit_spectral_pol: Optional[int] = None,
    fit_spectral_log_pol: Optional[int] = None,
    force_spectrum: Optional[str] = None,
    deconvolution_channels: Optional[int] = None,
    squared_channel_joining: bool = False,
    parallel_deconvolution: Optional[int] = None,
    deconvolution_threads: Optional[int] = None,
    restore: Optional[str] = None,
    restore_list: Optional[str] = None,
    beam_size: Optional[float] = None,
    beam_shape: Optional[str] = None,
    fit_beam: bool = False,
    no_fit_beam: bool = False,
    beam_fitting_size: Optional[float] = None,
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
        logger.warning("CAUTION - square channel joining and multiscale is unstable!")

    for key, value in arguments.items():
        if isinstance(value, bool):
            if value:
                command += f" -{key.replace('_', '-')}"
        elif value:
            if "ws_" in key:  # Catch for ws_continue command
                key.replace("ws_", "")
            command += f" -{key.replace('_', '-')} {value}"
    command += f" {' '.join(mslist)}"
    return command
