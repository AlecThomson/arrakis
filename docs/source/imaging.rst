Imaging
-------

.. attention::

    MeasurementSets produced by the ASKAPsoft pipeline need modification before using tools like WSClean. This can be done using `FixMS <https://fixms.readthedocs.io/>`_, which is called internally by *Arrakis*.

*Arrakis* provides an interface to the `WSClean <https://wsclean.readthedocs.io/en/latest/>`_ imaging software, with convencince functions for imaging multiple ASKAP beams simultaneously. There are two main interfaces for running the imaging pipeline:

The `spice_image` CLI and API
===================================

.. attention::

   This will only run using a sequential Prefect task runner. i.e. Only one beam will be imaged at a time.
   See either the Python API below, or the `spice_process` usage further below for parallel imaging.


This can be run using:

.. code-block::

    $ spice_image -h
    usage: spice_image [-h] [--temp_dir_wsclean TEMP_DIR_WSCLEAN] [--temp_dir_images TEMP_DIR_IMAGES] [--psf_cutoff PSF_CUTOFF] [--robust ROBUST] [--nchan NCHAN] [--pols POLS] [--size SIZE] [--scale SCALE] [--mgain MGAIN] [--niter NITER] [--nmiter NMITER] [--auto_mask AUTO_MASK]
                       [--auto_threshold AUTO_THRESHOLD] [--local_rms] [--local_rms_window LOCAL_RMS_WINDOW] [--force_mask_rounds FORCE_MASK_ROUNDS] [--gridder {direct-ft,idg,wgridder,tuned-wgridder,wstacking}] [--taper TAPER] [--minuv MINUV] [--parallel PARALLEL] [--purge] [--mpi]
                       [--multiscale] [--multiscale_scale_bias MULTISCALE_SCALE_BIAS] [--multiscale_scales MULTISCALE_SCALES] [--absmem ABSMEM] [--make_residual_cubes] [--ms_glob_pattern MS_GLOB_PATTERN] [--data_column DATA_COLUMN] [--no_mf_weighting] [--skip_fix_ms]
                       [--num_beams NUM_BEAMS] [--disable_pol_local_rms] [--disable_pol_force_mask_rounds] [--hosted-wsclean HOSTED_WSCLEAN | --local_wsclean LOCAL_WSCLEAN]
                       msdir datadir


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

        Arrkis imager


    options:
      -h, --help            show this help message and exit
      --hosted-wsclean HOSTED_WSCLEAN
                            Docker or Singularity image for wsclean (default: docker://alecthomson/wsclean:latest)
      --local_wsclean LOCAL_WSCLEAN
                            Path to local wsclean Singularity image (default: None)

    imaging arguments:
      msdir                 Directory containing MS files
      --temp_dir_wsclean TEMP_DIR_WSCLEAN
                            Temporary directory for WSClean to store intermediate files (default: None)
      --temp_dir_images TEMP_DIR_IMAGES
                            Temporary directory for to store intermediate image files (default: None)
      --psf_cutoff PSF_CUTOFF
                            Cutoff for smoothing in units of arcseconds.  (default: None)
      --robust ROBUST
      --nchan NCHAN
      --pols POLS
      --size SIZE
      --scale SCALE
      --mgain MGAIN
      --niter NITER
      --nmiter NMITER
      --auto_mask AUTO_MASK
      --auto_threshold AUTO_THRESHOLD
      --local_rms
      --local_rms_window LOCAL_RMS_WINDOW
      --force_mask_rounds FORCE_MASK_ROUNDS
      --gridder {direct-ft,idg,wgridder,tuned-wgridder,wstacking}
      --taper TAPER
      --minuv MINUV
      --parallel PARALLEL
      --purge               Purge intermediate files (default: False)
      --mpi                 Use MPI (default: False)
      --multiscale          Use multiscale clean (default: False)
      --multiscale_scale_bias MULTISCALE_SCALE_BIAS
                            The multiscale scale bias term provided to wsclean.  (default: None)
      --multiscale_scales MULTISCALE_SCALES
                            The scales used in the multiscale clean.  (default: 0,2,4,8,16,32,64,128)
      --absmem ABSMEM       Absolute memory limit in GB (default: None)
      --make_residual_cubes
                            Create residual cubes as well as cubes from restored images.  (default: False)
      --ms_glob_pattern MS_GLOB_PATTERN
                            The pattern used to search for measurement sets.  (default: scienceData*_averaged_cal.leakage.ms)
      --data_column DATA_COLUMN
                            Which column in the measurement set to image.  (default: CORRECTED_DATA)
      --no_mf_weighting     Do not use multi-frequency weighting.  (default: False)
      --skip_fix_ms         Do not apply the ASKAP MS corrections from the package fixms.  (default: False)
      --num_beams NUM_BEAMS
                            Number of beams to image (default: 36)
      --disable_pol_local_rms
                            Disable local RMS for polarisation images (default: False)
      --disable_pol_force_mask_rounds
                            Disable force mask rounds for polarisation images (default: False)

    workdir arguments:
      datadir               Directory to create/find full-size images and 'cutout' directory


You may instead prefer to use the Python API, which is more flexible and allows for parallel imaging. You will need to set up your own Prefect task-runner for this. Here is a (very) minimal example:

.. code-block:: python

    from prefect.task_runners import SequentialTaskRunner
    from arrakis.imager import main as imager_flow


    def main():
        task_runner = SequentialTaskRunner()
        imager_flow.with_options(task_runner=task_runner)(...)  # Add your arguments here


You can find the full list of arguments in the API docs here: :py:mod:`arrakis.imager.main`.


The `spice_process` CLI
=====================================

It is also possible to run just the imaging part of the pipeline using a the `spice_process` command line tool, as described in :ref:`Running the pipeline`. You will need to invoke the argument `--imager_only`, along with the other imaging arguments. This will run the imaging pipeline in parallel, using the Dask task runner defined in your config file of choice. Here is an example pipeline config for only imaging:

.. code-block:: yaml

    # SB8593.yaml
    imager_only: true
    ms_glob_pattern: 'scienceData_SB8593_RACS_1347-37A.beam*_averaged_cal.leakage.split.ms'
    imager_dask_config: petrichor.yaml
    mgain: 0.7
    force_mask_rounds: 8
    nmiter: 15
    niter: 500000
    local_rms: true
    auto_mask: 4
    local_rms_window: 60
    auto_threshold: 1
    size: 6144
    scale: 2.5
    robust: -0.5
    pols: IQU
    gridder: wgridder
    minuv: 200
    local_wsclean: wsclean_force_mask.sif
    multiscale: true
    multiscale_scale_bias: 0.7
    multiscale_scales: "0,2,4,8,16,32,64,128"
    purge: false
    absmem: 100
    nchan: 36
    psf_cutoff: 30
    skip_fix_ms: false
    data_column: CORRECTED_DATA
    disable_pol_local_rms: true
    disable_pol_force_mask_rounds: false
    temp_dir_images: /dev/shm
    temp_dir_wsclean: /dev/shm

You would then run the pipeline using:

.. code-block:: bash

    spice_process \
        --config SB8593.yaml \
        /path/to/ms/files/ \
        /path/to/work/dir/ \
        RACS_1347-37A
