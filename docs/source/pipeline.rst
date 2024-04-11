Running the pipeline
--------------------
So you're ready to run the pipeline? Make sure you've completed the :ref:`installation` and :ref:`Getting started` steps first.

The Arrakis pipeline requires 36 calibrated MeasurementSets, one per ASKAP beam. You can obtain these from the Observatory (via `CASDA <https://research.csiro.au/casda/>`_) or produce them yourself with a pipline like `Flint <https://github.com/tjgalvin/flint>`_. You'll need to have the visibilities stored in a single 'working' directory.

:code:`spice_process` and :code:`spice_field` orchestrate the pipeline flow using `Prefect <https://prefect.io>`_ and `Dask <https://dask.org>`_. These script calls the other :code:`arrakis` modules to do the work. You can control which modules run in the configuration of :code:`spice_process` or :code:`spice_field`. :code:`spice_process` operates on the level of a single RACS fields, whereas :code:`spice_field` merges multiple fields togther. You will need to run :code:`spice_process` on at least two fields before calling :code:`spice_field`. After running :code:`spice_process` or :code:`spice_field` you can run :code:`spice_cat` to produce a just a catalogue from the database values.

Details of each module can be found in the API documentation. But broadly the stages are:
    * Imaging - Create image cubes from visibilities using `WSClean <https://wsclean.readthedocs.io/>`_. This will also convolve the cubes to a common spatial resolution.

    * Cutout - Finds the position of the source in the image cubes and cuts out a square region around it.

    * LINMOS - Applies the primary beam and leakage correction to the cutout beam cubes, and then mosaics each into a single cube for each source per field.

    * FRion - Applies time-independent ionospheric Faraday rotation to the mosaicked cubes using `FRion <https://frion.readthedocs.io/en/latest/index.html/>`_.

    * RM synthesis - Extracts 1D spectra for each component of each source and runs RM synthesis using `RM-tools <https://github.com/CIRADA-Tools/RM-Tools>`_.

    * RM-CLEAN - Runs RM-CLEAN on the extracted 1D spectra using `RM-tools <https://github.com/CIRADA-Tools/RM-Tools>`_.

    * Catalogue - Queries the database for a given field and constructs a polarisation catalogue for each component.

    * Clean up - Create a tarball of the the cutouts, and remove beam cubes.

.. rst-class::  clear-both

----

With an initalised database you can call the pipeline on a single field: ::

    (arrakis310) $ spice_process -h
    usage: spice_process [-h] [--dask_config DASK_CONFIG] [--imager_dask_config IMAGER_DASK_CONFIG] [--imager_only] [--skip_imager] [--skip_cutout] [--skip_linmos]
                         [--skip_frion] [--skip_rmsynth] [--skip_rmclean] [--skip_cat] [--skip_cleanup] [--sbid SBID] [-s STOKESLIST [STOKESLIST ...]] [-e EPOCH] [-v]
                         [--host host] [--username USERNAME] [--password PASSWORD] [--limit LIMIT] [--database] [--temp_dir_wsclean TEMP_DIR_WSCLEAN]
                         [--temp_dir_images TEMP_DIR_IMAGES] [--psf_cutoff PSF_CUTOFF] [--robust ROBUST] [--nchan NCHAN] [--pols POLS] [--size SIZE] [--scale SCALE]
                         [--mgain MGAIN] [--niter NITER] [--nmiter NMITER] [--auto_mask AUTO_MASK] [--auto_threshold AUTO_THRESHOLD] [--local_rms]
                         [--local_rms_window LOCAL_RMS_WINDOW] [--force_mask_rounds FORCE_MASK_ROUNDS] [--gridder {direct-ft,idg,wgridder,tuned-wgridder,wstacking}]
                         [--taper TAPER] [--minuv MINUV] [--parallel PARALLEL] [--purge] [--mpi] [--multiscale] [--multiscale_scale_bias MULTISCALE_SCALE_BIAS]
                         [--multiscale_scales MULTISCALE_SCALES] [--absmem ABSMEM] [--make_residual_cubes] [--ms_glob_pattern MS_GLOB_PATTERN] [--data_column DATA_COLUMN]
                         [--no_mf_weighting] [--skip_fix_ms] [--hosted-wsclean HOSTED_WSCLEAN | --local_wsclean LOCAL_WSCLEAN] [-p PAD] [-d] [--holofile HOLOFILE]
                         [--yanda YANDA] [--yanda_image YANDA_IMAGE] [--ionex_server IONEX_SERVER] [--ionex_prefix IONEX_PREFIX] [--ionex_formatter IONEX_FORMATTER]
                         [--ionex_proxy_server IONEX_PROXY_SERVER] [--ionex_predownload] [--dimension DIMENSION] [--save_plots] [--rm_verbose] [--ion] [--tt0 TT0]
                         [--tt1 TT1] [--validate] [--own_fit] [--weight_type WEIGHT_TYPE] [--fit_function FIT_FUNCTION] [--fit_rmsf] [--phi_max PHI_MAX] [--dphi DPHI]
                         [--n_samples N_SAMPLES] [--poly_ord POLY_ORD] [--no_stokes_i] [--show_plots] [--not_rmsf] [--debug] [--cutoff CUTOFF] [--max_iter MAX_ITER]
                         [--gain GAIN] [--window WINDOW] [--leakage_degree LEAKAGE_DEGREE] [--leakage_bins LEAKAGE_BINS] [--leakage_snr LEAKAGE_SNR] [--write OUTFILE]
                         [--overwrite] [--config CONFIG]
                         datadir field msdir


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

        Arrakis pipeline.

        Before running make sure to start a session of mongodb e.g.
            $ mongod --dbpath=/path/to/database --bind_ip $(hostname -i)



    options:
      -h, --help            show this help message and exit
      --hosted-wsclean HOSTED_WSCLEAN
                            Docker or Singularity image for wsclean (default: docker://alecthomson/wsclean:latest)
      --local_wsclean LOCAL_WSCLEAN
                            Path to local wsclean Singularity image (default: None)
      --config CONFIG       Config file path (default: None)

    pipeline arguments:
      --dask_config DASK_CONFIG
                            Config file for Dask SlurmCLUSTER. (default: None)
      --imager_dask_config IMAGER_DASK_CONFIG
                            Config file for Dask SlurmCLUSTER. (default: None)
      --imager_only         Only run the imager component of the pipeline.  (default: False)
      --skip_imager         Skip imaging stage [False]. (default: False)
      --skip_cutout         Skip cutout stage [False]. (default: False)
      --skip_linmos         Skip LINMOS stage [False]. (default: False)
      --skip_frion          Skip cleanup stage [False]. (default: False)
      --skip_rmsynth        Skip RM Synthesis stage [False]. (default: False)
      --skip_rmclean        Skip RM-CLEAN stage [False]. (default: False)
      --skip_cat            Skip catalogue stage [False]. (default: False)
      --skip_cleanup        Skip cleanup stage [False]. (default: False)

    workdir arguments:
      datadir               Directory to create/find full-size images and 'cutout' directory

    generic arguments:
      field                 Name of field (e.g. RACS_2132-50).
      --sbid SBID           SBID of observation. (default: None)
      -s STOKESLIST [STOKESLIST ...], --stokes STOKESLIST [STOKESLIST ...]
                            List of Stokes parameters to image (default: ['I', 'Q', 'U'])
      -e EPOCH, --epoch EPOCH
                            Epoch of observation. (default: 0)
      -v                    Verbose output. (default: False)
      --host host           Host of mongodb (probably $hostname -i). (default: None)
      --username USERNAME   Username of mongodb. (default: None)
      --password PASSWORD   Password of mongodb. (default: None)
      --limit LIMIT         Limit the number of islands to process. (default: None)
      --database            Add data to MongoDB. (default: False)

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

    cutout arguments:
      -p PAD, --pad PAD     Number of beamwidths to pad around source [3]. (default: 3)
      -d, --dryrun          Do a dry-run [False]. (default: False)

    linmos arguments:
      --holofile HOLOFILE   Path to holography image (default: None)
      --yanda YANDA         Yandasoft version to pull from DockerHub [1.3.0]. (default: 1.3.0)
      --yanda_image YANDA_IMAGE
                            Path to an existing yandasoft singularity container image.  (default: None)

    frion arguments:
      --ionex_server IONEX_SERVER
                            IONEX server (default: ftp://ftp.aiub.unibe.ch/CODE/)
      --ionex_prefix IONEX_PREFIX
      --ionex_formatter IONEX_FORMATTER
                            IONEX formatter. (default: ftp.aiub.unibe.ch)
      --ionex_proxy_server IONEX_PROXY_SERVER
                            Proxy server. (default: None)
      --ionex_predownload   Pre-download IONEX files. (default: False)

    common rm arguments:
      --dimension DIMENSION
                            How many dimensions for RMsynth '1d' or '3d'. (default: 1d)
      --save_plots          save the plots. (default: False)
      --rm_verbose          Verbose RMsynth/RMClean. (default: False)

    rm-synth arguments:
      --ion                 Use ionospheric-corrected data. (default: False)
      --tt0 TT0             TT0 MFS image -- will be used for model of Stokes I -- also needs --tt1. (default: None)
      --tt1 TT1             TT1 MFS image -- will be used for model of Stokes I -- also needs --tt0. (default: None)
      --validate            Run on Stokes I. (default: False)
      --own_fit             Use own Stokes I fit function. (default: False)
      --weight_type WEIGHT_TYPE
                            weighting (inverse) 'variance' or 'uniform' (all 1s). (default: variance)
      --fit_function FIT_FUNCTION
                            Stokes I fitting function: 'linear' or 'log' polynomials. (default: log)
      --fit_rmsf            Fit a Gaussian to the RMSF (default: False)
      --phi_max PHI_MAX     Absolute max Faraday depth sampled (in rad/m^2) (overrides NSAMPLES). (default: None)
      --dphi DPHI           Width of Faraday depth channel. (default: None)
      --n_samples N_SAMPLES
                            Number of samples across the FWHM RMSF. (default: 5)
      --poly_ord POLY_ORD   polynomial order to fit to I spectrum. (default: 3)
      --no_stokes_i         ignore the Stokes I spectrum. (default: False)
      --show_plots          show the plots. (default: False)
      --not_rmsf            Skip calculation of RMSF? (default: False)
      --debug               turn on debugging messages & plots. (default: False)

    rm-clean arguments:
      --cutoff CUTOFF       CLEAN cutoff (+ve = absolute, -ve = sigma). (default: -3)
      --max_iter MAX_ITER   maximum number of CLEAN iterations. (default: 10000)
      --gain GAIN           CLEAN loop gain. (default: 0.1)
      --window WINDOW       Further CLEAN in mask to this threshold. (default: None)

    catalogue arguments:
      --leakage_degree LEAKAGE_DEGREE
                            Degree of leakage polynomial fit. (default: 4)
      --leakage_bins LEAKAGE_BINS
                            Number of bins for leakage fit. (default: 16)
      --leakage_snr LEAKAGE_SNR
                            SNR cut for leakage fit. (default: 30.0)
      --write OUTFILE       File to save table to. (default: None)

    cleanup arguments:
      --overwrite           Overwrite existing tarball (default: False)

    Args that start with '--' can also be set in a config file (/scratch3/projects/spiceracs/arrakis/arrakis/.default_config.yaml or specified via --config). Config file
    syntax allows: key=value, flag=true, stuff=[a,b,c] (for details, see syntax at https://goo.gl/R74nmi). In general, command-line values override config file values which
    override defaults.


You can optionally pass a configuration file (with the :code:`--config` argument) to set the options you prefer. An example file in contained in :file:`arrakis/.default_config.yaml`:

.. code-block:: yaml

  # options:
  hosted-wsclean: docker://alecthomson/wsclean:latest # Docker or Singularity image for wsclean (default: docker://alecthomson/wsclean:latest)
  local_wsclean: null # Path to local wsclean Singularity image (default: None)

  # pipeline arguments:
  dask_config: null # Config file for Dask SlurmCLUSTER. (default: None)
  imager_dask_config: null #Config  file for Dask SlurmCLUSTER. (default: None)
  imager_only: false # Only run the imager component of the pipeline.  (default: False)
  skip_imager: false #Skip imaging stage [False]. (default: False)
  skip_cutout: false #Skip cutout stage [False]. (default: False)
  skip_linmos: false #Skip LINMOS stage [False]. (default: False)
  skip_frion: false #Skip cleanup stage [False]. (default: False)
  skip_rmsynth: false #Skip RM Synthesis stage [False]. (default: False)
  skip_rmclean: false #Skip RM-CLEAN stage [False]. (default: False)
  skip_cat: false #Skip catalogue stage [False]. (default: False)
  skip_cleanup: false #Skip cleanup stage [False]. (default: False)

  # generic null arguments:
  sbid: null #SBID of observation. (default: None)
  stokes: # List of Stokes parameters to image (default: ['I', 'Q', 'U'])
    - I
    - Q
    - U
  epoch: 0 # Epoch of observation. (default: 0)
  host: null # Host of mongodb (probably $hostname -i). (default: None)
  username: null # Username of mongodb. (default: None)
  password: # Password of mongodb. (default: None)
  limit: null # Limit the number of islands to process. (default: None)
  database: false # Add data to MongoDB. (default: False)

  # imaging arguments:
  temp_dir_wsclean: null # Temporary directory for WSClean to store intermediate files (default: None)
  temp_dir_images: null # Temporary directory for to store intermediate image files (default: None)
  psf_cutoff: null # Cutoff for smoothing in units of arcseconds.  (default: None)
  robust: -0.5 # ROBUST
  nchan: 36 # NCHAN
  pols: IQU # POLS
  size: 6144 # SIZE
  scale: 2.5 # SCALE
  mgain: 0.7 # MGAIN
  niter: 500_000 # NITER
  nmiter: 15 # NMITER
  auto_mask: 4 # AUTO_MASK
  auto_threshold: 1 # AUTO_THRESHOLD
  local_rms: true #
  local_rms_window: 60 # LOCAL_RMS_WINDOW
  force_mask_rounds: 8 # FORCE_MASK_ROUNDS
  gridder: wgridder # {direct-ft,idg,wgridder,tuned-wgridder,wstacking}
  taper: null # TAPER
  minuv: 200 # MINUV
  parallel: null # PARALLEL
  mpi: false #                 Use MPI (default: False)
  purge: false # Purge intermediate files (default: False)
  multiscale: false # Use multiscale clean (default: False)
  multiscale_scale_bias: null # The multiscale scale bias term provided to wsclean.  (default: None)
  multiscale_scales: 0,2,4,8,16,32,64,12 # The scales used in the multiscale clean.  (default: 0,2,4,8,16,32,64,128)
  absmem: null # ABSMEM       Absolute memory limit in GB (default: None)
  make_residual_cubes: false # Create residual cubes as well as cubes from restored images.  (default: False)
  ms_glob_pattern: scienceData*_averaged_cal.leakage.ms # The pattern used to search for measurement sets.  (default: scienceData*_averaged_cal.leakage.ms)
  data_column: CORRECTED_DATA # Which column in the measurement set to image.  (default: CORRECTED_DATA)
  no_mf_weighting: false # Do not use multi-frequency weighting.  (default: False)
  skip_fix_ms: false # Do not apply the ASKAP MS corrections from the package fixms.  (default: False)

  # cutout arguments:
  pad: 3 # Number of beamwidths to pad around source [3]. (default: 3)
  dryrun: false # Do a dry-run [False]. (default: False)

  # linmos null arguments:
  holofile: null #Path to holography image (default: None)
  yanda: 1.3.0 # Yandasoft version to pull from DockerHub [1.3.0]. (default: 1.3.0)
  yanda_image: null #Path to an existing yandasoft singularity container image.  (default: None)

  # frion arguments:
  ionex_server: ftp://ftp.aiub.unibe.ch/CODE/ # IONEX server (default: ftp://ftp.aiub.unibe.ch/CODE/)
  ionex_prefix: codg # IONEX_PREFIX
  ionex_formatter: null # IONEX formatter. (default: ftp.aiub.unibe.ch)
  ionex_proxy_server: null # Proxy server. (default: None)
  ionex_predownload: false # Pre-download IONEX files. (default: False)

  # common rm arguments:
  dimension: 1d # How many dimensions for RMsynth '1d' or '3d'. (default: 1d)
  save_plots: false #          save the plots. (default: False)
  rm_verbose: false #          Verbose RMsynth/RMClean. (default: False)

  # rm-synth arguments:
  ion: false # Use ionospheric-corrected data. (default: False)
  tt0: null # TT0 MFS image -- will be used for model of Stokes I -- also needs --tt1. (default: None)
  tt1: null # TT1 MFS image -- will be used for model of Stokes I -- also needs --tt0. (default: None)
  validate: false # Run on Stokes I. (default: False)
  own_fit: false # Use own Stokes I fit function. (default: False)
  weight_type: # weighting (inverse) 'variance' or 'uniform' (all 1s). (default: variance)
  fit_function: # Stokes I fitting function: 'linear' or 'log' polynomials. (default: log)
  fit_rmsf: false # Fit a Gaussian to the RMSF (default: False)
  phi_max: null # Absolute max Faraday depth sampled (in rad/m^2) (overrides NSAMPLES). (default: None)
  dphi: null # Width of Faraday depth channel. (default: None)
  n_samples: # Number of samples across the FWHM RMSF. (default: 5)
  poly_ord: # polynomial order to fit to I spectrum. (default: 3)
  no_stokes_i: false # ignore the Stokes I spectrum. (default: False)
  show_plots: false # show the plots. (default: False)
  not_rmsf: false # Skip calculation of RMSF? (default: False)
  debug: false # turn on debugging messages & plots. (default: False)

  # rm-clean arguments:
  cutoff: -8 # CLEAN cutoff (+ve = absolute, -ve = sigma). (default: -3)
  max_iter: 10000 # maximum number of CLEAN iterations. (default: 10000)
  gain: 0.1 # CLEAN loop gain. (default: 0.1)
  window: null # Further CLEAN in mask to this threshold. (default: None)

  # catalogue arguments:
  leakage_degree: 4 # Degree of leakage polynomial fit. (default: 4)
  leakage_bins: 16 # Number of bins for leakage fit. (default: 16)
  leakage_snr: 30 # SNR cut for leakage fit. (default: 30.0)
  write: null # File to save table to. (default: None)

  # cleanup arguments:
  overwrite: false # Overwrite existing tarball (default: False)



For extra information you can refer to the API:

* :py:mod:`arrakis.process_spice`

Similarly, you can merge multiple fields togther using: ::

    (arrakis310) $ spice_region -h
    usage: spice_region [-h] [--dask_config DASK_CONFIG] [--skip_frion] [--skip_rmsynth] [--skip_rmclean] [--skip_cat] [--skip_cleanup] [--merge_name MERGE_NAME]
                        [--fields FIELDS [FIELDS ...]] [--datadirs DATADIRS [DATADIRS ...]] [--output_dir OUTPUT_DIR] [-e EPOCH] [--host host] [--username USERNAME]
                        [--password PASSWORD] [--holofile HOLOFILE] [--yanda YANDA] [--yanda_image YANDA_IMAGE] [--dimension DIMENSION] [--save_plots] [--rm_verbose]
                        [--ion] [--tt0 TT0] [--tt1 TT1] [--validate] [--own_fit] [--weight_type WEIGHT_TYPE] [--fit_function FIT_FUNCTION] [--fit_rmsf] [--phi_max PHI_MAX]
                        [--dphi DPHI] [--n_samples N_SAMPLES] [--poly_ord POLY_ORD] [--no_stokes_i] [--show_plots] [--not_rmsf] [--debug] [--cutoff CUTOFF]
                        [--max_iter MAX_ITER] [--gain GAIN] [--window WINDOW] [--leakage_degree LEAKAGE_DEGREE] [--leakage_bins LEAKAGE_BINS] [--leakage_snr LEAKAGE_SNR]
                        [--write OUTFILE] [--overwrite] [--config CONFIG]


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

        Arrakis regional pipeline.

        Before running make sure to start a session of mongodb e.g.
            $ mongod --dbpath=/path/to/database --bind_ip $(hostname -i)



    options:
      -h, --help            show this help message and exit
      --config CONFIG       Config file path (default: None)

    pipeline arguments:
      --dask_config DASK_CONFIG
                            Config file for Dask SlurmCLUSTER. (default: None)
      --skip_frion          Skip cleanup stage [False]. (default: False)
      --skip_rmsynth        Skip RM Synthesis stage [False]. (default: False)
      --skip_rmclean        Skip RM-CLEAN stage [False]. (default: False)
      --skip_cat            Skip catalogue stage [False]. (default: False)
      --skip_cleanup        Skip cleanup stage [False]. (default: False)

    merge arguments:
      --merge_name MERGE_NAME
                            Name of the merged region (default: None)
      --fields FIELDS [FIELDS ...]
                            RACS fields to mosaic - e.g. RACS_2132-50A. (default: None)
      --datadirs DATADIRS [DATADIRS ...]
                            Directories containing cutouts (in subdir outdir/cutouts).. (default: None)
      --output_dir OUTPUT_DIR
                            Path to save merged data (in output_dir/merge_name/cutouts) (default: None)
      -e EPOCH, --epoch EPOCH
                            Epoch of observation. (default: 0)
      --host host           Host of mongodb (probably $hostname -i). (default: None)
      --username USERNAME   Username of mongodb. (default: None)
      --password PASSWORD   Password of mongodb. (default: None)

    linmos arguments:
      --holofile HOLOFILE   Path to holography image (default: None)
      --yanda YANDA         Yandasoft version to pull from DockerHub [1.3.0]. (default: 1.3.0)
      --yanda_image YANDA_IMAGE
                            Path to an existing yandasoft singularity container image.  (default: None)

    common rm arguments:
      --dimension DIMENSION
                            How many dimensions for RMsynth '1d' or '3d'. (default: 1d)
      --save_plots          save the plots. (default: False)
      --rm_verbose          Verbose RMsynth/RMClean. (default: False)

    rm-synth arguments:
      --ion                 Use ionospheric-corrected data. (default: False)
      --tt0 TT0             TT0 MFS image -- will be used for model of Stokes I -- also needs --tt1. (default: None)
      --tt1 TT1             TT1 MFS image -- will be used for model of Stokes I -- also needs --tt0. (default: None)
      --validate            Run on Stokes I. (default: False)
      --own_fit             Use own Stokes I fit function. (default: False)
      --weight_type WEIGHT_TYPE
                            weighting (inverse) 'variance' or 'uniform' (all 1s). (default: variance)
      --fit_function FIT_FUNCTION
                            Stokes I fitting function: 'linear' or 'log' polynomials. (default: log)
      --fit_rmsf            Fit a Gaussian to the RMSF (default: False)
      --phi_max PHI_MAX     Absolute max Faraday depth sampled (in rad/m^2) (overrides NSAMPLES). (default: None)
      --dphi DPHI           Width of Faraday depth channel. (default: None)
      --n_samples N_SAMPLES
                            Number of samples across the FWHM RMSF. (default: 5)
      --poly_ord POLY_ORD   polynomial order to fit to I spectrum. (default: 3)
      --no_stokes_i         ignore the Stokes I spectrum. (default: False)
      --show_plots          show the plots. (default: False)
      --not_rmsf            Skip calculation of RMSF? (default: False)
      --debug               turn on debugging messages & plots. (default: False)

    rm-clean arguments:
      --cutoff CUTOFF       CLEAN cutoff (+ve = absolute, -ve = sigma). (default: -3)
      --max_iter MAX_ITER   maximum number of CLEAN iterations. (default: 10000)
      --gain GAIN           CLEAN loop gain. (default: 0.1)
      --window WINDOW       Further CLEAN in mask to this threshold. (default: None)

    catalogue arguments:
      --leakage_degree LEAKAGE_DEGREE
                            Degree of leakage polynomial fit. (default: 4)
      --leakage_bins LEAKAGE_BINS
                            Number of bins for leakage fit. (default: 16)
      --leakage_snr LEAKAGE_SNR
                            SNR cut for leakage fit. (default: 30.0)
      --write OUTFILE       File to save table to. (default: None)

    cleanup arguments:
      --overwrite           Overwrite existing tarball (default: False)

    Args that start with '--' can also be set in a config file (.default_config.cfg or specified via --config). Config file syntax allows: key=value, flag=true,
    stuff=[a,b,c] (for details, see syntax at https://goo.gl/R74nmi). In general, command-line values override config file values which override defaults.


* :py:mod:`arrakis.process_region`

Helper scripts (mostly for bespoke purposes) are available on the commandline. See the API reference for more details.
