Running the pipeline
--------------------
So you're ready to run the pipeline? Make sure you've completed the :ref:`installation` and :ref:`Getting started` steps first.

The Arrakis pipeline requires 36 calibrated MeasurementSets, one per ASKAP beam. You can obtain these from the Observatory (via `CASDA <https://research.csiro.au/casda/>`_) or produce them yourself with a pipline like `Flint <https://github.com/tjgalvin/flint>`_. You'll need to have the visibilities stored in a single 'working' directory.

:code:`spice_process` and :code:`spice_field` orchestrate the pipeline flow using `Prefect <https://prefect.io>`_ and `Dask <https://dask.org>`_. These script calls the other :code:`arrakis` modules to do the work. You can control which modules run in the configuration of :code:`spice_process` or :code:`spice_field`. :code:`spice_process` operates on the level of a single RACS fields, whereas :code:`spice_field` merges multiple fields togther. You will need to run :code:`spice_process` on at least two fields before calling :code:`spice_field`. After running :code:`spice_process` or :code:`spice_field` you can run :code:`spice_cat` to produce a just a catalogue from the database values.

Details of each module can be found in the API documentation. But broadly the stages are:
    * Imaging - Create image cubes from visibilities using `WSClean <https://wsclean.readthedocs.io/>`_. This will also convolve the cubes to a common spatial resolution.

    * Cutout - Finds the position of the source in the image cubes and cuts out a square region around it.

    * LINMOS - Applies the primary beam and leakage correction to the cutout beam cubes, and then mosaics each into a single cube for each source per field.

    * Clean up - Remove the cutout beam cubes from the cutouts directory.

    * FRion - Applies time-independent ionospheric Faraday rotation to the mosaicked cubes using `FRion <https://frion.readthedocs.io/en/latest/index.html/>`_.

    * RM synthesis - Extracts 1D spectra for each component of each source and runs RM synthesis using `RM-tools <https://github.com/CIRADA-Tools/RM-Tools>`_.

    * RM-CLEAN - Runs RM-CLEAN on the extracted 1D spectra using `RM-tools <https://github.com/CIRADA-Tools/RM-Tools>`_.

    * Catalogue - Queries the database for a given field and constructs a polarisation catalogue for each component.

.. rst-class::  clear-both

----

With an initalised database you can call the pipeline on a single field: ::

    (spice) $ spice_process -h
    usage: spice_process [-h] [--psf_cutoff PSF_CUTOFF] [--robust ROBUST] [--nchan NCHAN] [--pols POLS] [--size SIZE] [--scale SCALE] [--mgain MGAIN]
                        [--niter NITER] [--nmiter NMITER] [--auto_mask AUTO_MASK] [--auto_threshold AUTO_THRESHOLD] [--local_rms]
                        [--local_rms_window LOCAL_RMS_WINDOW] [--force_mask_rounds FORCE_MASK_ROUNDS]
                        [--gridder {direct-ft,idg,wgridder,tuned-wgridder,wstacking}] [--taper TAPER] [--minuv MINUV] [--parallel PARALLEL] [--purge] [--mpi]
                        [--multiscale] [--multiscale_scale_bias MULTISCALE_SCALE_BIAS] [--absmem ABSMEM] [--make_residual_cubes]
                        [--ms_glob_pattern MS_GLOB_PATTERN] [--data_column DATA_COLUMN] [--skip_fix_ms]
                        [--hosted-wsclean HOSTED_WSCLEAN | --local_wsclean LOCAL_WSCLEAN] [--config CONFIG] [--epoch EPOCH] [--host HOST]
                        [--username USERNAME] [--password PASSWORD] [--use_mpi] [--port_forward PORT_FORWARD [PORT_FORWARD ...]] [--dask_config DASK_CONFIG]
                        [--imager_dask_config IMAGER_DASK_CONFIG] [--holofile HOLOFILE] [--yanda YANDA] [--yanda_image YANDA_IMAGE] [--imager_only]
                        [--skip_imager] [--skip_cutout] [--skip_linmos] [--skip_cleanup] [--skip_frion] [--skip_rmsynth] [--skip_rmclean] [--skip_cat] [-v]
                        [-vw] [-p PAD] [--dryrun] [--dimension DIMENSION] [-m] [--tt0 TT0] [--tt1 TT1] [--validate] [--limit LIMIT] [--own_fit] [-sp]
                        [-w WEIGHTTYPE] [--fit_function FIT_FUNCTION] [-t] [-l PHIMAX_RADM2] [-d DPHI_RADM2] [-s NSAMPLES] [-o POLYORD] [-i] [--showPlots]
                        [-R] [-rmv] [-D] [-c CUTOFF] [-n MAXITER] [-g GAIN] [--window WINDOW] [--ionex_server IONEX_SERVER] [--ionex_prefix IONEX_PREFIX]
                        [--ionex_proxy_server IONEX_PROXY_SERVER] [--ionex_formatter IONEX_FORMATTER] [--ionex_predownload] [--outfile OUTFILE]
                        msdir outdir field


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


    Arrakis pipeline. Before running make sure to start a session of mongodb e.g. $ mongod
    --dbpath=/path/to/database --bind_ip $(hostname -i)

    positional arguments:
    field                 Name of field (e.g. 2132-50A).

    optional arguments:
    -h, --help            show this help message and exit
    --hosted-wsclean HOSTED_WSCLEAN
                            Docker or Singularity image for wsclean [docker://alecthomson/wsclean:latest] (default: docker://alecthomson/wsclean:latest)
    --local_wsclean LOCAL_WSCLEAN
                            Path to local wsclean Singularity image (default: None)
    --config CONFIG       Config file path (default: None)
    --epoch EPOCH         Epoch to read field data from (default: 0)
    --host HOST           Host of mongodb (probably $hostname -i). (default: None)
    --username USERNAME   Username of mongodb. (default: None)
    --password PASSWORD   Password of mongodb. (default: None)
    --use_mpi             Use Dask-mpi to parallelise -- must use srun/mpirun to assign resources. (default: False)
    --port_forward PORT_FORWARD [PORT_FORWARD ...]
                            Platform to fowards dask port [None]. (default: None)
    --dask_config DASK_CONFIG
                            Config file for Dask SlurmCLUSTER. (default: None)
    --imager_dask_config IMAGER_DASK_CONFIG
                            Config file for Dask SlurmCLUSTER. (default: None)
    --holofile HOLOFILE   Path to holography image (default: None)
    --yanda YANDA         Yandasoft version to pull from DockerHub [1.3.0]. (default: 1.3.0)
    --yanda_image YANDA_IMAGE
                            Path to an existing yandasoft singularity container image. (default: None)

    imaging arguments:
    msdir                 Directory containing MS files
    outdir                Directory to output images
    --psf_cutoff PSF_CUTOFF
                            Cutoff for smoothing in units of arcseconds. (default: None)
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
                            The multiscale scale bias term provided to wsclean. (default: None)
    --absmem ABSMEM       Absolute memory limit in GB (default: None)
    --make_residual_cubes
                            Create residual cubes as well as cubes from restored images. (default: False)
    --ms_glob_pattern MS_GLOB_PATTERN
                            The pattern used to search for measurement sets. (default: scienceData*_averaged_cal.leakage.ms)
    --data_column DATA_COLUMN
                            Which column in the measurement set to image. (default: CORRECTED_DATA)
    --skip_fix_ms         Do not apply the ASKAP MS corrections from the package fixms. (default: False)

    pipeline flow options:
    --imager_only         Only run the imager component of the pipeline. (default: False)
    --skip_imager         Skip imaging stage [False]. (default: False)
    --skip_cutout         Skip cutout stage [False]. (default: False)
    --skip_linmos         Skip LINMOS stage [False]. (default: False)
    --skip_cleanup        Skip cleanup stage [False]. (default: False)
    --skip_frion          Skip cleanup stage [False]. (default: False)
    --skip_rmsynth        Skip RM Synthesis stage [False]. (default: False)
    --skip_rmclean        Skip RM-CLEAN stage [False]. (default: False)
    --skip_cat            Skip catalogue stage [False]. (default: False)

    output options:
    -v, --verbose         Verbose output [False]. (default: False)
    -vw, --verbose_worker
                            Verbose worker output [False]. (default: False)

    cutout arguments:
    -p PAD, --pad PAD     Number of beamwidths to pad around source [5]. (default: 5)
    --dryrun              Do a dry-run [False]. (default: False)

    RM-synth/CLEAN arguments:
    --dimension DIMENSION
                            How many dimensions for RMsynth [1d] or '3d'. (default: 1d)
    -m, --database        Add RMsynth data to MongoDB [False]. (default: False)
    --tt0 TT0             TT0 MFS image -- will be used for model of Stokes I -- also needs --tt1. (default: None)
    --tt1 TT1             TT1 MFS image -- will be used for model of Stokes I -- also needs --tt0. (default: None)
    --validate            Run on RMsynth Stokes I [False]. (default: False)
    --limit LIMIT         Limit number of sources [All]. (default: None)
    --own_fit             Use own Stokes I fit function [False]. (default: False)

    RM-tools arguments:
    -sp, --savePlots      save the plots [False]. (default: False)
    -w WEIGHTTYPE, --weightType WEIGHTTYPE
                            weighting [variance] (all 1s) or 'uniform'. (default: variance)
    --fit_function FIT_FUNCTION
                            Stokes I fitting function: 'linear' or ['log'] polynomials. (default: log)
    -t, --fitRMSF         Fit a Gaussian to the RMSF [False] (default: False)
    -l PHIMAX_RADM2, --phiMax_radm2 PHIMAX_RADM2
                            Absolute max Faraday depth sampled (overrides NSAMPLES) [Auto]. (default: None)
    -d DPHI_RADM2, --dPhi_radm2 DPHI_RADM2
                            Width of Faraday depth channel [Auto]. (default: None)
    -s NSAMPLES, --nSamples NSAMPLES
                            Number of samples across the FWHM RMSF. (default: 5)
    -o POLYORD, --polyOrd POLYORD
                            polynomial order to fit to I spectrum [3]. (default: 3)
    -i, --noStokesI       ignore the Stokes I spectrum [False]. (default: False)
    --showPlots           show the plots [False]. (default: False)
    -R, --not_RMSF        Skip calculation of RMSF? [False] (default: False)
    -rmv, --rm_verbose    Verbose RMsynth/CLEAN [False]. (default: False)
    -D, --debug           turn on debugging messages & plots [False]. (default: False)
    -c CUTOFF, --cutoff CUTOFF
                            CLEAN cutoff (+ve = absolute, -ve = sigma) [-3]. (default: -3)
    -n MAXITER, --maxIter MAXITER
                            maximum number of CLEAN iterations [10000]. (default: 10000)
    -g GAIN, --gain GAIN  CLEAN loop gain [0.1]. (default: 0.1)
    --window WINDOW       Further CLEAN in mask to this threshold [False]. (default: None)
    --ionex_server IONEX_SERVER
                            IONEX server [ftp://ftp.aiub.unibe.ch/CODE/]. (default: ftp://ftp.aiub.unibe.ch/CODE/)
    --ionex_prefix IONEX_PREFIX
                            IONEX prefix. (default: codg)
    --ionex_proxy_server IONEX_PROXY_SERVER
                            Proxy server [None]. (default: None)
    --ionex_formatter IONEX_FORMATTER
                            IONEX formatter [None]. (default: None)
    --ionex_predownload   Pre-download IONEX files [False]. (default: False)

    catalogue arguments:
    --outfile OUTFILE     File to save table to [None]. (default: None)

    Args that start with '--' (eg. --psf_cutoff) can also be set in a config file (.default_config.txt or specified via --config). Config file syntax allows:
    key=value, flag=true, stuff=[a,b,c] (for details, see syntax at https://goo.gl/R74nmi). If an arg is specified in more than one place, then commandline
    values override config file values which override defaults.


You can optionally pass a configuration file (with the :code:`--config` argument) to set the options you prefer. An example file in contained in :file:`arrakis/.default_config.txt`:

.. code-block::

    # Arrakis default config
    [General options]
    # host: # Host of mongodb.
    # username: # Username of mongodb.
    # password: # Password of mongodb.
    port: 8787 # Port to run Dask dashboard on.
    # port_forward: # Platform to fowards dask port [None].
    # dask_config: # Config file for Dask SlurmCLUSTER.
    # holofile:
    yanda: "1.3.0"

    [Flow options]
    skip_cutout: False
    skip_linmos: False
    skip_cleanup: False
    skip_rmsynth: False
    skip_rmclean: False
    skip_cat: False

    [Output options]
    verbose: True # Verbose output
    verbose_worker: False # Verbose worker output

    [Cutout options]
    pad: 3 # Number of beamwidths to pad around source
    dryrun: False # Do a dry-run

    [RM-synth/CLEAN options]
    dimension: 1d # How many dimensions for RMsynth [1d] or '3d'
    database: True # Add RMsynth data to MongoDB
    # tt0: # TT0 MFS image -- will be used for model of Stokes I -- also needs --tt1.
    # tt1: # TT1 MFS image -- will be used for model of Stokes I -- also needs --tt0.
    validate: False # Run on RMsynth Stokes I
    # limit: # Limit number of sources
    own_fit: False # Use own fit for RMsynth Stokes I

    [RM-tools options]
    savePlots: False # save the plots
    weightType: variance # weighting [uniform] (all 1s) or 'variance'
    fit_function: log # Stokes I fitting function
    fitRMSF: True # Fit a gaussian to the RMSF
    # phiMax_radm2: # Absolute max Faraday depth sampled (overrides NSAMPLES)
    # dPhi_radm2: # Width of Faraday depth channel
    nSamples: 5 # Number of samples across the FWHM RMSF
    polyOrd: -3 # polynomial order to fit to I spectrum
    noStokesI: False # ignore the Stokes I spectrum
    showPlots: False # show the plots
    not_RMSF: False # Skip calculation of RMSF
    rm_verbose: False # Verbose RMsynth/CLEAN
    debug: False # turn on debugging messages & plots
    cutoff: -3 # CLEAN cutoff (+ve = absolute, -ve = sigma)
    maxIter: 10000 # maximum number of CLEAN iterations
    gain: 0.1 # CLEAN loop gain
    window: False # use a windowed CLEAN

    [Catalogue options]
    # outfile: # File to save table to
    format: fits # Format for output file


For extra information you can refer to the API:

* :py:mod:`arrakis.process_spice`

Similarly, you can merge multiple fields togther using: ::

    (spice) $ spice_region -h
    usage: spice_region [-h] [--config CONFIG] [--merge_name MERGE_NAME] [--fields FIELDS [FIELDS ...]] [--datadirs DATADIRS [DATADIRS ...]]
                        [--output_dir OUTPUT_DIR] [--epoch EPOCH] [--host HOST] [--username USERNAME] [--password PASSWORD] [--use_mpi]
                        [--port_forward PORT_FORWARD [PORT_FORWARD ...]] [--dask_config DASK_CONFIG] [--yanda YANDA] [--skip_merge] [--skip_rmsynth]
                        [--skip_rmclean] [--skip_cat] [-v] [--debugger] [-vw] [--dimension DIMENSION] [-m] [--tt0 TT0] [--tt1 TT1] [--validate]
                        [--limit LIMIT] [--own_fit] [-sp] [-w WEIGHTTYPE] [--fit_function FIT_FUNCTION] [-t] [-l PHIMAX_RADM2] [-d DPHI_RADM2] [-s NSAMPLES]
                        [-o POLYORD] [-i] [--showPlots] [-R] [-rmv] [-D] [-c CUTOFF] [-n MAXITER] [-g GAIN] [--window WINDOW] [--outfile OUTFILE]


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


    Arrakis regional pipeline. Before running make sure to start a session of mongodb e.g. $
    mongod --dbpath=/path/to/database --bind_ip $(hostname -i)

    optional arguments:
    -h, --help            show this help message and exit
    --config CONFIG       Config file path (default: None)
    --merge_name MERGE_NAME
                            Name of the merged region (default: None)
    --fields FIELDS [FIELDS ...]
                            RACS fields to mosaic - e.g. 2132-50A. (default: None)
    --datadirs DATADIRS [DATADIRS ...]
                            Directories containing cutouts (in subdir outdir/cutouts).. (default: None)
    --output_dir OUTPUT_DIR
                            Path to save merged data (in output_dir/merge_name/cutouts) (default: None)
    --epoch EPOCH         Epoch to read field data from (default: 0)
    --host HOST           Host of mongodb (probably $hostname -i). (default: None)
    --username USERNAME   Username of mongodb. (default: None)
    --password PASSWORD   Password of mongodb. (default: None)
    --use_mpi             Use Dask-mpi to parallelise -- must use srun/mpirun to assign resources. (default: False)
    --port_forward PORT_FORWARD [PORT_FORWARD ...]
                            Platform to fowards dask port [None]. (default: None)
    --dask_config DASK_CONFIG
                            Config file for Dask SlurmCLUSTER. (default: None)
    --yanda YANDA         Yandasoft version to pull from DockerHub [1.3.0]. (default: 1.3.0)

    pipeline flow options:
    --skip_merge          Skip merge stage [False]. (default: False)
    --skip_rmsynth        Skip RM Synthesis stage [False]. (default: False)
    --skip_rmclean        Skip RM-CLEAN stage [False]. (default: False)
    --skip_cat            Skip catalogue stage [False]. (default: False)

    output options:
    -v, --verbose         Verbose output [False]. (default: False)
    --debugger            Debug output [False]. (default: False)
    -vw, --verbose_worker
                            Verbose worker output [False]. (default: False)

    RM-synth/CLEAN arguments:
    --dimension DIMENSION
                            How many dimensions for RMsynth [1d] or '3d'. (default: 1d)
    -m, --database        Add RMsynth data to MongoDB [False]. (default: False)
    --tt0 TT0             TT0 MFS image -- will be used for model of Stokes I -- also needs --tt1. (default: None)
    --tt1 TT1             TT1 MFS image -- will be used for model of Stokes I -- also needs --tt0. (default: None)
    --validate            Run on RMsynth Stokes I [False]. (default: False)
    --limit LIMIT         Limit number of sources [All]. (default: None)
    --own_fit             Use own Stokes I fit function [False]. (default: False)

    RM-tools arguments:
    -sp, --savePlots      save the plots [False]. (default: False)
    -w WEIGHTTYPE, --weightType WEIGHTTYPE
                            weighting [variance] (all 1s) or 'uniform'. (default: variance)
    --fit_function FIT_FUNCTION
                            Stokes I fitting function: 'linear' or ['log'] polynomials. (default: log)
    -t, --fitRMSF         Fit a Gaussian to the RMSF [False] (default: False)
    -l PHIMAX_RADM2, --phiMax_radm2 PHIMAX_RADM2
                            Absolute max Faraday depth sampled (overrides NSAMPLES) [Auto]. (default: None)
    -d DPHI_RADM2, --dPhi_radm2 DPHI_RADM2
                            Width of Faraday depth channel [Auto]. (default: None)
    -s NSAMPLES, --nSamples NSAMPLES
                            Number of samples across the FWHM RMSF. (default: 5)
    -o POLYORD, --polyOrd POLYORD
                            polynomial order to fit to I spectrum [3]. (default: 3)
    -i, --noStokesI       ignore the Stokes I spectrum [False]. (default: False)
    --showPlots           show the plots [False]. (default: False)
    -R, --not_RMSF        Skip calculation of RMSF? [False] (default: False)
    -rmv, --rm_verbose    Verbose RMsynth/CLEAN [False]. (default: False)
    -D, --debug           turn on debugging messages & plots [False]. (default: False)
    -c CUTOFF, --cutoff CUTOFF
                            CLEAN cutoff (+ve = absolute, -ve = sigma) [-3]. (default: -3)
    -n MAXITER, --maxIter MAXITER
                            maximum number of CLEAN iterations [10000]. (default: 10000)
    -g GAIN, --gain GAIN  CLEAN loop gain [0.1]. (default: 0.1)
    --window WINDOW       Further CLEAN in mask to this threshold [False]. (default: None)

    catalogue arguments:
    --outfile OUTFILE     File to save table to [None]. (default: None)

    Args that start with '--' (eg. --merge_name) can also be set in a config file (.default_field_config.txt or specified via --config). Config file syntax
    allows: key=value, flag=true, stuff=[a,b,c] (for details, see syntax at https://goo.gl/R74nmi). If an arg is specified in more than one place, then
    commandline values override config file values which override defaults.

* :py:mod:`arrakis.process_region`

Helper scripts (mostly for bespoke purposes) are available on the commandline. See the API reference for more details.
