Running the pipeline
--------------------
So you're ready to run the pipeline? The Arrakis pipeline picks up from whatever imaging routine you used to create the image cubes. You'll need to have image cubes in a working directory all convolved to a common spatial resolution. You should also make sure the names of the image cubes are consistent with the naming convention used in the ASKAPsoft pipline i.e. ::

    working_directory/image.restored.{i,q,u}*contcube*beam{00..36}.conv.fits

:code:`spice_process` and :code:`spice_field` orchestrate the pipeline flow using `Prefect <https://prefect.io>`_ and `Dask <https://dask.org>`_. These script calls the other :code:`arrakis` modules to do the work. You can control which modules run in the configuration of :code:`spice_process` or :code:`spice_field`. :code:`spice_process` operates on the level of a single RACS fields, whereas :code:`spice_field` merges multiple fields togther. You will need to run :code:`spice_process` on at least two fields before calling :code:`spice_field`. After running :code:`spice_process` or :code:`spice_field` you can run :code:`spice_cat` to produce a just a catalogue from the database values.

.. image:: flow.png
    :alt: Left floating image
    :class: with-shadow float-left

Details of each module can be found in the API documentation. But broadly the stages are:
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
    usage: spice_process [-h] [--config CONFIG] [--host HOST] [--username USERNAME] [--password PASSWORD] [--use_mpi] [--port_forward PORT_FORWARD [PORT_FORWARD ...]] [--dask_config DASK_CONFIG]
                        [--holofile HOLOFILE] [--yanda YANDA] [--skip_cutout] [--skip_linmos] [--skip_cleanup] [--skip_frion] [--skip_rmsynth] [--skip_rmclean] [--skip_cat] [-v] [-vw] [-p PAD]
                        [--dryrun] [--dimension DIMENSION] [-m] [--tt0 TT0] [--tt1 TT1] [--validate] [--limit LIMIT] [--own_fit] [-sp] [-w WEIGHTTYPE] [--fit_function FIT_FUNCTION] [-t]
                        [-l PHIMAX_RADM2] [-d DPHI_RADM2] [-s NSAMPLES] [-o POLYORD] [-i] [--showPlots] [-R] [-rmv] [-D] [-c CUTOFF] [-n MAXITER] [-g GAIN] [--window WINDOW] [--outfile OUTFILE]
                        field datadir


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



    positional arguments:
    field                 Name of field (e.g. 2132-50A).
    datadir               Directory containing data cubes in FITS format.

    optional arguments:
    -h, --help            show this help message and exit
    --config CONFIG       Config file path
    --host HOST           Host of mongodb (probably $hostname -i).
    --username USERNAME   Username of mongodb.
    --password PASSWORD   Password of mongodb.
    --use_mpi             Use Dask-mpi to parallelise -- must use srun/mpirun to assign resources.
    --port_forward PORT_FORWARD [PORT_FORWARD ...]
                            Platform to fowards dask port [None].
    --dask_config DASK_CONFIG
                            Config file for Dask SlurmCLUSTER.
    --holofile HOLOFILE   Path to holography image
    --yanda YANDA         Yandasoft version to pull from DockerHub [1.3.0].

    pipeline flow options:
    --skip_cutout         Skip cutout stage [False].
    --skip_linmos         Skip LINMOS stage [False].
    --skip_cleanup        Skip cleanup stage [False].
    --skip_frion          Skip cleanup stage [False].
    --skip_rmsynth        Skip RM Synthesis stage [False].
    --skip_rmclean        Skip RM-CLEAN stage [False].
    --skip_cat            Skip catalogue stage [False].

    output options:
    -v, --verbose         Verbose output [False].
    -vw, --verbose_worker
                            Verbose worker output [False].

    cutout arguments:
    -p PAD, --pad PAD     Number of beamwidths to pad around source [5].
    --dryrun              Do a dry-run [False].

    RM-synth/CLEAN arguments:
    --dimension DIMENSION
                            How many dimensions for RMsynth [1d] or '3d'.
    -m, --database        Add RMsynth data to MongoDB [False].
    --tt0 TT0             TT0 MFS image -- will be used for model of Stokes I -- also needs --tt1.
    --tt1 TT1             TT1 MFS image -- will be used for model of Stokes I -- also needs --tt0.
    --validate            Run on RMsynth Stokes I [False].
    --limit LIMIT         Limit number of sources [All].
    --own_fit             Use own Stokes I fit function [False].

    RM-tools arguments:
    -sp, --savePlots      save the plots [False].
    -w WEIGHTTYPE, --weightType WEIGHTTYPE
                            weighting [variance] (all 1s) or 'uniform'.
    --fit_function FIT_FUNCTION
                            Stokes I fitting function: 'linear' or ['log'] polynomials.
    -t, --fitRMSF         Fit a Gaussian to the RMSF [False]
    -l PHIMAX_RADM2, --phiMax_radm2 PHIMAX_RADM2
                            Absolute max Faraday depth sampled (overrides NSAMPLES) [Auto].
    -d DPHI_RADM2, --dPhi_radm2 DPHI_RADM2
                            Width of Faraday depth channel [Auto].
    -s NSAMPLES, --nSamples NSAMPLES
                            Number of samples across the FWHM RMSF.
    -o POLYORD, --polyOrd POLYORD
                            polynomial order to fit to I spectrum [3].
    -i, --noStokesI       ignore the Stokes I spectrum [False].
    --showPlots           show the plots [False].
    -R, --not_RMSF        Skip calculation of RMSF? [False]
    -rmv, --rm_verbose    Verbose RMsynth/CLEAN [False].
    -D, --debug           turn on debugging messages & plots [False].
    -c CUTOFF, --cutoff CUTOFF
                            CLEAN cutoff (+ve = absolute, -ve = sigma) [-3].
    -n MAXITER, --maxIter MAXITER
                            maximum number of CLEAN iterations [10000].
    -g GAIN, --gain GAIN  CLEAN loop gain [0.1].
    --window WINDOW       Further CLEAN in mask to this threshold [False].

    catalogue arguments:
    --outfile OUTFILE     File to save table to [None].

    Args that start with '--' (eg. --host) can also be set in a config file (.default_config.txt or specified via --config). Config file syntax allows: key=value, flag=true, stuff=[a,b,c] (for
    details, see syntax at https://goo.gl/R74nmi). If an arg is specified in more than one place, then commandline values override config file values which override defaults.


You can optionally pass a configuration file (with the :code:`--config` argument) to set the options you prefer. An example file in contained in :file:`arrakis/.default_config.txt`.

For extra information you can refer to the API:

* :py:mod:`arrakis.process_spice`

Similarly, you can merge multiple fields togther using: ::

    (spice) $ spice_region -h
    usage: spice_region [-h] [--config CONFIG] [--merge_name MERGE_NAME] [--fields FIELDS [FIELDS ...]] [--datadirs DATADIRS [DATADIRS ...]] [--output_dir OUTPUT_DIR] [--host HOST]
                        [--username USERNAME] [--password PASSWORD] [--use_mpi] [--port_forward PORT_FORWARD [PORT_FORWARD ...]] [--dask_config DASK_CONFIG] [--yanda YANDA] [--skip_merge]
                        [--skip_rmsynth] [--skip_rmclean] [--skip_cat] [-v] [--debugger] [-vw] [--dimension DIMENSION] [-m] [--tt0 TT0] [--tt1 TT1] [--validate] [--limit LIMIT] [--own_fit] [-sp]
                        [-w WEIGHTTYPE] [--fit_function FIT_FUNCTION] [-t] [-l PHIMAX_RADM2] [-d DPHI_RADM2] [-s NSAMPLES] [-o POLYORD] [-i] [--showPlots] [-R] [-rmv] [-D] [-c CUTOFF] [-n MAXITER]
                        [-g GAIN] [--window WINDOW] [--outfile OUTFILE]


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



    optional arguments:
    -h, --help            show this help message and exit
    --config CONFIG       Config file path
    --merge_name MERGE_NAME
                            Name of the merged region
    --fields FIELDS [FIELDS ...]
                            RACS fields to mosaic - e.g. 2132-50A.
    --datadirs DATADIRS [DATADIRS ...]
                            Directories containing cutouts (in subdir outdir/cutouts)..
    --output_dir OUTPUT_DIR
                            Path to save merged data (in output_dir/merge_name/cutouts)
    --host HOST           Host of mongodb (probably $hostname -i).
    --username USERNAME   Username of mongodb.
    --password PASSWORD   Password of mongodb.
    --use_mpi             Use Dask-mpi to parallelise -- must use srun/mpirun to assign resources.
    --port_forward PORT_FORWARD [PORT_FORWARD ...]
                            Platform to fowards dask port [None].
    --dask_config DASK_CONFIG
                            Config file for Dask SlurmCLUSTER.
    --yanda YANDA         Yandasoft version to pull from DockerHub [1.3.0].

    pipeline flow options:
    --skip_merge          Skip merge stage [False].
    --skip_rmsynth        Skip RM Synthesis stage [False].
    --skip_rmclean        Skip RM-CLEAN stage [False].
    --skip_cat            Skip catalogue stage [False].

    output options:
    -v, --verbose         Verbose output [False].
    --debugger            Debug output [False].
    -vw, --verbose_worker
                            Verbose worker output [False].

    RM-synth/CLEAN arguments:
    --dimension DIMENSION
                            How many dimensions for RMsynth [1d] or '3d'.
    -m, --database        Add RMsynth data to MongoDB [False].
    --tt0 TT0             TT0 MFS image -- will be used for model of Stokes I -- also needs --tt1.
    --tt1 TT1             TT1 MFS image -- will be used for model of Stokes I -- also needs --tt0.
    --validate            Run on RMsynth Stokes I [False].
    --limit LIMIT         Limit number of sources [All].
    --own_fit             Use own Stokes I fit function [False].

    RM-tools arguments:
    -sp, --savePlots      save the plots [False].
    -w WEIGHTTYPE, --weightType WEIGHTTYPE
                            weighting [variance] (all 1s) or 'uniform'.
    --fit_function FIT_FUNCTION
                            Stokes I fitting function: 'linear' or ['log'] polynomials.
    -t, --fitRMSF         Fit a Gaussian to the RMSF [False]
    -l PHIMAX_RADM2, --phiMax_radm2 PHIMAX_RADM2
                            Absolute max Faraday depth sampled (overrides NSAMPLES) [Auto].
    -d DPHI_RADM2, --dPhi_radm2 DPHI_RADM2
                            Width of Faraday depth channel [Auto].
    -s NSAMPLES, --nSamples NSAMPLES
                            Number of samples across the FWHM RMSF.
    -o POLYORD, --polyOrd POLYORD
                            polynomial order to fit to I spectrum [3].
    -i, --noStokesI       ignore the Stokes I spectrum [False].
    --showPlots           show the plots [False].
    -R, --not_RMSF        Skip calculation of RMSF? [False]
    -rmv, --rm_verbose    Verbose RMsynth/CLEAN [False].
    -D, --debug           turn on debugging messages & plots [False].
    -c CUTOFF, --cutoff CUTOFF
                            CLEAN cutoff (+ve = absolute, -ve = sigma) [-3].
    -n MAXITER, --maxIter MAXITER
                            maximum number of CLEAN iterations [10000].
    -g GAIN, --gain GAIN  CLEAN loop gain [0.1].
    --window WINDOW       Further CLEAN in mask to this threshold [False].

    catalogue arguments:
    --outfile OUTFILE     File to save table to [None].

    Args that start with '--' (eg. --merge_name) can also be set in a config file (.default_field_config.txt or specified via --config). Config file syntax allows: key=value, flag=true,
    stuff=[a,b,c] (for details, see syntax at https://goo.gl/R74nmi). If an arg is specified in more than one place, then commandline values override config file values which override defaults.

* :py:mod:`arrakis.process_region`

Helper scripts (mostly for bespoke purposes) are available on the commandline. See the API reference for more details.