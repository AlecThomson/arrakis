# SPICE-RACS
**S**pectra and **P**olarization **I**n **C**utouts of **E**xtragalactic sources from **R**ACS

Scripts for processing polarized RACS data products.


## Installation
First, ensure you have anaconda / miniconda installed.

After cloning this repo, please run:
```
cd spiceracs/
conda env create
```
This will install the python dependencies and the command-line scrips into a conda environment called `spice`, which can be activated by:
```
conda activate spice
```

### Parallelisation
Dask is used for parallelisation and job submission. The pipeline is currently configured for the `galaxy` supercomputer at Pawsey, which uses a Slurm job manager. This configuration lives in YAML file in `spiceracs/dask/jobqueue.yaml`. Either use this as is for galaxy, or add your own configuration by editing the file (see the dask [docs](https://jobqueue.dask.org/en/latest/configuration.html)). Then export this environment varible to let dask find the file:
```
export DASK_ROOT_CONFIG=/path/to/spiceracs/dask
```
You may wish to simply put this line in your `~/.bashrc`. If you use a machine with scheduler other than Slurm, you'll also need to edit the `cluster` variable in `spiceracs/processSPICE.py`. See [dask-jobqueue](https://jobqueue.dask.org/) for the available alternatives and documentation.

## Getting started
First you must initialise the SPICE-RACS database using MongoDB and the RACS catalogues. For example, you can start mongo using (for NUMA systems):
```
host=$(hostname -i)
database=/path/to/your/database/dir
mkdir $database
numactl --interleave=all mongod --dbpath=$database --bind_ip $host >> /dev/null &
```
Then run the initialisation script:
```
(spice) $ initSPICE -h
usage: initSPICE [-h] [-i ISLANDCAT] [-c COMPCAT] [-v] [-l] host

    
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

    
    SPICE-RACS Initialisation:
    
    Create MongoDB database from RACS catalogues.

    Before running make sure to start a session of mongodb e.g.
        $ mongod --dbpath=/path/to/database --bind_ip $(hostname -i)

    

positional arguments:
  host          Host of mongodb (probably $hostname -i).

optional arguments:
  -h, --help    show this help message and exit
  -i ISLANDCAT  Master island RACS catalogue.
  -c COMPCAT    Master component RACS catalogue.
  -v            Verbose output [False].
  -l            Load catalogue into database [False].
```

## Running the pipeline
With an initalised database you can call the pipeline:
```
(spice) $ processSPICE -h
usage: processSPICE [-h] [--config CONFIG] [-v] [-vw] [-p PAD] [--dryrun]
                    [--dimension DIMENSION] [-m] [--validate] [--limit LIMIT]
                    [-sp] [-w WEIGHTTYPE] [-t] [-l PHIMAX_RADM2]
                    [-d DPHI_RADM2] [-s NSAMPLES] [-o POLYORD] [-i]
                    [--showPlots] [-R] [-rmv] [-D] [-c CUTOFF] [-n MAXITER]
                    [-g GAIN] [--outfile OUTFILE] [-f FORMAT]
                    field datadir host

    
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

    
    SPICE-RACS pipeline.

    Before running make sure to start a session of mongodb e.g.
        $ mongod --dbpath=/path/to/database --bind_ip $(hostname -i)

    

positional arguments:
  field                 Name of field (e.g. 2132-50A).
  datadir               Directory containing data cubes in FITS format.
  host                  Host of mongodb (probably $hostname -i).

optional arguments:
  -h, --help            show this help message and exit
  --config CONFIG       Config file path

output options:
  -v, --verbose         Verbose output [True].
  -vw, --verbose_worker
                        Verbose worker output [False].

cutout arguments:
  -p PAD, --pad PAD     Number of beamwidths to pad around source [5].
  --dryrun              Do a dry-run [False].

RM-synth/CLEAN arguments:
  --dimension DIMENSION
                        How many dimensions for RMsynth [1d] or '3d'.
  -m, --database        Add RMsynth data to MongoDB [False].
  --validate            Run on RMsynth Stokes I [False].
  --limit LIMIT         Limit number of sources [All].

RM-tools arguments:
  -sp, --savePlots      save the plots [False].
  -w WEIGHTTYPE, --weightType WEIGHTTYPE
                        weighting [variance] (all 1s) or 'variance'.
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

catalogue arguments:
  --outfile OUTFILE     File to save table to [None].
  -f FORMAT, --format FORMAT
                        Format for output file [None].

Args that start with '--' (eg. -v) can also be set in a config file
(.default_config.txt or specified via --config). Config file syntax allows:
key=value, flag=true, stuff=[a,b,c] (for details, see syntax at
https://goo.gl/R74nmi). If an arg is specified in more than one place, then
commandline values override config file values which override defaults.
```

You can optionally pass a configuration file (with the `--config` argument) to set the options you prefer. An example file in contained in `spiceracs/.defailt_config.txt`.

## Acknowledging
### Third-party software
This package utilises a number of third-party libraries. Please acknowledge these, as appropriate, if you use these tools for your research.

List of third party libraries:
* [Numpy](https://numpy.org/)
* [SciPy](https://www.scipy.org/)
* [Matplotlib](https://matplotlib.org/)
* [Astropy](https://www.astropy.org/)
* [MongoDB](https://www.mongodb.com/) / [pymongo](https://api.mongodb.com/python/current/) 
* [Dask](https://dask.org/)
* [Prefect](https://www.prefect.io/)
* [RM-Tools](https://github.com/CIRADA-Tools/RM)
* [RMTable](https://github.com/Cameron-Van-Eck/RMTable)
* [Spectral-Cube](https://spectral-cube.readthedocs.io/)
* [tqdm](https://tqdm.github.io/) 
* [ConfigArgParse](https://github.com/bw2/ConfigArgParse) 
