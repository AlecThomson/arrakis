# SPICE-RACS
**S**pectra and **P**olarization **I**n **C**utouts of **E**xtragalactic **S**ources from **R**ACS

Scripts for processing polarized RACS data products.


## Installation
```
pip install spiceracs
```
Or, after cloning this repo, please run:
```
cd SPICERACS/
pip install -e .
```
This will install the python dependencies and the command-line scrips. Please note: this package requires `python >= 3.6`.

### Third-party tools
Additionally, this package also requires `RM-Tools`, `java >= 8`, and `mongodb` to be installed. Details on how to install these can be found at their respective websites:
* [RM-Tools](https://github.com/CIRADA-Tools/RM)
* [Java](https://www.java.com/en/download/)
* [MongoDB](https://www.mongodb.com/what-is-mongodb)
* [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/)

### Parallelisation
If you wish to use MPI for parallelisation, please also make sure some implementation (such as [OpenMPI](https://www.open-mpi.org/)) is also installed. If MPI is not available, then python multiprocessing can be used. Please note that multiprocessing is not capable of going across nodes in an HPC environment. 

## Getting started
It is recommended to use the command-line tools. Alternatively, scripts for processing RACS data are located in `spiceracs/`. 

Bash scripts which run everything together on PBS or SBATCH systems are located in `submit/`.

To keep track of the many files, and associated metadata, these scripts use MongoDB. Additionally, the scripts are intended to be run in order (pipeline-style). 

Currently, the order to run each script is:

0. `mongod --dbpath /path/to/database/ &` -- This initialises MongoDB in a directory of your choosing.
1. `spicecutout` or `spiceracs/cutout.py` -- Produce cubelets from a RACS field using a Selavy table.
2. The datacubes can be removed from disk, if required.
Optional:
    3. `spiceunresolved` or `spiceracs/unresolved.py` -- Find unresolved sources from a Selavy catalogue.
    4. `spicemoments` or `spiceracs/moments_oncuts.py` -- Make Faraday moment maps for Farnes+ (2018) method.
    5. `spicepolfind` or `spiceracs/polfind.py`-- Find polarized sources in a RACS field using the Farnes+ (2018) method.
6. `spicermsynth` or `spiceracs/rmsynth_oncuts.py` -- Run RM synthesis on cutouts.
7. `spicermclean` or `spiceracs/rmclean_oncuts.py` -- Run RM-CLEAN on cutouts.
8.  `spicemakecat` or `spiceracs/makecat.py` -- Make RM catalogue from results.

### Tips for running on Pawsey/Galaxy
#### Multiprocessing
If you want to test on a single node using multiprocessing for parallelisation, run the pipeline using the following:

`salloc -n 1` -- request a single node on galaxy

`srun -n 1 -c 20 spice... OPTIONS --ncores 20` -- This will use 20 cores (the maximum available on a single galaxy node). Without the addition of `srun` the processes will run on a single physical core, and will be much slower.

#### MPI
If you want to use MPI, there can be issues with the communication between the executable `mpirun` and your installation of `mpi4py`. It is reccommended to use the module/local installation of MPI instead. To ensure all the libraries you need, it is easiest to make a conda virtual environment:
```
conda create -n spice python=3.6
conda activate spice
pip install spiceracs
```
Now when you want to run the task in parallel with `N` cores:
```
module load python/3.6.3
module load numpy
module load astropy
module load scipy
module load mpi4py
module load pandas
conda activate spice

cd /path/to/your/dir

srun -n N spice... OPTIONS --mpi
```
Make sure you have requested enough nodes for the number of cores you want!

## Acknowledging
### Third-party software
This package utilises a number of third-party libraries. Please acknowledge these, as appropriate, if you use these tools for your research.

List of third party libraries:
* [Numpy](https://numpy.org/)
* [Matplotlib](https://matplotlib.org/)
* [Astropy](https://www.astropy.org/)
* [RM-Tools](https://github.com/CIRADA-Tools/RM)
* [Spectral-Cube](https://spectral-cube.readthedocs.io/)
* [tqdm](https://tqdm.github.io/) 
* [MongoDB](https://www.mongodb.com/) / [pymongo](https://api.mongodb.com/python/current/) 
* [Schwimmbad](https://schwimmbad.readthedocs.io/)
* [AegeanTools](https://github.com/PaulHancock/Aegean)
* [STILTS](http://www.star.bristol.ac.uk/~mbt/stilts/)
* [pandas](https://pandas.pydata.org/)
* [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/)