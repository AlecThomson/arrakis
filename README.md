# SPICE-RACS
**S**pectra and **P**olarization **I**n **C**utouts of **E**xtragalactic **S**ources from **R**ACS

Scripts for processing polarized RACS data products.


## Installation
After cloning this repo, please run:
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

### Parallelisation
If you wish to use MPI for parallelisation, please also make sure some implementation (such as [OpenMPI](https://www.open-mpi.org/)) is also installed. If MPI is not available, then python multiprocessing can be used. Please note that multiprocessing is not capable of going across nodes in an HPC environment. 

## Getting started
It is recommended to use the command-line tools. Alternatively, scripts for processing RACS data are located in `spiceracs/`. 

Bash scripts which run everything together on PBS or SBATCH systems are located in `submit/`.

To keep track of the many files, and associated metadata, these scripts use MongoDB. Additionally, the scripts are intended to be run in order (pipeline-style). 

Currently, the order to run each script is:

0. `mongod --dbpath /path/to/database/ &` -- This initialises MongoDB in a directory of your choosing.
1. The following stages require the datacubes to be on disk.
    1. `spicecutout` or `spiceracs/cutout.py` -- Produce cubelets from a RACS field using a Selavy table.
    2. `spiceunresolved` or `spiceracs/unresolved.py` -- Find unresolved sources from a Selavy catalogue.
    3. `spicemoments` or `spiceracs/moments.py` -- Make Faraday moment maps for Farnes+ (2018) method.
2. The datacubes can be removed from disk, if required.
3. `spicepolfind` or `spiceracs/polfind.py`-- Find polarized sources in a RACS field using the Farnes+ (2018) method.
4. `spiceracs/rmsynth_oncuts.py` -- Run RM synthesis on unresolved, polarized sources.
5. ?????
6. Profit

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
