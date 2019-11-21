# SPICE-RACS
**S**pectra and **P**olarization **I**n **C**utouts of **E**xtragalactic **S**ources from **R**ACS

Scripts for processing polarized RACS data products.



## Getting started
Scripts for processing RACS data are located in `spiceracs/`. Bash scripts which run everything together on PBS or SBATCH systems are located in `submit/`.

To keep track of the many files, and associated metadata, these scripts use MongoDB. Additionally, the scripts are intended to be run in order (pipeline-style). 

Currently, the order to run each script is:

0. `mongod --dbpath /path/to/database/ &` -- This initialises MongoDB in a directory of your choosing.
1. `spiceracs/cutout.py` -- Produce cubelets from a RACS field using a Selavy table.
2. `spiceracs/unresolved.py` -- Find unresolved sources from a Selavy catalogue.
3. `spiceracs/polfind.py` -- Find polarized sources in a RACS field using the Farnes+ (2018) method.
4. `spiceracs/rmsynth_oncuts.py` -- Run RM synthesis on unresolved, polarized sources.
5. ?????
6. Profit
