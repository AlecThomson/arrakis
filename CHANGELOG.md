# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [2.1.0] - 2023-12-14

### Added
- Full `prefect`  flows and tasks

### Removed
- `dask` delayed etc.


## [2.0.0] - 2023-10-16

### Added
- `pre-commit` hooks for autoformatting checks
- `imager` script for sweet imaging
- `prefect` backbone

## [1.0.0] - 2023-06-15

### Added

- This changelog!
- `scripts/tar_cubelets.py` and CLI hook
- `makecat.py`: Added `flag_blended_components` to identify and flag blended components. Adds `is_blended_flag`, `N_blended`, `blend_ratio` to the catalogue.
- Proper logging module

### Fixed

- `columns_possum.py`: Add new Stokes I fit flags and UCDs (plus others) and descriptions
- `scripts/casda_prepare.py`: Refactor to make considated products and make CASDA happy
- `scripts/fix_dr1_cat.py`: Added extra columns that needed to be fixed in DR1 e.g. sbid, start_time
- Typing in various places

### Changed

- Renamed `scripts/helloworld.py` to `scripts/hellow_mpi_world.py`
- `makecat.py`: Added `compute_local_rm_flag` function
- `rmsynth_oncuts.py` Added new Stokes I fit flags
- `utils.py`: Refactored Stokes I fitting to use dicts to track values
- Use local installs of customised packages

### Removed

- `submit/casda_pre_prep.sh`
- `submit/casda_prep_test.sh`
- ASKAP RACS database as a submodule (changes how `init_databse.py` ingests data)


## [0.2.0] - 2019-12-01

The before times...
