# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [2.2.0] - 2024-04-11
### What's Changed

* Allow SBID to passed as an argument
    * This will enable a 'single field mode'
    * Database queries / updates changes to support this
* Unified ArgParse mode
    * Much easier argument parsing
    * Now reused amongst modules
* Fixes to typing
    * Much better use of `pathlib.Path` and `pandas`
* SBID by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/62

**Full Changelog**: https://github.com/AlecThomson/arrakis/compare/v2.1.7...v2.1.8

## [2.1.7] - 2024-04-03
### What's Changed

Updated documentation.

**Full Changelog**: https://github.com/AlecThomson/arrakis/compare/v2.1.6...v2.1.7

## [2.1.6] - 2024-04-01
### What's Changed
* Image docs by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/60


**Full Changelog**: https://github.com/AlecThomson/arrakis/compare/v2.1.5...v2.1.6

## [2.1.5] - 2024-03-27
### What's Changed
* [pre-commit.ci] pre-commit autoupdate by @pre-commit-ci in https://github.com/AlecThomson/arrakis/pull/56
* Tempdir by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/58
* Clean by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/59


**Full Changelog**: https://github.com/AlecThomson/arrakis/compare/v2.1.4...v2.1.5

## [2.1.4] - 2024-03-20
### What's Changed
* DR2 preparation improvements by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/57

**Full Changelog**: https://github.com/AlecThomson/arrakis/compare/v2.1.3...v2.1.4

## [2.1.3] - 2024-02-05

### What's Changed
* [pre-commit.ci] pre-commit autoupdate by @pre-commit-ci in https://github.com/AlecThomson/arrakis/pull/54
* Noise fix by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/55

**Full Changelog**: https://github.com/AlecThomson/arrakis/compare/v2.1.2...v2.1.3

## [2.1.2] - 2024-01-17

Hotfix updates.

## Fixed
- Base Python version bumped to 3.10 for docs
- Enable `-no-mf-weighting` in WSClean
- Ensure FixMS skipped is used in pipeline
- Pass data in correct order to RM-Tools

**Full Changelog**: https://github.com/AlecThomson/arrakis/compare/v2.1.1...v2.1.2

## [2.1.1] - 2024-01-16

### What's Changed
* [pre-commit.ci] pre-commit autoupdate by @pre-commit-ci in https://github.com/AlecThomson/arrakis/pull/48
* Create dependabot.yml by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/50
* [pre-commit.ci] pre-commit autoupdate by @pre-commit-ci in https://github.com/AlecThomson/arrakis/pull/51
* [pre-commit.ci] pre-commit autoupdate by @pre-commit-ci in https://github.com/AlecThomson/arrakis/pull/52
* In-memory cutouts by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/53

### New Contributors
* @pre-commit-ci made their first contribution in https://github.com/AlecThomson/arrakis/pull/48

**Full Changelog**: https://github.com/AlecThomson/arrakis/compare/v2.1.0...v2.1.1

## [2.1.0] - 2023-12-14

### What's Changed
* Fix weights by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/46
* described the scalability of prefect/postgres by @tjgalvin in https://github.com/AlecThomson/arrakis/pull/44
* Migrate the whole shebang to Prefect by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/49


**Full Changelog**: https://github.com/AlecThomson/arrakis/compare/v2.0.0...v2.1.0

### Added
- Full `prefect` flows and tasks
- More `prefect` documentation

### Fixed
- Incorrect image weighting for LINMOS

### Removed
- `dask` delayed etc.


## [2.0.0] - 2023-10-16

### Added
- `pre-commit` hooks for autoformatting checks
- `imager` script for sweet imaging
- `prefect` backbone

## [1.0.0] - 2023-06-15
Corresponds with DR1 paper.

### What's Changed
* Dev by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/1
* Simplify readme by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/2
* Don't need these by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/3
* Clean up tests by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/4
* Fix scripts and RM-synth by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/5
* Merge dev by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/7
* Dev by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/9
* Remove RACS db by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/14
* Arrakis by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/16
* Fix RACS db pathing by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/15
* Docs update by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/18
* Dev by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/23
* Dev by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/24
* Dev by @AlecThomson in https://github.com/AlecThomson/arrakis/pull/25

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
