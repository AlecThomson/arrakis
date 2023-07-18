#!/usr/bin/env python
"""Coordinate utilities"""

import copy
import dataclasses
import functools
import json
import logging
import os
import shlex
import stat
import subprocess
import time
import warnings
from dataclasses import asdict, dataclass, make_dataclass
from functools import partial
from glob import glob
from itertools import zip_longest
from os import name
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import astropy.units as u
import dask.array as da
import dask.distributed as distributed
import numpy as np
import pymongo
from astropy.coordinates import SkyCoord
from astropy.coordinates.angles import dms_tuple, hms_tuple
from astropy.io import fits
from astropy.stats import akaike_info_criterion_lsq
from astropy.table import Table
from astropy.utils.exceptions import AstropyWarning
from astropy.wcs import WCS
from casacore.tables import table
from casatasks import listobs
from dask import delayed
from dask.delayed import Delayed
from dask.distributed import Client, get_client
from distributed.client import futures_of
from distributed.diagnostics.progressbar import ProgressBar
from distributed.utils import LoopRunner, is_kernel
from FRion.correct import find_freq_axis
from prefect_dask import get_dask_client
from pymongo.collection import Collection
from scipy.optimize import curve_fit
from scipy.stats import normaltest
from spectral_cube import SpectralCube
from spectral_cube.utils import SpectralCubeWarning
from tornado.ioloop import IOLoop
from tqdm.auto import tqdm, trange

from arrakis.logger import logger

warnings.filterwarnings(action="ignore", category=SpectralCubeWarning, append=True)
warnings.simplefilter("ignore", category=AstropyWarning)


def deg_to_hms(deg: float) -> hms_tuple:
    """Convert degree to hms without astropy.

    Args:
        deg (float): Decimal degrees

    Returns:
        hms_tuple: HMS, like coord.ra.hms
    """
    h_per_d = 24 / 360
    hours = deg * h_per_d
    hour = float(int(hours))
    minutes = (hours - hour) * 60
    minute = float(int(minutes))
    seconds = (minutes - minute) * 60
    return hms_tuple(hour, minute, seconds)


def deg_to_dms(deg: float) -> dms_tuple:
    """Convert degree to hms without astropy.

    Args:
        deg (float): Decimal degrees

    Returns:
        hms_tuple: DMS, like coord.dec.dms
    """
    degree = float(int(deg))
    minutes = (deg - degree) * 60
    minute = float(int(minutes))
    seconds = (minutes - minute) * 60
    return dms_tuple(degree, minute, seconds)


def coord_to_string(coord: SkyCoord) -> Tuple[str, str]:
    """Convert coordinate to string without astropy

    Args:
        coord (SkyCoord): Coordinate

    Returns:
        Tuple[str,str]: Tuple of RA string, Dec string
    """
    ra = coord.ra
    dec = coord.dec

    ra_hms = deg_to_hms(ra.value)
    dec_dms = deg_to_dms(dec.value)

    ra_str = f"{ra_hms.h:02.0f}:{ra_hms.m:02.0f}:{ra_hms.s:06.3f}"
    dec_str = f"{dec_dms.d:02.0f}:{abs(dec_dms.m):02.0f}:{abs(dec_dms.s):05.2f}"
    return ra_str, dec_str
