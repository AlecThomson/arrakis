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

print = functools.partial(print, flush=True)
