#!/usr/bin/env python3
"""Typing utilities"""

from typing import TypeVar

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.units import Quantity
from rmtable import RMTable

ArrayLike = TypeVar(
    "ArrayLike", np.ndarray, pd.Series, pd.DataFrame, SkyCoord, Quantity
)
TableLike = TypeVar("TableLike", RMTable, Table)
