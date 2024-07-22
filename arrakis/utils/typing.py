"""Typing utilities"""

from __future__ import annotations

from pathlib import Path
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
PathLike = TypeVar("PathLike", str, Path)
T = TypeVar("T")
