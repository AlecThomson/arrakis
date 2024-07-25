"""Typing utilities."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.units import Quantity
from rmtable import RMTable

ArrayLike = np.ndarray | pd.Series | pd.DataFrame | SkyCoord | Quantity
TableLike = RMTable | Table
PathLike = str | Path
