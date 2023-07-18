#!/usr/bin/env python
"""Coordinate utilities"""

import warnings
from typing import Tuple

from astropy.coordinates import SkyCoord
from astropy.coordinates.angles import dms_tuple, hms_tuple
from astropy.utils.exceptions import AstropyWarning
from spectral_cube.utils import SpectralCubeWarning

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
