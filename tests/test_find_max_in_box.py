"""Tests for the Stokes I peak finder."""

from __future__ import annotations

import numpy as np
import pytest
from arrakis.rmsynth_oncuts import find_max_in_box
from astropy.io import fits
from astropy.wcs import WCS


def _write_cube(path, data: np.ndarray) -> WCS:
    """Write a (freq, dec, ra) cube with a simple SIN projection."""
    wcs = WCS(naxis=3)
    wcs.wcs.ctype = ["RA---SIN", "DEC--SIN", "FREQ"]
    wcs.wcs.crval = [180.0, -30.0, 1.0e9]
    wcs.wcs.crpix = [6.0, 6.0, 1.0]
    wcs.wcs.cdelt = [-0.001, 0.001, 1.0e6]
    fits.PrimaryHDU(data=data, header=wcs.to_header()).writeto(path)
    return wcs.celestial


@pytest.fixture
def cube(tmp_path):
    data = np.ones((3, 11, 11))
    data[:, 6, 4] = 10.0  # peak at (row=6, col=4)
    data[:, 5, 5] = np.nan  # blanked pixel at the target position
    wcs = _write_cube(tmp_path / "cube.i.fits", data)
    return tmp_path, wcs


def test_find_max_ignores_nan(cube) -> None:
    outdir, wcs = cube
    target = wcs.pixel_to_world(5, 5)
    peak = find_max_in_box(
        coord=target, box_size=5, field_dict={"i_file": "cube.i.fits"}, outdir=outdir
    )
    assert peak.separation(wcs.pixel_to_world(4, 6)).arcsec < 1e-6


def test_find_max_all_nan_keeps_target(tmp_path) -> None:
    wcs = _write_cube(tmp_path / "cube.i.fits", np.full((3, 11, 11), np.nan))
    target = wcs.pixel_to_world(5, 5)
    peak = find_max_in_box(
        coord=target, box_size=3, field_dict={"i_file": "cube.i.fits"}, outdir=tmp_path
    )
    assert peak.separation(target).arcsec < 1e-6


def test_find_max_missing_key(tmp_path) -> None:
    with pytest.raises(ValueError, match="i_file"):
        find_max_in_box(coord=None, box_size=3, field_dict={}, outdir=tmp_path)
