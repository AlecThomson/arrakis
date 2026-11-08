"""Tests for observation times in makecat."""

from __future__ import annotations

import astropy.units as u
import mongomock
import numpy as np
import pytest
from arrakis.makecat import (
    get_integration_time,
    get_sbid_observation,
    get_sbid_start_time,
)
from astropy.units import cds
from rmtable import RMTable

SECONDS_PER_DAY = 86400.0

# A re-observed field (two SBIDs, only the later one SELECT=1) and a
# neighbouring field. SCAN_START is MJD in seconds, as in the field database.
FIELD_DOCS = [
    {
        "FIELD_NAME": "RACS_1213-25",
        "SBID": 100,
        "SELECT": 0,
        "SCAN_START": 59000.25 * SECONDS_PER_DAY,
        "SCAN_TINT": 900.0,
    },
    {
        "FIELD_NAME": "RACS_1213-25",
        "SBID": 200,
        "SELECT": 1,
        "SCAN_START": 59500.5 * SECONDS_PER_DAY,
        "SCAN_TINT": 1039.0,
    },
    {
        "FIELD_NAME": "RACS_1237-25",
        "SBID": 300,
        "SELECT": 1,
        "SCAN_START": 59100.75 * SECONDS_PER_DAY,
        "SCAN_TINT": 840.0,
    },
]


@pytest.fixture
def field_col() -> mongomock.collection.Collection:
    """In-memory stand-in for the field database collection."""
    col = mongomock.MongoClient().db.fields
    col.insert_many([dict(d) for d in FIELD_DOCS])
    return col


def _make_cat(tile_ids: list[str], start_times: list[float]) -> RMTable:
    cat = RMTable()
    cat.add_column(tile_ids, name="tile_id")
    cat.add_column(np.array(start_times) * cds.MJD, name="start_time")
    return cat


def test_get_sbid_observation_missing(field_col) -> None:
    with pytest.raises(ValueError, match="SBID 999"):
        get_sbid_observation(field_col, 999)


@pytest.mark.parametrize(
    ("sbid", "expected_mjd"), [(100, 59000.25), (200, 59500.5), (300, 59100.75)]
)
def test_get_sbid_start_time(field_col, sbid: int, expected_mjd: float) -> None:
    assert get_sbid_start_time(field_col, sbid).mjd == pytest.approx(expected_mjd)


def test_integration_time_repeat_sbid(field_col) -> None:
    """A repeat (SELECT=0) SBID gets its own SCAN_TINT, not the latest SBID's."""
    cat = _make_cat(["RACS_1213-25"] * 3, [0.0] * 3)
    tints = get_integration_time.fn(cat, field_col, sbid=100)
    assert u.allclose(tints, [900.0] * 3 * u.s)


def test_integration_time_sbid_ignores_tile(field_col) -> None:
    """With an SBID, rows catalogued in other tiles get that SBID's SCAN_TINT."""
    cat = _make_cat(["RACS_1213-25", "RACS_1237-25"], [0.0, 0.0])
    tints = get_integration_time.fn(cat, field_col, sbid=300)
    assert u.allclose(tints, [840.0, 840.0] * u.s)


def test_integration_time_no_sbid_uses_select(field_col) -> None:
    """Without an SBID, each row's tile uses its SELECT=1 observation."""
    cat = _make_cat(["RACS_1213-25", "RACS_1237-25A"], [0.0, 0.0])
    tints = get_integration_time.fn(cat, field_col)
    assert u.allclose(tints, [1039.0, 840.0] * u.s)


def test_start_time_overwrite_and_epoch(field_col) -> None:
    """Overwriting start_time as in makecat.main gives a consistent epoch."""
    # Source finder start times from the catalogued (other) observation
    cat = _make_cat(["RACS_1213-25", "RACS_1237-25"], [59500.5, 59100.75])
    sbid = 100
    cat["start_time"] = get_sbid_start_time(field_col, sbid).mjd
    assert cat["start_time"].unit == cds.MJD
    assert np.allclose(cat["start_time"], 59000.25)

    tints = get_integration_time.fn(cat, field_col, sbid=sbid)
    epoch = cat["start_time"] + tints / 2
    expected = 59000.25 + 450.0 / SECONDS_PER_DAY
    assert np.allclose(epoch.to_value(cds.MJD), expected)
