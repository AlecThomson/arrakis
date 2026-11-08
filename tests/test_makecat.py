"""Tests for observation times in makecat."""

from __future__ import annotations

import astropy.units as u
import numpy as np
import pytest
from arrakis.makecat import (
    get_integration_time,
    get_sbid_observation,
    get_sbid_start_time,
    set_sbid_columns,
)
from astropy.units import cds
from rmtable import RMTable


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


def test_set_sbid_columns_and_epoch(field_col) -> None:
    """Processing SBID 100 overwrites the catalogued observation's values."""
    # Rows catalogued in the latest SBID of this field and in another field
    cat = _make_cat(["RACS_1213-25", "RACS_1237-25"], [59500.5, 59100.75])
    cat.add_column(["200", "300"], name="sbid")
    cat = set_sbid_columns(cat, field_col, field="RACS_1213-25", sbid=100)

    assert list(cat["sbid"]) == ["100", "100"]
    assert list(cat["tile_id"]) == ["RACS_1213-25"] * 2
    assert list(cat["field_name"]) == ["RACS_1213-25"] * 2
    assert cat["start_time"].unit == cds.MJD
    assert np.allclose(cat["start_time"], 59000.25)

    # Same steps as makecat.main
    tints = get_integration_time.fn(cat, field_col, sbid=100)
    epoch = cat["start_time"] + tints / 2
    expected = 59000.25 + (450.0 * u.s).to_value(u.day)
    assert np.allclose(epoch.to_value(cds.MJD), expected)


def test_set_sbid_columns_no_truncation(field_col) -> None:
    """Longer SBIDs and field names are not cut to the old column width."""
    field_col.insert_one(
        {
            "FIELD_NAME": "RACS_1213-25A",
            "SBID": 50230,
            "SELECT": 1,
            "SCAN_START": (60000.0 * u.day).to_value(u.s),
            "SCAN_TINT": 900.0,
        }
    )
    cat = _make_cat(["RACS_1213-25"], [59500.5])
    cat.add_column(["8584"], name="sbid")
    cat = set_sbid_columns(cat, field_col, field="RACS_1213-25A", sbid=50230)

    assert list(cat["sbid"]) == ["50230"]
    assert list(cat["tile_id"]) == ["RACS_1213-25A"]
