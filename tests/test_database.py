"""Tests for field database lookups."""

from __future__ import annotations

import pytest
from arrakis.utils.database import get_field_observation, get_sbid_observation


def test_get_sbid_observation(field_col) -> None:
    assert get_sbid_observation(field_col, 100)["SCAN_TINT"] == 900.0


def test_get_sbid_observation_missing(field_col) -> None:
    with pytest.raises(ValueError, match="SBID 999"):
        get_sbid_observation(field_col, 999)


def test_field_observation_repeat_sbid(field_col) -> None:
    """A repeat (SELECT=0) SBID gets its own entry, not the SELECT=1 one."""
    obs = get_field_observation(field_col, "RACS_1213-25", sbid=100)
    assert obs["SBID"] == 100
    assert obs["SCAN_TINT"] == 900.0


def test_field_observation_no_sbid_uses_select(field_col) -> None:
    obs = get_field_observation(field_col, "RACS_1213-25")
    assert obs["SBID"] == 200


def test_field_observation_sbid_wrong_field(field_col) -> None:
    with pytest.raises(ValueError, match="not RACS_1213-25"):
        get_field_observation(field_col, "RACS_1213-25", sbid=300)


def test_field_observation_multiple_select(field_col) -> None:
    field_col.update_one({"SBID": 100}, {"$set": {"SELECT": 1}})
    with pytest.raises(ValueError, match="More than one SELECT=1"):
        get_field_observation(field_col, "RACS_1213-25")


def test_field_observation_no_select_uses_highest_sbid(field_col) -> None:
    field_col.update_many({"FIELD_NAME": "RACS_1213-25"}, {"$set": {"SELECT": 0}})
    assert get_field_observation(field_col, "RACS_1213-25")["SBID"] == 200


def test_field_observation_missing_field(field_col) -> None:
    with pytest.raises(ValueError, match="RACS_0000-00 not found"):
        get_field_observation(field_col, "RACS_0000-00")
