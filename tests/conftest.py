"""Shared test fixtures."""

from __future__ import annotations

import astropy.units as u
import mongomock
import pytest


@pytest.fixture
def field_col() -> mongomock.collection.Collection:
    """In-memory field database: a re-observed field (two SBIDs, only the
    later one SELECT=1) and a neighbouring field. SCAN_START is MJD in
    seconds, as in the real field database."""
    docs = [
        ("RACS_1213-25", 100, 0, 59000.25, 900.0),
        ("RACS_1213-25", 200, 1, 59500.5, 1039.0),
        ("RACS_1237-25", 300, 1, 59100.75, 840.0),
    ]
    col = mongomock.MongoClient().db.fields
    col.insert_many(
        [
            {
                "FIELD_NAME": field_name,
                "SBID": sbid,
                "SELECT": select,
                "SCAN_START": (start_mjd * u.day).to_value(u.s),
                "SCAN_TINT": tint,
            }
            for field_name, sbid, select, start_mjd, tint in docs
        ]
    )
    return col
