"""Shared test fixtures."""

from __future__ import annotations

import mongomock
import pytest

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
