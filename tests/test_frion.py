"""Tests for the ionospheric Faraday rotation stage."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
import pytest
from arrakis import frion
from arrakis.utils.fitsutils import fix_header
from astropy.io import fits
from astropy.time import Time, TimeDelta

FIELD = "RACS_TEST"
ISLAND = "SRC1"
N_CHAN = 8


class Cutouts:
    """Synthetic cutouts with offline IONEX data."""

    def __init__(self, cutdir: Path, beam: dict[str, Any], island: dict[str, Any]):
        self.cutdir = cutdir
        self.beam = beam
        self.island = island
        self.start_time = Time("2020-01-08T02:00:00", scale="utc")
        self.end_time = self.start_time + TimeDelta(900 * u.s)
        self.freq = np.linspace(800e6, 1100e6, N_CHAN)


def make_header() -> fits.Header:
    """Cutout header with frequency on the last FITS axis."""
    header = fits.Header()
    for key, val in {
        "NAXIS": 4,
        "NAXIS1": 16,
        "NAXIS2": 16,
        "NAXIS3": 1,
        "NAXIS4": N_CHAN,
        "CTYPE1": "RA---SIN",
        "CRVAL1": 255.0,
        "CDELT1": -0.001,
        "CRPIX1": 8.0,
        "CUNIT1": "deg",
        "CTYPE2": "DEC--SIN",
        "CRVAL2": -30.0,
        "CDELT2": 0.001,
        "CRPIX2": 8.0,
        "CUNIT2": "deg",
        "CTYPE3": "STOKES",
        "CRVAL3": 1.0,
        "CDELT3": 1.0,
        "CRPIX3": 1.0,
        "CTYPE4": "FREQ",
        "CRVAL4": 800e6,
        "CDELT4": 300e6 / (N_CHAN - 1),
        "CRPIX4": 1.0,
        "CUNIT4": "Hz",
        "SPECSYS": "TOPOCENT",
    }.items():
        header[key] = val
    return header


@pytest.fixture
def cutouts(tmp_path: Path) -> Cutouts:
    ionex_dir = tmp_path / "IONEXdata"
    ionex_dir.mkdir()
    for ionex_file in (Path(__file__).parent / "data" / "ionex").glob("*.Z"):
        shutil.copy(ionex_file, ionex_dir)
        # FRion looks for the uncompressed, upper case file
        uncompressed = ionex_dir / ionex_file.with_suffix("").name.upper()
        with uncompressed.open("wb") as f:
            subprocess.run(["gzip", "-dc", str(ionex_file)], stdout=f, check=True)

    cutdir = tmp_path / "cutouts"
    island_dir = cutdir / ISLAND
    island_dir.mkdir(parents=True)
    rng = np.random.default_rng(42)
    shape = (N_CHAN, 1, 16, 16)
    for stokes in "iqu":
        data = rng.normal(size=shape).astype(np.float32)
        fits.writeto(island_dir / f"{ISLAND}.{stokes}.fits", data, make_header())

    beam = {
        "Source_ID": ISLAND,
        "beams": {
            FIELD: {
                f"{stokes}_file": f"{ISLAND}/{ISLAND}.{stokes}.fits" for stokes in "iqu"
            }
        },
    }
    island = {"Source_ID": ISLAND, "RA": 255.0, "Dec": -30.0}
    return Cutouts(cutdir=cutdir, beam=beam, island=island)


def predict(cutouts: Cutouts) -> frion.Prediction:
    return frion.predict_worker.fn(
        island=cutouts.island,
        field=FIELD,
        beam=cutouts.beam,
        start_time=cutouts.start_time,
        end_time=cutouts.end_time,
        freq=cutouts.freq,
        cutdir=cutouts.cutdir,
        server="ftp://ftp.aiub.unibe.ch/CODE/",
        prefix="codg",
        formatter="ftp.aiub.unibe.ch",
    )


def test_parser_defaults() -> None:
    args = frion.frion_parser().parse_args([])
    assert args.ionex_server == "ftp://ftp.aiub.unibe.ch/CODE/"
    assert args.ionex_prefix == "codg"
    assert args.ionex_formatter == "ftp.aiub.unibe.ch"
    assert args.ionex_proxy_server is None
    assert args.ionex_predownload is False


def test_parser_flags() -> None:
    args = frion.frion_parser().parse_args(
        [
            "--ionex_server",
            "https://example.org/ionex",
            "--ionex_prefix",
            "jplg",
            "--ionex_formatter",
            "cddis.gsfc.nasa.gov",
            "--ionex_proxy_server",
            "proxy.example.org",
            "--ionex_predownload",
        ]
    )
    assert args.ionex_server == "https://example.org/ionex"
    assert args.ionex_prefix == "jplg"
    assert args.ionex_formatter == "cddis.gsfc.nasa.gov"
    assert args.ionex_proxy_server == "proxy.example.org"
    assert args.ionex_predownload is True


def test_index_beams() -> None:
    beams = [{"Source_ID": "A"}, {"Source_ID": "B"}]
    assert frion.index_beams.fn({"Source_ID": "B"}, beams) == {"Source_ID": "B"}


def test_predict_worker(cutouts: Cutouts) -> None:
    prediction = predict(cutouts)

    predict_file = Path(prediction.predict_file)
    assert predict_file.exists()
    assert predict_file.parent == cutouts.cutdir / ISLAND
    assert predict_file.stem == f"{ISLAND}_ion"

    assert prediction.update._filter == {"Source_ID": ISLAND}
    ion = prediction.update._doc["$set"]["frion"]
    assert set(ion) == {"times", "RMs", "theta_real", "theta_imag"}

    # MJD seconds, within one 300 s step of the scan
    times = np.array(ion["times"])
    assert times.min() >= cutouts.start_time.mjd * 86400 - 300
    assert times.max() <= cutouts.end_time.mjd * 86400 + 300
    assert len(ion["RMs"]) == len(times)

    theta = np.array(ion["theta_real"]) + 1j * np.array(ion["theta_imag"])
    assert theta.shape == (N_CHAN,)
    assert np.all(np.abs(theta) > 0.9)
    assert np.all(np.abs(theta) <= 1)


def test_predict_worker_reference(cutouts: Cutouts) -> None:
    """Values from FRion 1.1.5 with the bundled CODE maps."""
    ion = predict(cutouts).update._doc["$set"]["frion"]
    assert np.mean(ion["RMs"]) == pytest.approx(-1.19, abs=0.01)
    theta = np.array(ion["theta_real"]) + 1j * np.array(ion["theta_imag"])
    assert np.angle(theta[0]) == pytest.approx(-0.334, abs=0.01)
    assert np.angle(theta[-1]) == pytest.approx(-0.177, abs=0.01)


def test_correct_worker(cutouts: Cutouts) -> None:
    prediction = predict(cutouts)
    update = frion.correct_worker.fn(
        beam=cutouts.beam,
        outdir=cutouts.cutdir,
        field=FIELD,
        prediction=prediction,
        island=cutouts.island,
    )

    q_ion = f"{ISLAND}/{ISLAND}.q.ion.fits"
    u_ion = f"{ISLAND}/{ISLAND}.u.ion.fits"
    assert update._filter == {"Source_ID": ISLAND}
    assert update._doc == {
        "$set": {
            f"beams.{FIELD}.q_file_ion": q_ion,
            f"beams.{FIELD}.u_file_ion": u_ion,
        }
    }

    ion = prediction.update._doc["$set"]["frion"]
    theta = np.array(ion["theta_real"]) + 1j * np.array(ion["theta_imag"])
    q_in = fits.getdata(cutouts.cutdir / ISLAND / f"{ISLAND}.q.fits")
    u_in = fits.getdata(cutouts.cutdir / ISLAND / f"{ISLAND}.u.fits")
    q_out = fits.getdata(cutouts.cutdir / q_ion)
    u_out = fits.getdata(cutouts.cutdir / u_ion)
    assert q_out.shape == q_in.shape
    expected = (q_in + 1j * u_in) / theta[:, None, None, None]
    np.testing.assert_allclose(q_out + 1j * u_out, expected, rtol=1e-5, atol=1e-6)


def test_serial_loop(cutouts: Cutouts) -> None:
    result = frion.serial_loop.fn(
        island=cutouts.island,
        field=FIELD,
        beam=cutouts.beam,
        start_time=cutouts.start_time,
        end_time=cutouts.end_time,
        freq_hz_array=cutouts.freq,
        cutdir=cutouts.cutdir,
        ionex_server="ftp://ftp.aiub.unibe.ch/CODE/",
        ionex_prefix="codg",
        ionex_proxy_server=None,
        ionex_formatter="ftp.aiub.unibe.ch",
        ionex_predownload=False,
    )
    assert isinstance(result, frion.FrionResults)
    assert "frion" in result.prediction.update._doc["$set"]
    assert f"beams.{FIELD}.q_file_ion" in result.correction._doc["$set"]


def test_fix_header_moves_freq_axis() -> None:
    original = make_header()
    original["CTYPE3"], original["CTYPE4"] = "FREQ", "STOKES"
    cutout = make_header()
    fixed = fix_header(cutout, original)
    assert fixed["CTYPE3"] == "FREQ"
    assert fixed["CTYPE4"] == "STOKES"


def test_fix_header_same_axes() -> None:
    header = make_header()
    assert fix_header(header, header) == header
