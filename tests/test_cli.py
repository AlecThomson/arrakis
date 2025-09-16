"""Tests for CLI."""

from __future__ import annotations

import subprocess as sp


def run_cli_help(cmd: str) -> bool:
    res = sp.run([cmd, "--help"], check=True)
    return res.returncode == 0


def test_cli_init() -> None:
    assert run_cli_help("spice_init")


def test_cli_process() -> None:
    assert run_cli_help("spice_process")


def test_cli_region() -> None:
    assert run_cli_help("spice_region")


def test_cli_cat() -> None:
    assert run_cli_help("spice_cat")


def test_cli_image() -> None:
    assert run_cli_help("spice_image")
