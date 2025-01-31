"""Tests for CLI."""

from __future__ import annotations

import subprocess
import unittest


class test_cli(unittest.TestCase):
    def test_cli_init(self):
        """Tests that the CLI `spice_init` runs."""
        res = subprocess.run(["spice_init", "--help"], check=True)
        self.assertEqual(res.returncode, 0)

    def test_cli_process(self):
        """Tests that the CLI `spice_process` runs."""
        res = subprocess.run(["spice_process", "--help"], check=True)
        self.assertEqual(res.returncode, 0)


if __name__ == "__main__":
    unittest.main()
