#!/usr/bin/env python
"""Tests for CLI."""

import subprocess
import unittest

class test_cli(unittest.TestCase):
    def test_cli_spicecutout(self):
        """Tests that the CLI `spicecutout` runs."""
        res = subprocess.run(['spicecutout', '--help'])
        self.assertEqual(res.returncode, 0)

    def test_cli_spiceunresolved(self):
        """Tests that the CLI `spiceunresolved` runs."""
        res = subprocess.run(['spiceunresolved', '--help'])
        self.assertEqual(res.returncode, 0)

    def test_cli_spicepolfind(self):
        """Tests that the CLI `spicemoments` runs."""
        res = subprocess.run(['spicemoments', '--help'])
        self.assertEqual(res.returncode, 0)

if __name__ == '__main__':
    unittest.main()
