"""Tests for CLI."""

import subprocess
import unittest


class test_cli(unittest.TestCase):
    def test_cli_init(self):
        """Tests that the CLI `initSPICE` runs."""
        res = subprocess.run(['initSPICE', '--help'], check=True)
        self.assertEqual(res.returncode, 0)

    def test_cli_processSPICE(self):
        """Tests that the CLI `processSPICE` runs."""
        res = subprocess.run(['processSPICE', '--help'], check=True)
        self.assertEqual(res.returncode, 0)


if __name__ == '__main__':
    unittest.main()