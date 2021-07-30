"""Tests for importing modules."""

import unittest

class test_imports(unittest.TestCase):
    def test_imports(self):
        """Tests that package imports are working correctly."""
        # This is a bit of a weird test, but package imports
        # have not worked before.
        modules = [
            'spiceracs.cleanup',
            'spiceracs.columns_possum',
            'spiceracs.cutout',
            'spiceracs.get_seps',
            'spiceracs.init_database',
            'spiceracs.linmos',
            'spiceracs.makecat',
            'spiceracs.processSPICE',
            'spiceracs.rmclean_oncuts',
            'spiceracs.rmsynth_oncuts',
            'spiceracs.utils',
        ]
        for module in modules:
            __import__(module)


if __name__ == '__main__':
    unittest.main()