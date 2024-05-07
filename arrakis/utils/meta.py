#!/usr/bin/env python
"""Generic program utilities"""

import importlib
import warnings
from itertools import zip_longest

from astropy.utils.exceptions import AstropyWarning
import numpy as np
from spectral_cube.utils import SpectralCubeWarning

warnings.filterwarnings(action="ignore", category=SpectralCubeWarning, append=True)
warnings.simplefilter("ignore", category=AstropyWarning)


# From https://stackoverflow.com/questions/58065055/floor-and-ceil-with-number-of-decimals#:~:text=The%20function%20np.,a%20number%20with%20zero%20decimals.
def my_ceil(a, precision=0):
    return np.true_divide(np.ceil(a * 10**precision), 10**precision)


def my_floor(a, precision=0):
    return np.true_divide(np.floor(a * 10**precision), 10**precision)


# From https://stackoverflow.com/questions/1176136/convert-string-to-python-class-object
def class_for_name(module_name: str, class_name: str) -> object:
    """Returns a class object given a module name and class name

    Args:
        module_name (str): Module name
        class_name (str): Class name

    Returns:
        object: Class object
    """
    # load the module, will raise ImportError if module cannot be loaded
    m = importlib.import_module(module_name)
    # get the class, will raise AttributeError if class cannot be found
    c = getattr(m, class_name)
    return c


# stolen from https://stackoverflow.com/questions/32954486/zip-iterators-asserting-for-equal-length-in-python
def zip_equal(*iterables):
    sentinel = object()
    for combo in zip_longest(*iterables, fillvalue=sentinel):
        if sentinel in combo:
            raise ValueError("Iterables have different lengths")
        yield combo


def yes_or_no(question: str) -> bool:
    """Ask a yes or no question via input()

    Args:
        question (str): Question to ask

    Returns:
        bool: True for yes, False for no
    """
    while "Please answer 'y' or 'n'":
        reply = str(input(question + " (y/n): ")).lower().strip()
        if reply[:1] == "y":
            return True
        elif reply[:1] == "n":
            return False
        else:
            raise ValueError("Please answer 'y' or 'n'")
