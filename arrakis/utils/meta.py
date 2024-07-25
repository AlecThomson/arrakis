"""Generic program utilities."""

from __future__ import annotations

import importlib
import warnings
from itertools import zip_longest

import numpy as np
from astropy.utils.exceptions import AstropyWarning
from spectral_cube.utils import SpectralCubeWarning

warnings.filterwarnings(action="ignore", category=SpectralCubeWarning, append=True)
warnings.simplefilter("ignore", category=AstropyWarning)


# From https://stackoverflow.com/questions/58065055/floor-and-ceil-with-number-of-decimals#:~:text=The%20function%20np.,a%20number%20with%20zero%20decimals.
def my_ceil[T](a: T, precision=0) -> T:
    """Ceil a number to a given precision.

    Args:
        a (T): A numeric value to ceil
        precision (int, optional): Precision of ceil. Defaults to 0.

    Returns:
        T: The ceil of a number
    """
    return np.true_divide(np.ceil(a * 10**precision), 10**precision)


def my_floor[T](a: T, precision=0) -> T:
    """Floor a number to a given precision.

    Args:
        a (T): A numeric value to floor
        precision (int, optional): Precision of floor. Defaults to 0.

    Returns:
        T: The floor of a number.
    """
    return np.true_divide(np.floor(a * 10**precision), 10**precision)


# From https://stackoverflow.com/questions/1176136/convert-string-to-python-class-object
def class_for_name(module_name: str, class_name: str) -> object:
    """Returns a class object given a module name and class name.

    Args:
        module_name (str): Module name
        class_name (str): Class name

    Returns:
        object: Class object
    """
    # load the module, will raise ImportError if module cannot be loaded
    m = importlib.import_module(module_name)
    # get the class, will raise AttributeError if class cannot be found
    return getattr(m, class_name)


# stolen from https://stackoverflow.com/questions/32954486/zip-iterators-asserting-for-equal-length-in-python
def zip_equal(*iterables):
    """Zip iterables and assert they are the same length.

    Args:
        *iterables: Iterables to zip

    Yields:
        tuple: Zipped iterables

    Raises:
        ValueError: If iterables have different lengths
    """
    sentinel = object()
    for combo in zip_longest(*iterables, fillvalue=sentinel):
        if sentinel in combo:
            msg = "Iterables have different lengths"
            raise ValueError(msg)
        yield combo


def yes_or_no(question: str) -> bool:
    """Ask a yes or no question via input().

    Args:
        question (str): Question to ask

    Returns:
        bool: True for yes, False for no
    """
    while "Please answer 'y' or 'n'":
        reply = str(input(question + " (y/n): ")).lower().strip()
        if reply[:1] == "y":
            return True
        if reply[:1] == "n":
            return False

        msg = "Please answer 'y' or 'n'"
        raise ValueError(msg)
    return None
