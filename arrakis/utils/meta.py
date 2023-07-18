#!/usr/bin/env python
"""Generic program utilities"""


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
