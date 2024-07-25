"""Plotting utilities."""

from __future__ import annotations

import warnings

from astropy.utils.exceptions import AstropyWarning
from spectral_cube.utils import SpectralCubeWarning

from arrakis.logger import logger

warnings.filterwarnings(action="ignore", category=SpectralCubeWarning, append=True)
warnings.simplefilter("ignore", category=AstropyWarning)


def latexify(fig_width=None, fig_height=None, columns=1):
    """Set up matplotlib's RC params for LaTeX plotting.

    Call this before plotting a figure.

    Args:
        fig_width (float, optional): Figure width. Defaults to None.
        fig_height (float, optional): Figure height. Defaults to None.
        columns (int, optional): Number of columns. Defaults

    """
    from math import sqrt

    import matplotlib as mpl

    # code adapted from http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples
    # Width and max height in inches for IEEE journals taken from
    # computer.org/cms/Computer.org/Journal%20templates/transactions_art_guide.pdf

    assert columns in [1, 2]

    if fig_width is None:
        fig_width = 3.39 if columns == 1 else 6.9  # width in inches

    if fig_height is None:
        golden_mean = (sqrt(5) - 1.0) / 2.0  # Aesthetic ratio
        fig_height = fig_width * golden_mean  # height in inches

    MAX_HEIGHT_INCHES = 8.0
    if fig_height > MAX_HEIGHT_INCHES:
        logger.waning(
            "fig_height too large:"
            + fig_height
            + "so will reduce to"
            + MAX_HEIGHT_INCHES
            + "inches."
        )
        fig_height = MAX_HEIGHT_INCHES

    params = {
        "backend": "pdf",
        "axes.labelsize": 8,  # fontsize for x and y labels (was 10)
        "axes.titlesize": 8,
        "font.size": 8,  # was 10
        "legend.fontsize": 8,  # was 10
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "text.usetex": False,
        "figure.figsize": [fig_width, fig_height],
        "font.family": "serif",
    }

    mpl.rcParams.update(params)
