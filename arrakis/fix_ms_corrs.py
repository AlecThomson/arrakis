#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fix the correlation rotation of ASKAP MSs.

Converts the ASKAP standard correlations to the 'standard' correlations
This will make them compatible with most imagers (e.g. wsclean, CASA)

The old correlations are moved to a new column called 'DATA_ASKAP'
"""

import logging
from pathlib import Path

import astropy.units as u
import numpy as np
from casacore.tables import makecoldesc, table
from tqdm import tqdm

logger = logging.getLogger(__name__)


def get_pol_axis(ms: Path) -> u.Quantity:
    """
    Get the polarization axis from the ASKAP MS
    """
    with table((ms / "FEED").as_posix(), readonly=True, ack=False) as tf:
        ms_feed = tf.getcol("RECEPTOR_ANGLE") * u.rad
        pol_axes = -(ms_feed - 45.0 * u.deg)

    assert (ms_feed[:, 0] == ms_feed[0, 0]).all() & (
        ms_feed[:, 1] == ms_feed[0, 1]
    ).all(), "The RECEPTOR_ANGLE changes with time, please check the MS"
    return pol_axes[0, 0].to(u.deg)


def convert_correlations(correlations: np.ndarray, pol_axis: u.Quantity):
    """
    Convert ASKAP standard correlations to the 'standard' correlations

    Args:
        correlations (np.ndarray): The correlations from the MS. Has shape (NCOR, NCHAN, 4)
        pol_axis (float): The polarization axis angle of the MS

    Returns:
        np.ndarray: The correlations in the 'standard' format

    NOTES:
    In general, ASKAP forms Stokes I, Q, U, V as follows:
    ⎡I⎤   ⎡    1         0         0          1    ⎤ ⎡XX_a⎤
    ⎢ ⎥   ⎢                                        ⎥ ⎢    ⎥
    ⎢Q⎥   ⎢sin(2⋅θ)   cos(2⋅θ)  cos(2⋅θ)  -sin(2⋅θ)⎥ ⎢XY_a⎥
    ⎢ ⎥ = ⎢                                        ⎥⋅⎢    ⎥
    ⎢U⎥   ⎢-cos(2⋅θ)  sin(2⋅θ)  sin(2⋅θ)  cos(2⋅θ) ⎥ ⎢YX_a⎥
    ⎢ ⎥   ⎢                                        ⎥ ⎢    ⎥
    ⎣V⎦   ⎣    0       -1.0⋅i    1.0⋅i        0    ⎦ ⎣YY_a⎦

    Where theta is the polarization axis angle. In the common case of PA=-45deg -> theta=0deg, this becomes:
    ⎡I⎤   ⎡1     0       0    1⎤ ⎡XX_a⎤
    ⎢ ⎥   ⎢                    ⎥ ⎢    ⎥
    ⎢Q⎥   ⎢0     1       1    0⎥ ⎢XY_a⎥
    ⎢ ⎥ = ⎢                    ⎥⋅⎢    ⎥
    ⎢U⎥   ⎢-1    0       0    1⎥ ⎢YX_a⎥
    ⎢ ⎥   ⎢                    ⎥ ⎢    ⎥
    ⎣V⎦   ⎣0   -1.0⋅i  1.0⋅i  0⎦ ⎣YY_a⎦
                or
    ⎡I⎤   ⎡    XX_a + YY_a     ⎤
    ⎢ ⎥   ⎢                    ⎥
    ⎢Q⎥   ⎢    XY_a + YX_a     ⎥
    ⎢ ⎥ = ⎢                    ⎥
    ⎢U⎥   ⎢    -XX_a + YY_a    ⎥
    ⎢ ⎥   ⎢                    ⎥
    ⎣V⎦   ⎣-i⋅XY_a + 1.0⋅i⋅YX_a⎦

    However, most imagers (e.g. wsclean, CASA) expect
    ⎡I⎤   ⎡0.5    0       0    0.5 ⎤ ⎡XX_w⎤
    ⎢ ⎥   ⎢                        ⎥ ⎢    ⎥
    ⎢Q⎥   ⎢0.5    0       0    -0.5⎥ ⎢XY_w⎥
    ⎢ ⎥ = ⎢                        ⎥⋅⎢    ⎥
    ⎢U⎥   ⎢ 0    0.5     0.5    0  ⎥ ⎢YX_w⎥
    ⎢ ⎥   ⎢                        ⎥ ⎢    ⎥
    ⎣V⎦   ⎣ 0   -0.5⋅i  0.5⋅i   0  ⎦ ⎣YY_w⎦
                    or
    ⎡I⎤   ⎡  0.5⋅XX_w + 0.5⋅YY_w   ⎤
    ⎢ ⎥   ⎢                        ⎥
    ⎢Q⎥   ⎢  0.5⋅XX_w - 0.5⋅YY_w   ⎥
    ⎢ ⎥ = ⎢                        ⎥
    ⎢U⎥   ⎢  0.5⋅XY_w + 0.5⋅YX_w   ⎥
    ⎢ ⎥   ⎢                        ⎥
    ⎣V⎦   ⎣-0.5⋅i⋅XY_w + 0.5⋅i⋅YX_w⎦

    To convert between the two, we can use the following matrix:
    ⎡XX_w⎤   ⎡sin(2.0⋅θ) + 1    cos(2.0⋅θ)      cos(2.0⋅θ)    1 - sin(2.0⋅θ)⎤ ⎡XX_a⎤
    ⎢    ⎥   ⎢                                                              ⎥ ⎢    ⎥
    ⎢XY_w⎥   ⎢ -cos(2.0⋅θ)    sin(2.0⋅θ) + 1  sin(2.0⋅θ) - 1    cos(2.0⋅θ)  ⎥ ⎢XY_a⎥
    ⎢    ⎥ = ⎢                                                              ⎥⋅⎢    ⎥
    ⎢YX_w⎥   ⎢ -cos(2.0⋅θ)    sin(2.0⋅θ) - 1  sin(2.0⋅θ) + 1    cos(2.0⋅θ)  ⎥ ⎢YX_a⎥
    ⎢    ⎥   ⎢                                                              ⎥ ⎢    ⎥
    ⎣YY_w⎦   ⎣1 - sin(2.0⋅θ)   -cos(2.0⋅θ)     -cos(2.0⋅θ)    sin(2.0⋅θ) + 1⎦ ⎣YY_a⎦
    Where _w is the 'wsclean' format and _a is the 'ASKAP' format.

    In the case of PA=-45deg -> theta=0deg, this becomes:
    ⎡XX_w⎤   ⎡1   1   1   1⎤ ⎡XX_a⎤
    ⎢    ⎥   ⎢             ⎥ ⎢    ⎥
    ⎢XY_w⎥   ⎢-1  1   -1  1⎥ ⎢XY_a⎥
    ⎢    ⎥ = ⎢             ⎥⋅⎢    ⎥
    ⎢YX_w⎥   ⎢-1  -1  1   1⎥ ⎢YX_a⎥
    ⎢    ⎥   ⎢             ⎥ ⎢    ⎥
    ⎣YY_w⎦   ⎣1   -1  -1  1⎦ ⎣YY_a⎦
                or
    ⎡XX_w⎤   ⎡XX_a + XY_a + YX_a + YY_a ⎤
    ⎢    ⎥   ⎢                          ⎥
    ⎢XY_w⎥   ⎢-XX_a + XY_a - YX_a + YY_a⎥
    ⎢    ⎥ = ⎢                          ⎥
    ⎢YX_w⎥   ⎢-XX_a - XY_a + YX_a + YY_a⎥
    ⎢    ⎥   ⎢                          ⎥
    ⎣YY_w⎦   ⎣XX_a - XY_a - YX_a + YY_a ⎦


    """
    theta = (pol_axis + 45.0 * u.deg).to(u.rad).value
    correction_matrix = np.matrix(
        [
            [
                np.sin(2.0 * theta) + 1,
                np.cos(2.0 * theta),
                np.cos(2.0 * theta),
                1 - np.sin(2.0 * theta),
            ],
            [
                -np.cos(2.0 * theta),
                np.sin(2.0 * theta) + 1,
                np.sin(2.0 * theta) - 1,
                np.cos(2.0 * theta),
            ],
            [
                -np.cos(2.0 * theta),
                np.sin(2.0 * theta) - 1,
                np.sin(2.0 * theta) + 1,
                np.cos(2.0 * theta),
            ],
            [
                1 - np.sin(2.0 * theta),
                -np.cos(2.0 * theta),
                -np.cos(2.0 * theta),
                np.sin(2.0 * theta) + 1,
            ],
        ]
    )
    # This is a matrix multiplication broadcasted along the time and channel axes
    return np.einsum("ij,...j->...i", correction_matrix, correlations)


def get_data_chunk(ms: Path, chunksize: int, data_column: str = "DATA_ASKAP"):
    with table(ms.as_posix(), readonly=True, ack=False) as tab:
        data = tab.__getattr__(data_column)
        for i in range(0, len(data), chunksize):
            yield np.array(data[i : i + chunksize])


def get_nchunks(ms: Path, chunksize: int, data_column: str = "DATA_ASKAP"):
    with table(ms.as_posix(), readonly=True, ack=False) as tab:
        return int(np.ceil(len(tab.__getattr__(data_column)) / chunksize))


def main(
    ms: Path,
    chunksize: int = 10_000,
    data_column: str = "DATA",
):
    """
    Fix the correlations in the MS
    """
    # Open the MS, move the 'data_column' column to DATA_ASKAP
    try:
        with table(ms.as_posix(), readonly=False, ack=False) as tab:
            tab.renamecol(data_column, "DATA_ASKAP")
            tab.flush()
    except RuntimeError:
        logger.warning(
            f"""'DATA_ASKAP' column already exists in {ms}. Have you already run this?
            Using the existing 'DATA_ASKAP' column and continuing...
            """
        )

    # Get the polarization axis
    pol_axis = get_pol_axis(ms)
    logger.info(f"Polarization axis is {pol_axis}")

    # Get the data chunk by chunk and convert the correlations
    # then write them back to the MS in the 'data_column' column
    data_chunks = get_data_chunk(ms, chunksize, data_column="DATA_ASKAP")
    nchunks = get_nchunks(ms, chunksize, data_column="DATA_ASKAP")
    start_row = 0
    with table(ms.as_posix(), readonly=False, ack=False) as tab:
        desc = makecoldesc(data_column, tab.getcoldesc("DATA_ASKAP"))
        try:
            tab.addcols(desc)
        except RuntimeError:
            # Shouldn't ever happen...
            # Putting this here for interactive use when you might muck around with the MS
            logger.warning(f"Column {data_column} already exists in {ms}")
            pass
        for data_chunk in tqdm(data_chunks, total=nchunks):
            data_chunk_cor = convert_correlations(
                data_chunk,
                pol_axis,
            )
            tab.putcol(
                data_column,
                data_chunk_cor,
                startrow=start_row,
                nrow=len(data_chunk_cor),
            )
            tab.flush()
            start_row += len(data_chunk_cor)


def cli():
    import argparse

    parser = argparse.ArgumentParser(description="Fix the correlations in the MS")
    parser.add_argument("ms", type=str, help="The MS to fix")
    parser.add_argument(
        "--chunksize",
        type=int,
        default=1000,
        help="The chunksize to use when reading the MS",
    )
    parser.add_argument(
        "--data-column", type=str, default="DATA", help="The column to fix"
    )
    args = parser.parse_args()
    main(Path(args.ms), args.chunksize, args.data_column)


if __name__ == "__main__":
    cli()
