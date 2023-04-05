#!/usr/bin/env python3

# This file demonstrates a deconvolution approach similar to simple-deconvolution-example.py,
# but supports multiple frequencies. It thereby also demonstrates fitting the components
# to a requested spectrum.

# run this e.g. with:
# wsclean \
#  -python-deconvolution mf-deconvolution-example.py \
#  -fit-spectral-pol 2 -channels-out 8 -parallel-gridding 8 -join-channels \
#  -niter 1000 -save-first-residual -auto-threshold 3 -mgain 0.8 \
#  -interval 10 11 -size 1024 1024 -scale 1amin \
#  1052736496-averaged.ms/

import numpy as np
import numba as nb
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from typing import NamedTuple
from functools import partial
from time import time

class RMSynthParams(NamedTuple):
    phis: np.ndarray
    lsq: np.ndarray
    lsq_0: float
    ays: np.ndarray
    fwhm: float

def get_rmsynth_params(freqs: np.ndarray, weights: np.ndarray, nsamp: int=10) -> RMSynthParams:
    """Get parameters for RMSynth

    Args:
        freqs (np.ndarray): Frequencies in Hz
        weights (np.ndarray): Weights

    Returns:
        RMSynthParams: RMSynth parameters
    """
    speed_of_light=299792458.0
    lsq = (speed_of_light / freqs)**2
    Delta_lsq = lsq.max() - lsq.min()
    delta_lsq = np.diff(lsq[::-1]).min()

    # Brentjens & de Bruyn (2005)
    delta_phi = 3.8 / Delta_lsq
    phi_max = np.abs(np.sqrt(3) / delta_lsq)

    phi_step = int( np.floor(delta_phi) / nsamp)
    nphi = int(round(abs((phi_max - 0.0) / phi_step)) * 2.0 + 1.0)
    phi_start = - (nphi-1.0) * phi_step / 2.0
    phi_end = + (nphi-1.0) * phi_step / 2.0
    phis = np.linspace(phi_start, phi_end, nphi)

    lsq_0 = np.average(lsq, weights=weights)
    ays = lsq - lsq_0

    return RMSynthParams(
        phis=phis,
        lsq_0=lsq_0,
        ays=ays,
        lsq=lsq,
        fwhm=delta_phi,
    )


@nb.njit(parallel=True, fastmath=True,)
def rmsynth_1d(
        stokes_q: np.ndarray,
        stokes_u: np.ndarray,
        weights: np.ndarray,
        ays: np.ndarray,
        phis: np.ndarray,
) -> np.ndarray:
    """1D RMSynth

    Args:
        stokes_q (np.ndarray): Stokes Q
        stokes_u (np.ndarray): Stokes U
        weights (np.ndarray): Weights
        ays (np.ndarray): Lsq - Lsq_0
        phis (np.ndarray): Faraday depths

    Returns:
        np.ndarray: FDF
    """
    fdf = np.zeros(phis.shape, dtype='complex128')
    stokes_p = stokes_q + 1j * stokes_u
    kay = 1 / np.sum(weights)

    for i in nb.prange(len(phis)):
        phi = phis[i]
        arg = np.exp(-2.0j * phi * ays)
        fdf[i] = kay * (stokes_p * arg).sum()#(axis=0)
    return fdf

@nb.njit(parallel=True, fastmath=True,)
def rmsynth3d(
        stokes_q: np.ndarray,
        stokes_u: np.ndarray,
        weights: np.ndarray,
        ays: np.ndarray,
        phis: np.ndarray,
) -> np.ndarray:
    """1D RMSynth

    Args:
        stokes_q (np.ndarray): Stokes Q [nchan, y, x]
        stokes_u (np.ndarray): Stokes U [nchan, y, x]
        weights (np.ndarray): Weights [nchan]
        ays (np.ndarray): Lsq - Lsq_0 [nchan]
        phis (np.ndarray): Faraday depths [nphi]

    Returns:
        np.ndarray: FDF [nphi, y, x]
    """
    shape = (len(phis), stokes_q.shape[1], stokes_q.shape[2])
    fdf = np.zeros(shape, dtype='complex128')
    stokes_p = stokes_q + 1j * stokes_u
    kay = 1 / np.sum(weights)

    for i in nb.prange(len(phis)):
        phi = phis[i]
        arg = np.exp(-2.0j * phi * ays)
        arg_mat = arg.repeat(stokes_p.shape[1]).repeat(stokes_p.shape[2]).reshape(stokes_p.shape)
        fdf[i] = kay * (stokes_p * arg_mat).sum(axis=0)

    return fdf



def simple_clean(
        phis: np.ndarray,
        fdf: np.ndarray,
        ays: np.ndarray,
        fwhm: float,
    ) -> np.ndarray:
    """Simple clean

    Args:
        phis (np.ndarray): Faraday depths
        fdf (np.ndarray): FDF
        ays (np.ndarray): Lsq - Lsq_0
        fwhm (float): FWHM

    Returns:
        np.ndarray: Spectrum
    """
    def _gauss(x: np.ndarray, amp: float, mu: float, sigma: float) -> np.ndarray:
        """Gaussion function

        Args:
            x (np.ndarray): X values
            amp (float): Maximum value
            mu (float): Mean
            sigma (float): Standard deviation

        Returns:
            np.ndarray: Gaussian data
        """
        return amp * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))

    # Fit to PI
    fdf_p = np.abs(fdf)
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    partial_gauss = partial(_gauss, sigma=sigma)
    # Initial guess
    p0 = [np.nanmax(fdf_p), phis[np.nanargmax(fdf_p)]]
    coeff, var_matrix = curve_fit(partial_gauss, phis, fdf_p, p0=p0)
    peak_p, peak_phi = coeff

    q_interp = interp1d(
        x=phis,
        y=fdf.real,
    )
    u_interp = interp1d(
        x=phis,
        y=fdf.imag,
    )
    cc_scalar = q_interp(peak_phi) + 1j * u_interp(peak_phi)
    cc_vec = np.zeros(len(phis)).astype('complex128')
    # Update the CC at the peak
    cc_vec[np.argmin(np.abs(phis - peak_phi))] = cc_scalar
    quarr = np.sum(cc_vec[:, np.newaxis] * np.exp(2.0j * np.outer(phis, ays)),axis=0)
    spec = np.array([quarr.real, quarr.imag]).T

    return spec


def deconvolve(
        residual: np.ndarray,
        model: np.ndarray,
        psf: np.ndarray,
        meta: dict,
    ):
    if meta.channels == []:
        raise ValueError("No channels in meta")
    nchan, npol, height, width = residual.shape
    print(
        "Python deconvolve() function was called for "
        + f"{width} x {height} x {npol} (npol) x {nchan} (chan) dataset"
    )

    # residual and model are numpy arrays with dimensions nchan x npol x height x width
    # psf is a numpy array with dimensions nchan x height x width.

    # This file demonstrates a very simple deconvolution strategy, which doesn't
    # support multiple polarizations:
    if npol != 2:
        raise NotImplementedError("npol must be 2")
    # If there are channels missing (flagged), they will be set to NaN
    # They're here set to zero to make it easy to calculate the integrated image
    residual = np.nan_to_num(residual)

    # find the largest peak in the integrated image
    freqs = np.array([x.frequency for x in meta.channels])
    weights_3d = np.ones_like(freqs)
    params_3d = get_rmsynth_params(freqs, weights_3d, nsamp=3)
    tick = time()
    fdf_3d = rmsynth3d(
        stokes_q=residual[:, 0],
        stokes_u=residual[:, 1],
        weights=weights_3d,
        ays=params_3d.ays,
        phis=params_3d.phis,
    )
    tock = time()
    print(f"RMSynth3D took {tock - tick:0.2f} seconds")
    integrated_residual = np.sum(np.abs(fdf_3d), axis=0)
    peak_index = np.unravel_index(
        np.argmax(integrated_residual), integrated_residual.shape
    )
    peak_value = integrated_residual[peak_index]

    mgain_threshold = peak_value * (1.0 - meta.mgain)
    first_threshold = np.max(
        [meta.major_iter_threshold, meta.final_threshold, mgain_threshold]
    )

    print(
        f"Starting iteration {meta.iteration_number}, peak={peak_value}, first threshold={first_threshold}"
    )
    while (
        peak_value > first_threshold
        and meta.iteration_number < meta.max_iterations
    ):
        y = peak_index[0]
        x = peak_index[1]
        spectrum_complex = residual[:, :, y, x]
        if meta.iteration_number < 10:
            np.savetxt(f"spectrum_iter_{meta.iteration_number}.txt", spectrum_complex)
        weights = (spectrum_complex.sum(axis=1) > 0).astype(int)

        params = get_rmsynth_params(freqs, weights)
        fdf = rmsynth_1d(
            stokes_q=spectrum_complex[:, 0],
            stokes_u=spectrum_complex[:, 1],
            weights=weights,
            ays=params.ays,
            phis=params.phis,
        )
        model_spectrum = simple_clean(
            phis=params.phis,
            fdf=fdf,
            ays=params.ays,
            fwhm=params.fwhm,
        )
        if meta.iteration_number < 10:
            np.savetxt(f"model_spectrum_iter_{meta.iteration_number}.txt", model_spectrum)
        # Update the model
        model[:, :, y, x] += model_spectrum

        # Subtract the model from the residual
        psf_shift = (y + height // 2, x + width // 2)
        residual[:, 0] = (
            residual[:, 0]
            - np.roll(psf, psf_shift, axis=(1, 2))
            * model_spectrum[:, 0, np.newaxis, np.newaxis]
        )
        residual[:, 1] = (
            residual[:, 1]
            - np.roll(psf, psf_shift, axis=(1, 2))
            * model_spectrum[:, 1, np.newaxis, np.newaxis]
        )

        ######################
        # Update the residual
        tick = time()
        fdf_3d = rmsynth3d(
            stokes_q=residual[:, 0],
            stokes_u=residual[:, 1] ,
            weights=weights_3d,
            ays=params_3d.ays,
            phis=params_3d.phis,
        )
        tock = time()
        print(f"RMSynth3D took {tock - tick:0.2f} seconds")
        integrated_residual = np.sum(np.abs(fdf_3d), axis=0)
        peak_index = np.unravel_index(
            np.argmax(integrated_residual), integrated_residual.shape
        )
        peak_value = integrated_residual[peak_index]

        print(f"{peak_value=}")

        meta.iteration_number = meta.iteration_number + 1

    print(
        f"Stopped after iteration {meta.iteration_number}, peak={peak_value}"
    )

    # Fill a dictionary with values that wsclean expects:
    result = dict()
    result["residual"] = residual
    result["model"] = model
    result["level"] = peak_value
    result["continue"] = (
        peak_value > meta.final_threshold
        and meta.iteration_number < meta.max_iterations
    )

    print("Finished deconvolve()")
    return result
