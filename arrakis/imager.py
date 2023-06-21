#!/usr/bin/env python3
"""Arrkis imager"""
import argparse
import hashlib
import logging
import multiprocessing as mp
import os
import pickle
from glob import glob
from pathlib import Path
from subprocess import CalledProcessError
from typing import Any, Dict, List, NamedTuple, Optional, Tuple, Union

import astropy.units as u
import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.stats import mad_std
from astropy.table import Table
from astropy.wcs import WCS
from casatasks import vishead
from dask import compute, delayed, visualize
from dask.delayed import Delayed
from dask.distributed import Client, LocalCluster
from dask_mpi import initialize
from racs_tools import beamcon_2D
from radio_beam import Beam
from schwimmbad import SerialPool
from spython.main import Client as sclient
from tqdm.auto import tqdm

from arrakis import fix_ms_dir, fix_ms_corrs
from arrakis.logger import logger
from arrakis.utils import (
    beam_from_ms,
    field_idx_from_ms,
    logo_str,
    wsclean,
)


class ImageSet(NamedTuple):
    """Container to organise files related to t he imaging of a measurement set."""

    ms: Path
    prefix: str
    image_lists: Dict[str, List[str]]
    aux_lists: Optional[Dict[Tuple[str, str], List[str]]] = None


def get_wsclean(wsclean: Union[Path, str]) -> Path:
    """Pull wsclean image from dockerhub (or wherver).

    Args:
        version (str, optional): wsclean image tag. Defaults to "3.1".

    Returns:
        Path: Path to wsclean image.
    """
    sclient.load(str(wsclean))
    if isinstance(wsclean, str):
        return Path(sclient.pull(wsclean))
    return wsclean


def cleanup_imageset(purge: bool, image_set: ImageSet) -> None:
    """Delete images associated with an input ImageSet

    Args:
        purge (bool): Whether files will be deleted or skipped.
        image_set (ImageSet): Collection of files that will be removed.
    """
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    if not purge:
        logger.info("Not purging intermediate files")
        return

    for pol, image_list in image_set.image_lists.items():
        logger.critical(f"Removing {pol=} images for {image_set.ms}")
        for image in image_list:
            logger.critical(f"Removing {image}")
            try:
                os.remove(image)
            except FileNotFoundError:
                logger.critical(f"{image} not available for deletion. ")

    # The aux images are the same between the native images and the smoothed images,
    # they were just copied across directly without modification
    if image_set.aux_lists:
        logger.critical(f"Removing auxillary images. ")
        for (pol, aux), aux_list in image_set.aux_lists.items():
            for aux_image in aux_list:
                try:
                    logger.critical(f"Removing {aux_image}")
                    os.remove(aux_image)
                except FileNotFoundError:
                    logger.critical(f"{aux_image} not available for deletion. ")

    return


def get_prefix(
    ms: Path,
    out_dir: Path,
) -> Path:
    """Derive a consistent prefix style from a input MS name.

    Args:
        ms (Path): Path to a Measurement Set that a prefix will be derived from
        out_dir (Path): The final location that wsclean output data will be written to

    Returns:
        Path: The prefix, including the output directory name.
    """
    idx = field_idx_from_ms(ms.resolve(strict=True).as_posix())
    field = vishead(vis=ms.resolve(strict=True).as_posix(), mode="list")["field"][0][
        idx
    ]
    beam = beam_from_ms(ms.resolve(strict=True).as_posix())
    prefix = f"image.{field}.contcube.beam{beam:02}"
    return out_dir / prefix


@delayed()
def image_beam(
    ms: Path,
    field_idx: int,
    out_dir: Path,
    prefix: str,
    simage: Path,
    pols: str = "IQU",
    nchan: int = 36,
    scale: u.Quantity = 2.5 * u.arcsec,
    npix: int = 4096,
    join_polarizations: bool = True,
    join_channels: bool = True,
    squared_channel_joining: bool = True,
    mgain: float = 0.7,
    niter: int = 100_000,
    auto_mask: float = 3,
    force_mask_rounds: Optional[int] = None,
    auto_threshold: float = 1,
    gridder: Optional[str] = None,
    robust: float = -0.5,
    mem: float = 90,
    absmem: Optional[float] = None,
    taper: Optional[float] = None,
    reimage: bool = False,
    minuv_l: float = 0.0,
    parallel_deconvolution: Optional[int] = None,
    nmiter: Optional[int] = None,
    local_rms: bool = False,
    local_rms_window: Optional[float] = None,
    multiscale: bool = False,
    multiscale_scale_bias: Optional[float] = None,
) -> ImageSet:
    """Image a single beam"""
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    if not reimage:
        # Look for existing images
        # NOTE: The current format of the I polarisation is different to the expression below.
        # Raising an error for visibility.
        checkvals = np.array(
            [f"{prefix}-{i:04}-{s}-image.fits" for s in pols for i in range(nchan)]
        )
        checks = np.array([os.path.exists(f) for f in checkvals])
        raise ValueError("The reimage option is not properly supported. ")

    commands = []
    # Do any I cleaning separately
    do_stokes_I = "I" in pols
    if do_stokes_I:
        command = wsclean(
            mslist=[ms.resolve(strict=True).as_posix()],
            use_mpi=False,
            name=prefix,
            pol="I",
            verbose=True,
            channels_out=nchan,
            scale=f"{scale.to(u.arcsec).value}asec",
            size=f"{npix} {npix}",
            join_polarizations=False,  # Only do I
            join_channels=join_channels,
            squared_channel_joining=False,  # Dont want to square I
            mgain=mgain,
            niter=niter,
            auto_mask=auto_mask,
            force_mask_rounds=force_mask_rounds,
            auto_threshold=auto_threshold,
            gridder=gridder,
            weight=f"briggs {robust}",
            log_time=False,
            mem=mem,
            abs_mem=absmem,
            taper_gaussian=f"{taper}asec" if taper else None,
            field=field_idx,
            parallel_deconvolution=parallel_deconvolution,
            minuv_l=minuv_l,
            nmiter=nmiter,
            local_rms=local_rms,
            local_rms_window=local_rms_window,
            multiscale_scale_bias=multiscale_scale_bias,
            multiscale=multiscale,
        )
        commands.append(command)
        pols = pols.replace("I", "")

    if squared_channel_joining:
        logger.info("Using squared channel joining")
        logger.info("Reducing mask by sqrt(2) to account for this")
        auto_mask_reduce = np.round(auto_mask / (np.sqrt(2)), decimals=2)

        logger.info(f"auto_mask = {auto_mask}")
        logger.info(f"auto_mask_reduce = {auto_mask_reduce}")
    else:
        auto_mask_reduce = auto_mask

    if local_rms_window:
        local_rms_window = int(local_rms_window / 2)
        logger.info(f"Scaled local RMS window to {local_rms_window}.")
        
    command = wsclean(
        mslist=[ms.resolve(strict=True).as_posix()],
        use_mpi=False,
        name=prefix,
        pol=pols,
        verbose=True,
        channels_out=nchan,
        scale=f"{scale.to(u.arcsec).value}asec",
        size=f"{npix} {npix}",
        join_polarizations=join_polarizations,
        join_channels=join_channels,
        squared_channel_joining=squared_channel_joining,
        mgain=mgain,
        niter=niter,
        auto_mask=auto_mask_reduce,
        force_mask_rounds=force_mask_rounds,
        auto_threshold=auto_threshold,
        gridder=gridder,
        weight=f"briggs {robust}",
        log_time=False,
        mem=mem,
        abs_mem=absmem,
        taper_gaussian=f"{taper}asec" if taper else None,
        field=field_idx,
        parallel_deconvolution=parallel_deconvolution,
        minuv_l=minuv_l,
        nmiter=nmiter,
        local_rms=local_rms,
        local_rms_window=local_rms_window,
        multiscale=multiscale,
        multiscale_scale_bias=multiscale_scale_bias,
    )
    commands.append(command)

    root_dir = ms.parent

    for command in commands:
        logger.info(f"Running wsclean with command: {command}")
        try:
            output = sclient.execute(
                image=simage.resolve(strict=True).as_posix(),
                command=command.split(),
                bind=f"{out_dir}:{out_dir}, {root_dir.resolve(strict=True).as_posix()}:{root_dir.resolve(strict=True).as_posix()}",
                return_result=True,
                quiet=False,
                stream=True,
            )
            for line in output:
                logger.info(line.rstrip())
                # Catch divergence - look for the string 'KJy' in the output
                if "KJy" in line:
                    raise ValueError(
                        f"Detected divergence in wsclean output: {line.rstrip()}"
                    )
        except CalledProcessError as e:
            logger.error(f"Failed to run wsclean with command: {command}")
            logger.error(f"Stdout: {e.stdout}")
            logger.error(f"Stderr: {e.stderr}")
            logger.error(f"{e=}")
            raise e

    # Check rms of image to check for divergence
    if do_stokes_I:
        pols += "I"
    for pol in pols:
        mfs_image = (
            f"{prefix}-MFS-image.fits"
            if pol == "I"
            else f"{prefix}-MFS-{pol}-image.fits"
        )
        rms = mad_std(fits.getdata(mfs_image), ignore_nan=True)
        if rms > 1:
            # raise ValueError(f"RMS of {rms} is too high in image {mfs_image}, try imaging with lower mgain {mgain - 0.1}")
            logger.error(
                f"RMS of {rms} is too high in image {mfs_image}, try imaging with lower mgain {mgain - 0.1}"
            )

    # Get images
    image_lists = {}
    aux_lists = {}
    for pol in pols:
        imglob = (
            f"{prefix}*[0-9]-image.fits"
            if pol == "I"
            else f"{prefix}*[0-9]-{pol}-image.fits"
        )
        image_list = sorted(glob(imglob))
        image_lists[pol] = image_list

        logger.info(f"Found {len(image_list)} images for {pol=} {ms}.")

        for aux in ["model", "psf", "residual", "dirty"]:
            aux_list = (
                sorted(glob(f"{prefix}*[0-9]-{aux}.fits"))
                if pol == "I" or aux == "psf"
                else sorted(glob(f"{prefix}*[0-9]-{pol}-{aux}.fits"))
            )
            aux_lists[(pol, aux)] = aux_list

            logger.info(f"Found {len(aux_list)} images for {pol=} {aux=} {ms}.")

    logger.info("Constructing ImageSet")
    image_set = ImageSet(
        ms=ms, prefix=prefix, image_lists=image_lists, aux_lists=aux_lists
    )

    logger.debug(f"{image_set=}")

    return image_set


@delayed(nout=2)
def make_cube(
    pol: str, image_set: ImageSet, common_beam_pkl: Path, aux_mode: Optional[str] = None
) -> Tuple[Path, Path]:
    """Make a cube from the images"""
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    logger.info(f"Creating cube for {pol=} {image_set.ms=}")
    image_list = image_set.image_lists[pol]

    image_type = "restored" if aux_mode is None else aux_mode

    # First combine images into cubes
    freqs = []
    rmss = []
    for chan, image in enumerate(
        tqdm(
            image_list,
            desc="Reading channel image",
            leave=False,
        )
    ):
        # init cube
        if chan == 0:
            old_name = image
            old_header = fits.getheader(old_name)
            wcs = WCS(old_header)
            idx = 0
            for j, t in enumerate(
                wcs.axis_type_names[::-1]
            ):  # Reverse to match index order
                if t == "FREQ":
                    idx = j
                    break

            plane_shape = list(fits.getdata(old_name).shape)
            cube_shape = plane_shape.copy()
            cube_shape[idx] = len(image_list)

            data_cube = np.zeros(cube_shape)

            out_dir = os.path.dirname(old_name)
            old_base = os.path.basename(old_name)
            new_base = old_base
            b_idx = new_base.find("beam") + len("beam") + 2
            sub = new_base[b_idx:]
            new_base = new_base.replace(sub, ".conv.fits")
            new_base = new_base.replace("image", f"image.{image_type}.{pol.lower()}")
            new_name = os.path.join(out_dir, new_base)

        plane = fits.getdata(image) # Stokes do NOT need to be scaled so long as fix_ms_corrs applied!
        plane_rms = mad_std(plane, ignore_nan=True)
        rmss.append(plane_rms)
        data_cube[:, chan] = plane
        freq = WCS(image).spectral.pixel_to_world(0)
        freqs.append(freq.to(u.Hz).value)
    # Write out cubes
    freqs = np.array(freqs) * u.Hz
    rmss_arr = np.array(rmss) * u.Jy / u.beam
    assert np.diff(freqs).std() < 1e-6 * u.Hz, "Frequencies are not evenly spaced"
    new_header = old_header.copy()
    new_header["NAXIS"] = len(cube_shape)
    new_header["NAXIS3"] = len(freqs)
    new_header["CRPIX3"] = 1
    new_header["CRVAL3"] = freqs[0].value
    new_header["CDELT3"] = np.diff(freqs).mean().value
    new_header["CUNIT3"] = "Hz"

    tmp_header = new_header.copy()
    # Need to swap NAXIS 3 and 4 to make LINMOS happy - booo
    for a, b in ((3, 4), (4, 3)):
        new_header[f"CTYPE{a}"] = tmp_header[f"CTYPE{b}"]
        new_header[f"CRPIX{a}"] = tmp_header[f"CRPIX{b}"]
        new_header[f"CRVAL{a}"] = tmp_header[f"CRVAL{b}"]
        new_header[f"CDELT{a}"] = tmp_header[f"CDELT{b}"]
        new_header[f"CUNIT{a}"] = tmp_header[f"CUNIT{b}"]

    # Cube is currently STOKES, FREQ, RA, DEC - needs to be FREQ, STOKES, RA, DEC
    data_cube = np.moveaxis(data_cube, 1, 0)

    # Deserialise beam
    with open(common_beam_pkl, "rb") as f:
        common_beam = pickle.load(f)
    new_header = common_beam.attach_to_header(new_header)
    fits.writeto(new_name, data_cube, new_header, overwrite=True)
    logger.info(f"Written {new_name}")

    # Write out weights
    # Must be of the format:
    # #Channel Weight
    # 0 1234.5
    # 1 6789.0
    # etc.
    new_w_name = new_name.replace(
        f"image.{image_type}", f"weights.{image_type}"
    ).replace(".fits", ".txt")
    data = dict(
        Channel=np.arange(len(rmss_arr.value)),
        Weight=rmss_arr.value,
    )
    tab = Table(data)
    tab.write(new_w_name, format="ascii.commented_header", overwrite=True)

    return new_name, new_w_name


@delayed()
def get_beam(image_sets: List[ImageSet], cutoff: Optional[float]) -> Path:
    """Derive a common resolution across all images within a set of ImageSet

    Args:
        image_sets (List[ImageSet]): All input beam ImageSets that a common resolution will be derived for
        cuttoff (float, optional): The maximum major axis of the restoring beam that is allowed when
        searching for the lowest common beam. Images whose restoring beam's major acis is larger than
        this are ignored. Defaults to None.

    Returns:
        Path: Path to the pickled beam object
    """
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    # convert dict to list
    image_list = []
    for image_set in image_sets:
        for _, sub_image_list in image_set.image_lists.items():
            image_list.extend(sub_image_list)

    # Consistent hash between runs
    image_list = sorted(image_list)

    logger.info(f"The length of the image list is: {len(image_list)}")

    # Create a unique hash for the beam log filename
    image_hash = hashlib.md5("".join(image_list).encode()).hexdigest()

    common_beam, _ = beamcon_2D.getmaxbeam(files=image_list, cutoff=cutoff)

    logger.info(f"The common beam is: {common_beam=}")

    # serialise the beam
    common_beam_pkl = Path(f"beam_{image_hash}.pkl")

    with open(common_beam_pkl, "wb") as f:
        logger.info(f"Creating {common_beam_pkl}")
        pickle.dump(common_beam, f)

    return common_beam_pkl


@delayed()
def smooth_imageset(
    image_set: ImageSet,
    common_beam_pkl: Path,
    cutoff: Optional[float] = None,
    aux_mode: Optional[str] = None,
) -> ImageSet:
    """Smooth all images described within an ImageSet to a desired resolution

    Args:
        image_set (ImageSet): Container whose image_list will be convolved to common resolution
        common_beam_pkl (Path): Location of pickle file with beam description
        cutoff (Optional[float], optional): PSF cutoff passed to the beamcon_2D worker. Defaults to None.
        aux_model (Optional[str], optional): The image type in the `aux_lists` property of `image_set` that contains the images to smooth. If
        not set then the `image_lists` property of `image_set` is used. Defaults to None.
    Returns:
        ImageSet: A copy of `image_set` pointing to the smoothed images. Note the `aux_images` property is not carried forward.
    """
    # Smooth image
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    # Deserialise the beam
    with open(common_beam_pkl, "rb") as f:
        logger.info(f"Loading common beam from {common_beam_pkl}")
        common_beam = pickle.load(f)

    logger.info(f"Smooting {image_set.ms} images")

    images_to_smooth: Dict[str, List[str]]
    if aux_mode is None:
        images_to_smooth = image_set.image_lists
    else:
        logger.info(f"Extracting images for {aux_mode=}.")
        assert image_set.aux_lists is not None, f"{image_set=} has empty aux_lists."
        images_to_smooth = {
            pol: images
            for (pol, img_type), images in image_set.aux_lists.items()
            if aux_mode == img_type
        }

    sm_images = {}
    for pol, pol_images in images_to_smooth.items():
        logger.info(f"Smoothing {pol=} for {image_set.ms}")
        for img in pol_images:
            logger.info(f"Smoothing {img}")
            beamcon_2D.worker(
                file=img,
                outdir=None,
                new_beam=common_beam,
                conv_mode="robust",
                suffix="conv",
                cutoff=cutoff,
            )

        sm_images[pol] = [image.replace(".fits", ".conv.fits") for image in pol_images]

    return ImageSet(
        ms=image_set.ms,
        prefix=image_set.prefix,
        image_lists=sm_images,
    )


@delayed()
def cleanup(purge: bool, image_sets: List[ImageSet], ignore_files: List[Any]) -> None:
    """Utility to remove all images described by an collection of ImageSets. Internally
    called `cleanup_imageset`.

    Args:
        purge (bool): Whether files are actually removed or skipped.
        image_sets (List[ImageSet]): Collection of ImageSets that would be deleted
        ignore_files (List[Any]): Collection of items to ignore. Nothing is done with this
        and is purely used to exploit the dask dependency tracking.
    """
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    logger.warn(f"Ignoring files in {ignore_files=}. ")

    if not purge:
        logger.info("Not purging intermediate files")
        return

    for image_set in image_sets:
        cleanup_imageset(purge=purge, image_set=image_set)

    return


@delayed()
def fix_ms(ms: Path) -> Path:
    """Apply the corrections to the FEED table of a measurement set that
    is required for the ASKAP measurement sets.

    Args:
        ms (Path): Path to the measurement set to fix.

    Returns:
        Path: Path to the corrected measurement set.
    """
    fix_ms_dir.main(ms.resolve(strict=True).as_posix())
    return ms

@delayed()
def fix_ms_askap_corrs(ms: Path) -> Path:
    """Applies a correction to raw telescope polarisation products to rotate them
    to the wsclean espected form. This is essentially related to the third-axis of 
    ASKAP and reorientating its 'X' and 'Y's. 

    Args:
        ms (Path): Path to the measurement set to be corrected. 

    Returns:
        Path: Path of the measurementt set containing the corrections.
    """
    logger.info(f"Correcting {str(ms)} correlations for wsclean. ")
    
    fix_ms_corrs.main(ms=ms)

    return ms

def main(
    msdir: Path,
    out_dir: Path,
    cutoff: Optional[float] = None,
    robust: float = -0.5,
    pols: str = "IQU",
    nchan: int = 36,
    size: int = 6074,
    scale: float = 2.5,
    mgain: float = 0.8,
    niter: int = 100_000,
    auto_mask: float = 3,
    force_mask_rounds: Union[int, None] = None,
    auto_threshold: float = 1,
    taper: Union[float, None] = None,
    reimage: bool = False,
    purge: bool = False,
    minuv: float = 0.0,
    parallel_deconvolution: Optional[int] = None,
    gridder: Optional[str] = None,
    nmiter: Optional[int] = None,
    local_rms: bool = False,
    local_rms_window: Optional[float] = None,
    wsclean_path: Union[Path, str] = "docker://alecthomson/wsclean:latest",
    multiscale: Optional[bool] = None,
    multiscale_scale_bias: Optional[float] = None,
    absmem: Optional[float] = None,
    make_residual_cubes: Optional[bool] = False
):
    simage = get_wsclean(wsclean=wsclean_path)

    mslist = sorted(msdir.glob("scienceData*_averaged_cal.leakage.ms"))

    scale = scale * u.arcsecond

    assert (len(mslist) > 0) & (
        len(mslist) == 36
    ), f"Incorrect number of MS files found: {len(mslist)} / 36"

    logger.info(f"Will image {len(mslist)} MS files in {msdir} to {out_dir}")
    cleans = []  # type: List[Delayed]

    # Do this in serial since CASA gets upset
    prefixs = {}
    field_idxs = {}
    for ms in tqdm(mslist, "Getting metadata"):
        prefix = get_prefix(ms, out_dir)
        prefixs[ms] = prefix
        field_idxs[ms] = field_idx_from_ms(ms.resolve(strict=True).as_posix())

    # Image_sets will be a containter that represents the output wsclean image products
    # produced for each beam. A single ImageSet is a container for a single beam.
    image_sets = []
    for ms in mslist:
        logger.info(f"Imaging {ms}")
        # Apply Emil's fix for MSs feed centre
        ms_fix = fix_ms(ms)
        ms_fix = fix_ms_askap_corrs(ms=ms_fix)
        # Image with wsclean
        image_set = image_beam(
            ms=ms_fix,
            field_idx=field_idxs[ms],
            out_dir=out_dir,
            prefix=prefixs[ms].resolve(strict=False).as_posix(),
            simage=simage.resolve(strict=True),
            robust=robust,
            pols=pols,
            nchan=nchan,
            scale=scale,
            npix=size,
            mgain=mgain,
            niter=niter,
            auto_mask=auto_mask,
            force_mask_rounds=force_mask_rounds,
            auto_threshold=auto_threshold,
            taper=taper,
            reimage=reimage,
            minuv_l=minuv,
            parallel_deconvolution=parallel_deconvolution,
            gridder=gridder,
            nmiter=nmiter,
            local_rms=local_rms,
            local_rms_window=local_rms_window,
            multiscale=multiscale,
            multiscale_scale_bias=multiscale_scale_bias,
            absmem=absmem,
        )

        image_sets.append(image_set)

    # Compute the smallest beam that all images can be convolved to.
    # This requires all imaging rounds to be completed, so the total
    # set of ImageSets are first derived before this is called.
    common_beam_pkl = get_beam(image_sets=image_sets, cutoff=cutoff)

    cube_aux_modes = (None, "residual") if make_residual_cubes else (None, )

    # With the final beam each *image* in the ImageSet across IQU are
    # smoothed and then form the cube for each stokes.
    for image_set in image_sets:
        # Per loop containers since we are iterating over image modes
        smooth_image_sets = []
        clean_sm_image_sets = []
        for aux_mode in cube_aux_modes:
            # Smooth the *images* in an ImageSet across all Stokes. This
            # limits the number of workers to 36, i.e. this is operating
            # beamwise
            sm_image_set = smooth_imageset(
                image_set,
                common_beam_pkl=common_beam_pkl,
                cutoff=cutoff,
                aux_mode=aux_mode,
            )
            smooth_image_sets.append(sm_image_set)

            # Make a cube. This is operating across beams and stokes
            cube_images = [
                make_cube(
                    pol=pol,
                    image_set=sm_image_set,
                    common_beam_pkl=common_beam_pkl,
                    aux_mode=aux_mode,
                )
                for pol in pols
            ]

            # Clean up smoothed images files. Note the
            # ignore_files that is used to preserve the
            # dependency between dask tasks
            clean = cleanup(
                purge=purge,
                image_sets=[sm_image_set],
                ignore_files=cube_images,  # To keep the dask dependency tracking
            )
            clean_sm_image_sets.append(clean)

        # Now clean the original output images from wscean
        clean = cleanup(
            purge=purge,
            image_sets=[image_set],
            ignore_files=clean_sm_image_sets,  # To keep the dask dependency tracking
        )
        cleans.append(clean)

    # Trust nothing
    visualize(cleans, filename="compute_graph.pdf", optimize_graph=True, rankdir="LR")

    return compute(cleans)


def imager_parser(parent_parser: bool = False) -> argparse.ArgumentParser:
    """Return the argument parser for the imager routine.

    Args:
        parent_parser (bool, optional): Ensure the parser is configured so it can be added as a parent to a new parser. This will disables the -h/--help action from being generated. Defaults to False.

    Returns:
        argparse.ArgumentParser: Arguments required for the imager routine
    """

    # Help string to be shown using the -h option
    descStr = f"""
    {logo_str}

    {__doc__}
    """

    # Parse the command line options
    img_parser = argparse.ArgumentParser(
        add_help=not parent_parser,
        description=descStr,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser = img_parser.add_argument_group("imaging arguments")

    parser.add_argument(
        "msdir",
        type=Path,
        help="Directory containing MS files",
    )
    parser.add_argument(
        "outdir",
        type=Path,
        help="Directory to output images",
    )
    parser.add_argument(
        "--psf_cutoff",
        type=float,
        help="Cutoff for smoothing in units of arcseconds. ",
    )
    parser.add_argument(
        "--robust",
        type=float,
        default=-0.5,
    )
    parser.add_argument(
        "--nchan",
        type=int,
        default=36,
    )
    parser.add_argument(
        "--pols",
        type=str,
        default="IQU",
    )
    parser.add_argument(
        "--size",
        type=int,
        default=4096,
    )
    parser.add_argument(
        "--scale",
        type=float,
        default=2.5,
    )
    parser.add_argument(
        "--mgain",
        type=float,
        default=0.8,
    )
    parser.add_argument(
        "--niter",
        type=int,
        default=100_000,
    )
    parser.add_argument(
        "--nmiter",
        type=int,
        default=None,
    )
    parser.add_argument(
        "--auto_mask",
        type=float,
        default=3.0,
    )
    parser.add_argument(
        "--auto_threshold",
        type=float,
        default=1.0,
    )
    parser.add_argument(
        "--local_rms",
        action="store_true",
    )
    parser.add_argument(
        "--local_rms_window",
        type=float,
        default=None,
    )
    parser.add_argument(
        "--force_mask_rounds",
        type=int,
        default=None,
    )
    parser.add_argument(
        "--gridder",
        type=str,
        default=None,
        choices=["direct-ft", "idg", "wgridder", "tuned-wgridder", "wstacking"],
    )
    parser.add_argument(
        "--taper",
        type=float,
        default=None,
    )
    parser.add_argument(
        "--minuv",
        type=float,
        default=0.0,
    )
    parser.add_argument(
        "--parallel",
        type=int,
        default=None,
    )
    parser.add_argument(
        "--purge",
        action="store_true",
        help="Purge intermediate files",
    )
    parser.add_argument(
        "--mpi",
        action="store_true",
        help="Use MPI",
    )
    parser.add_argument(
        "--reimage",
        action="store_true",
        help="Force a new round of imaging. Otherwise, will skip if images already exist.",
    )
    parser.add_argument(
        "--multiscale",
        action="store_true",
        help="Use multiscale clean",
    )
    parser.add_argument(
        "--multiscale_scale_bias",
        type=float,
        default=None,
        help="The multiscale scale bias term provided to wsclean. ",
    )
    parser.add_argument(
        "--absmem",
        type=float,
        default=None,
        help="Absolute memory limit in GB",
    )

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--hosted-wsclean",
        type=str,
        default="docker://alecthomson/wsclean:latest",
        help="Docker or Singularity image for wsclean [docker://alecthomson/wsclean:latest]",
    )
    group.add_argument(
        "--local_wsclean",
        type=Path,
        default=None,
        help="Path to local wsclean Singularity image",
    )
    
    group.add_argument(
        "--make_residual_cubes",
        action='store_true',
        help="Create residual cubes as well as cubes from restored images. "
    )

    return img_parser


def cli():
    """Command-line interface"""
    parser = imager_parser()

    args = parser.parse_args()

    if args.mpi:
        initialize(interface="ipogif0")
        cluster = None

    else:
        cluster = LocalCluster(
            n_workers=1,
            threads_per_worker=1,
            # processes=False,
        )

    with Client(cluster) as client:
        logger.debug(f"{cluster=}")
        logger.debug(f"{client=}")
        main(
            msdir=args.msdir,
            out_dir=args.outdir,
            cutoff=args.cutoff,
            robust=args.robust,
            pols=args.pols,
            nchan=args.nchan,
            size=args.size,
            scale=args.scale,
            mgain=args.mgain,
            niter=args.niter,
            auto_mask=args.auto_mask,
            force_mask_rounds=args.force_mask_rounds,
            auto_threshold=args.auto_threshold,
            minuv=args.minuv,
            purge=args.purge,
            taper=args.taper,
            reimage=args.reimage,
            parallel_deconvolution=args.parallel,
            gridder=args.gridder,
            wsclean_path=Path(args.local_wsclean)
            if args.local_wsclean
            else args.hosted_wsclean,
            multiscale=args.multiscale,
        )


if __name__ == "__main__":
    cli()
