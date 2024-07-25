"""Pipeline and flow utility functions."""

from __future__ import annotations

import argparse
import base64
import logging
import shlex
import subprocess
import time
import warnings
from pathlib import Path
from uuid import UUID

import astropy.units as u
import dask.array as da
import numpy as np
from astropy.utils.exceptions import AstropyWarning
from dask import distributed
from dask.delayed import Delayed
from dask.distributed import get_client
from distributed.client import futures_of
from distributed.diagnostics.progressbar import ProgressBar
from distributed.utils import LoopRunner
from prefect import task
from prefect.artifacts import create_markdown_artifact
from spectral_cube.utils import SpectralCubeWarning
from tornado.ioloop import IOLoop
from tqdm.auto import tqdm, trange

from arrakis.logger import TqdmToLogger, UltimateHelpFormatter, logger

warnings.filterwarnings(action="ignore", category=SpectralCubeWarning, append=True)
warnings.simplefilter("ignore", category=AstropyWarning)

SUPPORTED_IMAGE_TYPES = ("png",)

TQDM_OUT = TqdmToLogger(logger, level=logging.INFO)

# Help string to be shown using the -h option
logo_str = """
    mmm   mmm   mmm   mmm   mmm
    )-(   )-(   )-(   )-(   )-(
   ( S ) ( P ) ( I ) ( C ) ( E )
   |   | |   | |   | |   | |   |
   |___| |___| |___| |___| |___|
    mmm     mmm     mmm     mmm
    )-(     )-(     )-(     )-(
   ( R )   ( A )   ( C )   ( S )
   |   |   |   |   |   |   |   |
   |___|   |___|   |___|   |___|

"""


# Stolen from Flint
@task(name="Upload image as artifact")
def upload_image_as_artifact_task(
    image_path: Path, description: str | None = None
) -> UUID:
    """Create and submit a markdown artifact tracked by prefect for an input image.

    Currently supporting png formatted images.

    The input image is converted to a base64 encoding, and embedded directly
    within the markdown string. Therefore, be mindful of the image size as this
    is tracked in the postgres database.

    Args:
        image_path (Path): Path to the image to upload
        description (Optional[str], optional): A description passed to the markdown artifact. Defaults to None.

    Returns:
        UUID: Generated UUID of the registered artifact
    """
    image_type = image_path.suffix.replace(".", "")
    assert image_path.exists(), f"{image_path} does not exist"
    assert (
        image_type in SUPPORTED_IMAGE_TYPES
    ), f"{image_path} has type {image_type}, and is not supported. Supported types are {SUPPORTED_IMAGE_TYPES}"

    with image_path.open("rb") as open_image:
        logger.info(f"Encoding {image_path} in base64")
        image_base64 = base64.b64encode(open_image.read()).decode()

    logger.info("Creating markdown tag")
    markdown = f"![{image_path.stem}](data:image/{image_type};base64,{image_base64})"

    logger.info("Registering artifact")
    image_uuid: UUID = create_markdown_artifact(
        markdown=markdown, description=description
    )

    return image_uuid


def workdir_arg_parser(parent_parser: bool = False) -> argparse.ArgumentParser:
    """Create a workdir parser.

    Args:
        parent_parser (bool, optional): If a parent parser. Defaults to False.

    Returns:
        argparse.ArgumentParser: The parser.
    """
    # Parse the command line options
    work_parser = argparse.ArgumentParser(
        add_help=not parent_parser,
        formatter_class=UltimateHelpFormatter,
    )
    parser = work_parser.add_argument_group("workdir arguments")
    parser.add_argument(
        "datadir",
        type=Path,
        help="Directory to create/find full-size images and 'cutout' directory",
    )

    return work_parser


def generic_parser(parent_parser: bool = False) -> argparse.ArgumentParser:
    """Create a generic parser.

    Args:
        parent_parser (bool, optional): If a parent parser. Defaults to False.

    Returns:
        argparse.ArgumentParser: The parser.
    """
    descStr = f"""
    {logo_str}
    Generic pipeline options

    """

    # Parse the command line options
    gen_parser = argparse.ArgumentParser(
        add_help=not parent_parser,
        description=descStr,
        formatter_class=UltimateHelpFormatter,
    )
    parser = gen_parser.add_argument_group("generic arguments")

    parser.add_argument(
        "field", metavar="field", type=str, help="Name of field (e.g. RACS_2132-50)."
    )

    parser.add_argument(
        "--sbid",
        type=int,
        default=None,
        help="SBID of observation.",
    )

    parser.add_argument(
        "-s",
        "--stokes",
        dest="stokeslist",
        nargs="+",
        type=str,
        default=["I", "Q", "U"],
        help="List of Stokes parameters to image",
    )

    parser.add_argument(
        "-e",
        "--epoch",
        type=int,
        default=0,
        help="Epoch of observation.",
    )

    parser.add_argument(
        "-v", dest="verbose", action="store_true", help="Verbose output."
    )
    parser.add_argument(
        "--host",
        metavar="host",
        type=str,
        default=None,
        help="Host of mongodb (probably $hostname -i).",
    )
    parser.add_argument(
        "--username", type=str, default=None, help="Username of mongodb."
    )

    parser.add_argument(
        "--password", type=str, default=None, help="Password of mongodb."
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Limit the number of islands to process.",
    )
    parser.add_argument(
        "--database", dest="database", action="store_true", help="Add data to MongoDB."
    )

    return gen_parser


def inspect_client(
    client: distributed.Client | None = None,
) -> tuple[str, int, int, u.Quantity, int, u.Quantity]:
    """_summary_.

    Args:
        client (Union[distributed.Client,None]): Dask client to inspect.
            if None, will use the default client.

    Returns:
        Tuple[ str, int, int, u.Quantity, float, u.Quantity ]: addr, nworkers,
            nthreads, memory, threads_per_worker, memory_per_worker
    """
    """Inspect a client"""
    if client is None:
        client = get_client()
    logger.debug(f"Client: {client}")
    info = client._scheduler_identity
    addr = info.get("address")
    workers = info.get("workers", {})
    nworkers = len(workers)
    nthreads = sum(w["nthreads"] for w in workers.values())
    memory = sum([w["memory_limit"] for w in workers.values()]) * u.byte
    threads_per_worker = nthreads // nworkers
    memory_per_worker = memory / nworkers
    return addr, nworkers, nthreads, memory, threads_per_worker, memory_per_worker


def chunk_dask(
    outputs: list,
    batch_size: int = 10_000,
    task_name="",
    progress_text="",
    verbose=True,
) -> list:
    """Run a task in chunks.

    Args:
        outputs (list): List of outputs to chunk
        batch_size (int, optional): Chunk size. Defaults to 10_000.
        task_name (str, optional): Name of task. Defaults to "".
        progress_text (str, optional): Description of task. Defaults to "".
        verbose (bool, optional): Verbose output. Defaults to True.

    Returns:
        list: Completed futures
    """
    client = get_client()
    chunk_outputs = []
    for i in trange(
        0, len(outputs), batch_size, desc=f"Chunking {task_name}", disable=(not verbose)
    ):
        outputs_chunk = outputs[i : i + batch_size]
        futures = client.persist(outputs_chunk)
        # dumb solution for https://github.com/dask/distributed/issues/4831
        if i == 0:
            logger.debug("I sleep!")
            time.sleep(10)
            logger.debug("I awake!")
        tqdm_dask(futures, desc=progress_text, disable=(not verbose), file=TQDM_OUT)
        chunk_outputs.extend(futures)
    return chunk_outputs


def delayed_to_da(list_of_delayed: list[Delayed], chunk: int | None = None) -> da.Array:
    """Convert list of delayed arrays to a dask array.

    Args:
        list_of_delayed (List[delayed]): List of delayed objects
        chunk (int, optional): Chunksize to use. Defaults to None.

    Returns:
        da.Array: Dask array
    """
    sample = list_of_delayed[0].compute()
    dim = (len(list_of_delayed), *sample.shape)
    c_dim = dim if chunk is None else (chunk, *sample.shape)
    darray_list = [
        da.from_delayed(lazy, dtype=sample.dtype, shape=sample.shape)
        for lazy in list_of_delayed
    ]
    return da.stack(darray_list, axis=0).reshape(dim).rechunk(c_dim)


# stolen from https://github.com/tqdm/tqdm/issues/278
class TqdmProgressBar(ProgressBar):
    """Tqdm for Dask."""

    def __init__(
        self,
        keys,
        scheduler=None,
        interval="100ms",
        loop=None,
        complete=True,
        start=True,
        **tqdm_kwargs,
    ):
        """Make a Tqdm progress bar.

        Args:
            keys (Any): Iterable of keys to track
            scheduler (Any | None, optional): scheduler. Defaults to None.
            interval (str, optional): update interval. Defaults to "100ms".
            loop (Any | None, optional): Loop. Defaults to None.
            complete (bool, optional): Complete. Defaults to True.
            start (bool, optional): Start. Defaults to True.

        kwargs:
            **tqdm_kwargs: Tqdm keyword arguments
        """
        super().__init__(keys, scheduler, interval, complete)
        self.tqdm = tqdm(keys, **tqdm_kwargs)
        self.loop = loop or IOLoop()

        if start:
            loop_runner = LoopRunner(self.loop)
            loop_runner.run_sync(self.listen)

    def _draw_bar(self, remaining, all, **kwargs):
        _ = kwargs
        update_ct = (all - remaining) - self.tqdm.n
        self.tqdm.update(update_ct)

    def _draw_stop(self, **kwargs):
        _ = kwargs
        self.tqdm.close()


def tqdm_dask(futures_in: distributed.Future, **kwargs) -> None:
    """Tqdm for Dask futures."""
    futures = futures_of(futures_in)
    if not isinstance(futures, (set, list)):
        futures = [futures]
    TqdmProgressBar(futures, **kwargs)


def port_forward(port: int, target: str) -> None:
    """Forward ports to local host.

    Args:
        port (int): port to forward
        target (str): Target host
    """
    logger.info(f"Forwarding {port} from localhost to {target}")
    cmd = f"ssh -N -f -R {port}:localhost:{port} {target}"
    command = shlex.split(cmd)
    _ = subprocess.Popen(command)


def cpu_to_use(max_cpu: int, count: int) -> int:
    """Find number of cpus to use.

    Find the right number of cpus to use when dividing up a task, such
    that there are no remainders.

    Args:
        max_cpu (int): Maximum number of cores to use for a process.
        count (int): Number of tasks.

    Returns:
        Maximum number of cores to be used that divides into the number

    """
    factors = []
    for i in range(1, count + 1):
        if count % i == 0:
            factors.append(i)
    factors_arr = np.array(factors)
    return np.max(factors_arr[factors_arr <= max_cpu])
