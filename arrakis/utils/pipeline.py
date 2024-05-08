#!/usr/bin/env python
"""Pipeline and flow utility functions"""

import argparse
import base64
import logging
import shlex
import subprocess
import time
import warnings
from pathlib import Path
from typing import List, Optional, Tuple, Union
from uuid import UUID

import astropy.units as u
import dask.array as da
import dask.distributed as distributed
import numpy as np
from astropy.utils.exceptions import AstropyWarning
from dask.delayed import Delayed
from dask.distributed import get_client
from distributed.client import futures_of
from distributed.diagnostics.progressbar import ProgressBar
from distributed.utils import LoopRunner
from prefect import task
from prefect.artifacts import create_markdown_artifact
from prefect_dask import get_dask_client
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
    image_path: Path, description: Optional[str] = None
) -> UUID:
    """Create and submit a markdown artifact tracked by prefect for an
    input image. Currently supporting png formatted images.

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

    with open(image_path, "rb") as open_image:
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


class performance_report_prefect:
    """Gather performance report from prefect_dask

    Basically stolen from:
        https://distributed.dask.org/en/latest/_modules/distributed/client.html#performance_report

    This creates a static HTML file that includes many of the same plots of the
    dashboard for later viewing.

    The resulting file uses JavaScript, and so must be viewed with a web
    browser.  Locally we recommend using ``python -m http.server`` or hosting
    the file live online.

    Parameters
    ----------
    filename: str, optional
        The filename to save the performance report locally

    stacklevel: int, optional
        The code execution frame utilized for populating the Calling Code section
        of the report. Defaults to `1` which is the frame calling ``performance_report_prefect``

    mode: str, optional
        Mode parameter to pass to :func:`bokeh.io.output.output_file`. Defaults to ``None``.

    storage_options: dict, optional
         Any additional arguments to :func:`fsspec.open` when writing to a URL.

    Examples
    --------
    >>> with performance_report_prefect(filename="myfile.html", stacklevel=1):
    ...     x.compute()
    """

    def __init__(
        self, filename="dask-report.html", stacklevel=1, mode=None, storage_options=None
    ):
        self.filename = filename
        # stacklevel 0 or less - shows dask internals which likely isn't helpful
        self._stacklevel = stacklevel if stacklevel > 0 else 1
        self.mode = mode
        self.storage_options = storage_options or {}

    async def __aenter__(self):
        self.start = time.time()
        with get_dask_client() as client:
            self.last_count = client.run_on_scheduler(
                lambda dask_scheduler: dask_scheduler.monitor.count
            )
            client.get_task_stream(start=0, stop=0)  # ensure plugin

    async def __aexit__(self, exc_type, exc_value, traceback, code=None):
        import fsspec

        with get_dask_client() as client:
            if code is None:
                code = client._get_computation_code(self._stacklevel + 1)
            data = await client.scheduler.performance_report(
                start=self.start, last_count=self.last_count, code=code, mode=self.mode
            )
            with fsspec.open(
                self.filename, mode="w", compression="infer", **self.storage_options
            ) as f:
                f.write(data)

    def __enter__(self):
        with get_dask_client() as client:
            client.sync(self.__aenter__)

    def __exit__(self, exc_type, exc_value, traceback):
        with get_dask_client() as client:
            code = client._get_computation_code(self._stacklevel + 1)
            client.sync(self.__aexit__, exc_type, exc_value, traceback, code=code)


def inspect_client(
    client: Union[distributed.Client, None] = None,
) -> Tuple[str, int, int, u.Quantity, int, u.Quantity]:
    """_summary_

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


def delayed_to_da(
    list_of_delayed: List[Delayed], chunk: Union[int, None] = None
) -> da.Array:
    """Convert list of delayed arrays to a dask array

    Args:
        list_of_delayed (List[delayed]): List of delayed objects
        chunk (int, optional): Chunksize to use. Defaults to None.

    Returns:
        da.Array: Dask array
    """
    sample = list_of_delayed[0].compute()
    dim = (len(list_of_delayed),) + sample.shape
    if chunk is None:
        c_dim = dim
    else:
        c_dim = (chunk,) + sample.shape
    darray_list = [
        da.from_delayed(lazy, dtype=sample.dtype, shape=sample.shape)
        for lazy in list_of_delayed
    ]
    darray = da.stack(darray_list, axis=0).reshape(dim).rechunk(c_dim)

    return darray


# stolen from https://github.com/tqdm/tqdm/issues/278
class TqdmProgressBar(ProgressBar):
    """Tqdm for Dask"""

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
        super(TqdmProgressBar, self).__init__(keys, scheduler, interval, complete)
        self.tqdm = tqdm(keys, **tqdm_kwargs)
        self.loop = loop or IOLoop()

        if start:
            loop_runner = LoopRunner(self.loop)
            loop_runner.run_sync(self.listen)

    def _draw_bar(self, remaining, all, **kwargs):
        update_ct = (all - remaining) - self.tqdm.n
        self.tqdm.update(update_ct)

    def _draw_stop(self, **kwargs):
        self.tqdm.close()


def tqdm_dask(futures_in: distributed.Future, **kwargs) -> None:
    """Tqdm for Dask futures"""
    futures = futures_of(futures_in)
    if not isinstance(futures, (set, list)):
        futures = [futures]
    TqdmProgressBar(futures, **kwargs)


def port_forward(port: int, target: str) -> None:
    """Forward ports to local host

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
