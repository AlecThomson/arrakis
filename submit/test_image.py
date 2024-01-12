#!/usr/bin/env python3
# SBATCH --output=/scratch2/tho822/spiceracs/pipe_test/test_image_%j.log
# SBATCH --error=/scratch2/tho822/spiceracs/pipe_test/test_image_%j.log
# SBATCH --time=1-00:00:00
# SBATCH --tasks=1
# SBATCH --cpus-per-task=1
# SBATCH --account=OD-217087
# SBATCH --qos=express

import logging
from pathlib import Path

import yaml
from astropy import units as u
from dask.distributed import Client
from dask_jobqueue import SLURMCluster
from spiceracs import imager
from spiceracs.logger import logger
from spiceracs.utils import port_forward

logger.setLevel(logging.INFO)


def main():
    with open(
        "/scratch2/projects/spiceracs/spiceracs/spiceracs/configs/petrichor.yaml"
    ) as f:
        config = yaml.safe_load(f)

    config["job_cpu"] = config["cores"]
    config["cores"] = 1
    config["processes"] = 1

    cluster = SLURMCluster(
        **config,
    )
    cluster.scale(72)
    logger.debug(f"Submitted scripts will look like: \n {cluster.job_script()}")

    client = Client(cluster)

    port = client.scheduler_info()["services"]["dashboard"]
    port_forward(port, "petrichor-i1")
    logger.info(client.scheduler_info()["services"])

    results = imager.main(
        msdir=Path("/scratch2/tho822/spiceracs/pipe_test"),
        out_dir=Path("/scratch2/tho822/spiceracs/pipe_test"),
        mgain=0.8,
        force_mask_rounds=8,
        nmiter=25,
        niter=50000000,
        local_rms=True,
        auto_mask=3.75,
        local_rms_window=60,
        auto_threshold=1,
        size=6144,
        scale=2.5 * u.arcsec,
        robust=-0.5,
        pols="IQU",
        gridder="wgridder",
        minuv=200,
        wsclean_path=Path("/scratch2/tho822/singularity_images/wsclean_force_mask.sif"),
        reimage=True,
        multiscale=False,
        # parallel_deconvolution=6144,
        absmem=float(config["memory"].replace("GB", "").replace("GiB", "")),
    )
    logs = client.get_worker_logs()


if __name__ == "__main__":
    main()
