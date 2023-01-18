#!/usr/bin/env python3
#SBATCH --output=/scratch2/tho822/spiceracs/RACS_1213-25A/test_image_%j.log
#SBATCH --error=/scratch2/tho822/spiceracs/RACS_1213-25A/test_image_%j.log
#SBATCH --time=3-00:00:00
#SBATCH --tasks=1
#SBATCH --account=OD-217087
#SBATCH --qos=express

import logging as log

import yaml
from dask.distributed import Client, LocalCluster
from dask_jobqueue import SLURMCluster
from IPython import embed

from spiceracs import imager
from spiceracs.utils import port_forward


def main():
    with open("/scratch2/tho822/spiceracs/spiceracs/spiceracs/configs/petrichor.yaml") as f:
        config = yaml.safe_load(f)

    config["job_cpu"] = config["cores"]
    config["cores"] = 1
    config["processes"] = 1

    cluster = SLURMCluster(
        **config,
    )
    cluster.scale(1)
    log.debug(f"Submitted scripts will look like: \n {cluster.job_script()}")
    # # exit()
    # cluster = LocalCluster(n_workers=10, threads_per_worker=1)
    # cluster.adapt(minimum=1, maximum=36)


    client = Client(cluster)

    port = client.scheduler_info()["services"]["dashboard"]
    port_forward(port, "petrichor-i1")
    log.info(client.scheduler_info()["services"])

    results = imager.main(
        msdir="/scratch2/tho822/spiceracs/RACS_1213-25A",
        out_dir="/scratch2/tho822/spiceracs/RACS_1213-25A",
        cutoff=25,
        pols="IQU",
        nchan=36,
        minuv=200,
        mgain=0.6,
        # parallel_deconvolution=3072,
        reimage=True,
    )

    # log.info(results)

if __name__ == "__main__":
    log.basicConfig(
        level=log.DEBUG,
        format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )
    main()