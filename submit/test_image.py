#!/usr/bin/env python3
#SBATCH --output=test_image_%j.log
#SBATCH --error=test_image_%j.log
#SBATCH --time=2-00:00:00
#SBATCH --tasks=2
#SBATCH --account=OD-217087

import logging as log

import yaml
from dask.distributed import Client
from dask_jobqueue import SLURMCluster
from IPython import embed
from spiceracs.utils import port_forward
from spiceracs import imager

def main():
    with open("/scratch2/tho822/spiceracs/spiceracs/spiceracs/configs/petrichor.yaml") as f:
        config = yaml.safe_load(f)

    config["job_cpu"] = config["cores"]
    config["cores"] = 1
    config["processes"] = 1

    cluster = SLURMCluster(
        **config,
    )
    log.debug(f"Submitted scripts will look like: \n {cluster.job_script()}")
    # exit()
    # cluster.scale(36)
    cluster.adapt(minimum=1, maximum=36)


    client = Client(cluster)

    port = client.scheduler_info()["services"]["dashboard"]
    port_forward(port, "petrichor-i1")
    log.info(client.scheduler_info()["services"])

    results = imager.main(
        msdir="/scratch2/tho822/spiceracs/RACS_1213-25A",
        out_dir="/scratch2/tho822/spiceracs/RACS_1213-25A",
        cutoff=25,
        taper=20,
        pols="IQU",
        nchan=36,
    )

    log.info(results)

if __name__ == "__main__":
    log.basicConfig(
        level=log.DEBUG,
        format="%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )
    main()