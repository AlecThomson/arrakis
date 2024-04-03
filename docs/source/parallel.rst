Parallelisation
---------------
The pipeline uses `Dask <https://www.dask.org/>`_ for parallelisation and job submission, and `Prefect <https://docs.prefect.io/latest/>`_ for pipeline orchestration. Specifically, a Prefect `DaskTaskRunner <https://prefecthq.github.io/prefect-dask/>`_ is created based on a supplied configuration file.

Any Dask cluster is supported, including `dask-jobque <https://jobqueue.dask.org/en/latest/>`_ for HPC schedulers such as Slurm. This allows the pipeline to be run on virtually any system. For example, to use the `SlurmCluster` from `dask-jobqueue`, set the following in your configuration file:

.. code-block:: yaml

    cluster_class: "dask_jobqueue.SLURMCluster"

Configuration is specicfied by a file written in `YAML <https://yaml.org/>`_. These are stored in :file:`arrakis/configs/`. Add your own configuration by adding and editing a configuration, and point the pipeline to the file (see the dask-jobqueue `docs <https://jobqueue.dask.org/en/latest/configuration.html#configuration>`_). Note that cluster configuration options should be specicfied under the `cluster_kwargs` section, and adaptative scaling options under the `adapt_kwargs` section (see examples below). For further reading on Dask's adaptive scaling, see `here <https://docs.dask.org/en/latest/adaptive.html>`_.

*Arrakis* supports two configurations to be supplied to the `spice_process` pipeline (see :ref:`Running the pipeline`) via the `--dask_config` and `--imager_dask_config` arguments. The former is used by the cutout pipeline, and the latter by the imager pipeline. Imaging typically requires more memory, and more CPUs per task, whereas the cutout pipeline requires high overall number of tasks. We provide two example configurations for CSIRO `petrichor` HPC cluster.

For the imaging pipeline, an example configuration file is:

.. code-block:: yaml

    # petrichor.yaml
    # Set up for Petrichor
    cluster_class: "dask_jobqueue.SLURMCluster"
    cluster_kwargs:
        cores: 16
        processes: 1
        name: 'spice-worker'
        memory: "128GiB"
        account: 'OD-217087'
        walltime: '1-00:00:00'
        job_extra_directives: ['--qos express']
        # interface for the workers
        interface: "ib0"
        log_directory: 'spice_logs'
        job_script_prologue: [
            'module load singularity',
            'unset SINGULARITY_BINDPATH'
        ]
        local_directory: $LOCALDIR
        silence_logs: 'info'
    adapt_kwargs:
        minimum_jobs: 1
        maximum_jobs: 36
        wait_count: 20
        target_duration: "300s"
        interval: "30s"

For the cutout pipeline, an example configuration file is:

.. code-block:: yaml

    # rm_petrichor.yaml
    # Set up for Petrichor
    cluster_class: "dask_jobqueue.SLURMCluster"
    cluster_kwargs:
        cores: 4
        processes: 4
        name: 'spice-worker'
        memory: "256GiB"
        account: 'OD-217087'
        walltime: '0-01:00:00'
        job_extra_directives: ['--qos express']
        # interface for the workers
        interface: "ib0"
        log_directory: 'spice_logs'
        job_script_prologue: [
            'module load singularity',
            'unset SINGULARITY_BINDPATH',
            'export OMP_NUM_THREADS=1'
        ]
        local_directory: $LOCALDIR
        silence_logs: 'info'
    adapt_kwargs:
        minimum: 108
        maximum: 512
        wait_count: 20
        target_duration: "5s"
        interval: "10s"
