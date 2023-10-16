Parallelisation
---------------
The pipeline uses Dask for parallelisation, and optionally for job submission. Dask be run either using either `dask-jobque <https://jobqueue.dask.org/en/latest/>`_ or `dask-mpi <http://mpi.dask.org/en/latest/>`_ for parallelisation. The latter requires a working version of the `mpi4py <https://mpi4py.readthedocs.io/en/latest/>`_ package. The pipeline currently contains configurations for the CSIRO `petrichor` supercomputer and for the `galaxy` and `magnus` supercomputers at the `Pawsey Centre <https://pawsey.org.au/>`_, which all use the Slurm job manager.

.. tip ::
    Note that mpi4py needs to point to the same MPI compiler used by the MPI executable. This can be tricky to find on some systems. If in doubt, get in touch with your local sysadmin.

Configuration is specicfied by a configuration file (written in YAML). These are stored in :file:`arrakis/configs/`. Add your own configuration by adding and editing a configuration, and point the pipeline to the file (see the dask-jobqueue `docs <https://jobqueue.dask.org/en/latest/configuration.html/>`_).

.. code-block:: yaml

    # Set up for Magnus
    cores: 24
    processes: 12
    name: 'spice-worker'
    memory: "60GB"
    project: 'ja3'
    queue: 'workq'
    walltime: '6:00:00'
    job_extra: ['-M magnus']
    # interface for the workers
    interface: "ipogif0"
    log_directory: 'spice_logs'
    env_extra: [
        'export OMP_NUM_THREADS=1',
        'source /home/$(whoami)/.bashrc',
        'conda activate spice'
    ]
    python: 'srun -n 1 -c 24 python'
    extra: [
        "--lifetime", "11h",
        "--lifetime-stagger", "5m",
    ]
    death_timeout: 300
    local_directory: '/dev/shm'
