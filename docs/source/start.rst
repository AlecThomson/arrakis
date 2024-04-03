Getting started
---------------

Pipeline database
=================

The Arrakis pipeline relies on a MongoDB database to store and pull the results of the analysis. It is a document-oriented database that stores data in JSON format. MongoDB runs as a server, either locally or on a remote server. However you run MongoDB, you'll just need to know the IP address of the server and provide it to the pipeline. You should also add a user and password to the database with full read/write access. The pipeline will use this user and password to connect to the database.

.. attention::

   Make sure not upload your MongoDB password publicly, especially if you use the pipeline configuration file to store the password.

For example, you can start mongo using (for NUMA systems like Pawsey): ::

    host=$(hostname -i) # Get the localhost IP
    database=/path/to/your/database/dir
    mkdir $database
    numactl --interleave=all mongod --dbpath=$database --bind_ip $host --fork
        --logpath mongod.log --auth

.. tip::
    It can be very convenient to run this database on a VM service like Pawsey's Nimbus cloud. You can then access the database from anywhere, and you don't need to worry about the database being deleted when you log out. This will require some network setup, such as opening the port for MongoDB (27017) on the VM. Get in touch with your local helpdesk if you need help with this.

RACS database
=============
.. attention::
    Presently, you will also need a copy of the RACS 'database'. You will need to run:
    ```
    git clone https://bitbucket.csiro.au/scm/askap_surveys/racs.git
    ```
    (This may take a while to clone...)
    Make a note of its path to pass to the `--database-path` argument of the pipeline.

    This will no longer be required when we migrate the data to a proper database.

We are currently exploring using the `VAST database <https://bitbucket.csiro.au/projects/ASKAP_SURVEYS/repos/vast/browse>`_ for processing VAST ASKAP observations.

Stokes I catalogue:
===================
Finally, the pipeline requires a Stokes I catalogue in order to construct cutouts. At present, we only support the RACS catalogue produced by `Hale et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021PASA...38...58H>`_. You can download the catalogue from the `CSIRO Data Access Portal <https://doi.org/10.25919/8zyw-5w85>`_. We are exploring ways to generalise the pipeline to accept other catalogues. If you have a catalogue you would like to use, the easiest way to get it into the pipeline is to convert it to the same format as the RACS catalogue.

Initialisation:
===============

When you have the above component in hand, run the initialisation script: ::

    (spice) $ spice_init -h
    usage: spice_init [-h] [-u USERNAME] [-p PASSWORD] [-d DATABASE_PATH]
                    [-i ISLANDCAT] [-c COMPCAT] [-v] [-l] [-f]
                    [-e EPOCHS [EPOCHS ...]]
                    host


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


    Arrakis
    Initialisation: Create MongoDB database from RACS catalogues. Before running
    make sure to start a session of mongodb e.g. $ mongod
    --dbpath=/path/to/database --bind_ip $(hostname -i)

    positional arguments:
    host                  Host of mongodb (probably $hostname -i).

    optional arguments:
    -h, --help            show this help message and exit
    -u USERNAME, --username USERNAME
                            Username of mongodb. (default: None)
    -p PASSWORD, --password PASSWORD
                            Password of mongodb. (default: None)
    -d DATABASE_PATH, --database-path DATABASE_PATH
                            Path to RACS database (i.e. 'askap_surveys/racs'
                            repo). (default: None)
    -i ISLANDCAT, --islandcat ISLANDCAT
                            Master island RACS catalogue. (default: None)
    -c COMPCAT, --compcat COMPCAT
                            Master component RACS catalogue. (default: None)
    -v, --verbose         Verbose output (default: False)
    -l, --load            Load catalogue into database. (default: False)
    -f, --field           Load field table into database. (default: False)
    -e EPOCHS [EPOCHS ...], --epochs EPOCHS [EPOCHS ...]
                            Epochs to load. (default: [0])

For extra information you can refer to the API:

* :py:mod:`arrakis.init_database`

Optional: Prefect setup:
========================
The pipeline is run using the `Prefect <https://docs.prefect.io/core/>`_ workflow management system. It is possible to set up a persistent Prefect server (called 'Orion') that will store information about pipeline runs, and even allow the configuration and triggering of pipeline execution. This is not required, but it can be useful for running the pipeline. If the pipeline is run without a server, a temporary server will be fired up and the logs etc. will be stored locally.

.. tip::
    If you already have a VM for running your MongoDB instance, it would be convenient to run Prefect on the same machine.

To set up a Prefect Server, fist install Prefect with `pip`. You will also need Postgres installed on the server to store the Prefect data. We recommend using Singularity to install and run. Below we provide two scripts for starting the Postgres and Prefect server on a remote machine:

.. tip::
    The default configuration connection settings of the Prefect and Postgres services may not be appropriate if you expect to run many workers concurrently. If the services are being overwhelmed you will see API and database time outs in the logs of these services. To avoid this it is recommended to update the `max_connection` and `shared_buffers` properties of Postgres, and set `WEB_CONCURRENCY` to enable multiple `uvicorn` workers.

.. tip::
    In each of the scripts below, you will need to set the password for the Postgres database. You can do this by setting the environment variable ``POSTGRES_PASS``. You will also need toset the hostname of the machine running the database with ``POSTGRES_ADDR``.

    Start each of these scripts in order as a background process.

.. code-block:: bash

    #!/bin/bash -l
    #start_postgres.sh

    help="This script will configure the prefect environment, and if requested start the
    postgres server necessary.

    Options:
        -s  - will attempt to start the postgres server from a singularity container
        -h  - will print this help page

    Usage:
    postgres_database.sh [-s | -h]
    "

    START_POSTGRES=1

    while getopts 'sh' arg; do
        case $arg in
        s)
            echo "Will attempt to start postgres server"
            START_POSTGRES=0
            ;;
        *)
            echo "$help"
            exit 1
            ;;
        esac
    done

    # Now set up some postgres values
    export POSTGRES_PASS='{SET YOUR PASSWORD HERE}'
    export POSTGRES_ADDR="{PUT YOUR HOSTNAME HERE}"
    export POSTGRES_USER='postgres'
    export POSTGRES_DB=orion
    export POSTGRES_SCRATCH=$(realpath $(pwd))

    if [[ $START_POSTGRES -eq 0 ]]
    then
        # Need singulaity, and to remove the badness of pawsey
        SINGULARITY_BINDPATH="$POSTGRES_SCRATCH"

        export SINGULARITY_BINDPATH


        if [[ ! -e "${POSTGRES_SCRATCH}/pgdata" ]]; then
            echo "Creating pgdata for the postgres server operation"
            mkdir pgdata
        fi

        if [[ ! -e postgres_latest.sif ]]
        then
            echo "Downloading the latest postgres docker container"
            singularity pull docker://postgres
        fi
        SINGULARITYENV_POSTGRES_PASSWORD="$POSTGRES_PASS" SINGULARITYENV_POSTGRES_DB="$POSTGRES_DB" SINGULARITYENV_PGDATA="$POSTGRES_SCRATCH/pgdata" \
            singularity run --cleanenv --bind "$POSTGRES_SCRATCH":/var postgres_latest.sif -c max_connections=1000 -c shared_buffers=1024MB
    fi

.. code-block:: bash

    #!/bin/bash -l
    #start_prefect.sh

    help="This script will configure the prefect environment, and if requested start the
    prefect server necessary.

    Options:
        -s - will atempt to start an prefect server
        -h  - will print this help page

    "

    START_ORION=1
    while getopts 'sho' arg; do
        case $arg in
        s)
        echo "Will attempt to start orion server"
        START_SERVER=0
        ;;
        *)
            echo "$help"
            exit 1
            ;;
        esac
    done

    # Now set up some postgres values
    export POSTGRES_PASS="{SET YOUR PASSWORD HERE}"
    export POSTGRES_ADDR="{PUT YOUR HOSTNAME HERE}"
    export POSTGRES_USER='postgres'
    export POSTGRES_DB=orion
    export POSTGRES_SCRATCH=$(pwd)

    export PREFECT_API_URL="http://${POSTGRES_ADDR}:4200/api"
    export PREFECT_SERVER_API_HOST="127.0.0.1"

    export PREFECT_API_DATABASE_CONNECTION_URL="postgresql+asyncpg://$POSTGRES_USER:$POSTGRES_PASS@$POSTGRES_ADDR:5432/$POSTGRES_DB"

    # This establishes a larger number of workers for prefect on the webserver (uvicorn under the hood)
    export WEB_CONCURRENCY=16
    # These can be tweaked to allow for more persistent data connections
    export PREFECT_SQLALCHEMY_POOL_SIZE=5
    export PREFECT_SQLALCHEMY_MAX_OVERFLOW=10
    PREFECT_HOME="$(pwd)/prefect"
    export PREFECT_HOME

    if [[ $START_SERVER -eq 0 ]]
    then
        prefect server start --host 0.0.0.0
    fi


Tips on adaptive scaling:
=========================

There can be strange failure modes when a prefect based workflow is being executed on a Dask task runner on a `dask_jobqueue.SLURMCluster` object with adaptive scaling enabled. Commonly, this presents as a previously completed taskrun restarting. Depending on the actual workflow, this may outright fail (e.g. if a data product that is expected has been removed), or may run perfectly fine (e.g. wsclean clobbering existing files and reimaging). Naturally, this is not behaviour that should be encouraged.

It appears as those the issue is related job stealing among a dask workers established in an adaptive scaling mode. The error mode is not entirely clear, but as workers are started (or shutdown, whether in a controlled manner or by SLURM itself restarting the node) the dask scheduler will attempt to rebalance work. For whatever reason, keys representing tasks are marked as needing to be repeated (perhaps as an attempt by dask believing it needs to recover data that was not persistent?) and are rescheduled.

The dask environment variables below are intended to try to limit these failure modes. These should be exported in the `sbatch` launch script before the python prefect / dask entry point.

.. code-block:: bash

    # See https://docs.dask.org/en/latest/configuration.html#distributed-scheduler
    # For more information on these variables
    export DASK_DISTRIBUTED__SCHEDULER__WORKER_SATURATION=0.01
    export DASK_DISTRIBUTED__SCHEDULER__WORK_STEALING=True
    export DASK_DISTRIBUTED__SCHEDULER__WORK_STEALING_INTERVAL="120s"
    export DASK_DISTRIBUTED__SCHEDULER__WORKER_TTL="3600s"
    export DASK_DISTRIBUTED__SCHEDULER__ALLOWED_FAILURES=100
    export DASK_DISTRIBUTED__WORKER__PROFILE__INTERVAL="10000ms"
    export DASK_DISTRIBUTED__WORKER__PROFILE__CYCLE="1000000ms"
    export DASK_DISTRIBUTED__COMM__TIMEOUTS__CONNECT="300s"
    export DASK_DISTRIBUTED__COMM__SOCKET_BACKLOG="16384"
    export DASK_DISTRIBUTED__COMM__TIMEOUTS__TCP="300s"
    export DASK_DISTRIBUTED__COMM__RETRY__COUNT=12


Additionally, these should be provided to the `.adapt` call that would automatically scale the dask cluster (whether in code directory or through an appropriate cluster YAML definition file).

.. code-block:: python

    {
        minimum: 2,
        maximum: 36,
        wait_count: 20,
        target_interval: "300s",
        interval: "30s"
    }
