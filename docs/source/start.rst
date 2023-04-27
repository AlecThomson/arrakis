Getting started
---------------
The SPICE-RACS pipeline relies on a MongoDB database to store and pull the results of the analysis. It is a document-oriented database that stores data in JSON format. MongoDB runs as a server, either locally or on a remote server. However you run MongoDB, you'll just need to know the IP address of the server and provide it to the pipeline. You should also add a user and password to the database with full read/write access. The pipeline will use this user and password to connect to the database.

.. attention::

   Make sure not upload your MongoDB password publicly, especially if you use the pipeline configuration file to store the password.

For example, you can start mongo using (for NUMA systems like Pawsey): ::

    host=$(hostname -i) # Get the localhost IP
    database=/path/to/your/database/dir
    mkdir $database
    numactl --interleave=all mongod --dbpath=$database --bind_ip $host --fork
        --logpath mongod.log --auth


.. attention::
    Presently, you will also need a copy of the RACS 'database'. You will need to run:
    ```
    git clone https://bitbucket.csiro.au/scm/askap_surveys/racs.git
    ```
    (This may take a while to clone...)
    Make a note of its path to pass to the `--database-path` argument of the pipeline.

    This will no longer be required when we migrate the data to a proper database.

Then run the initialisation script: ::

    (spice) $ spice_init -h
    usage: spice_init [-h] [-u USERNAME] [-p PASSWORD] [-d--database-path D__DATABASE_PATH] [-i ISLANDCAT] [-c COMPCAT] [-v] [-l] [-f] [-e EPOCH] host


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


        SPICE-RACS Initialisation:

        Create MongoDB database from RACS catalogues.

        Before running make sure to start a session of mongodb e.g.
            $ mongod --dbpath=/path/to/database --bind_ip $(hostname -i)



    positional arguments:
    host                  Host of mongodb (probably $hostname -i).

    optional arguments:
    -h, --help            show this help message and exit
    -u USERNAME, --username USERNAME
                            Username of mongodb.
    -p PASSWORD, --password PASSWORD
                            Password of mongodb.
    -d--database-path D__DATABASE_PATH
                            Path to RACS database (i.e. 'askap_surveys' repo).
    -i ISLANDCAT, --islandcat ISLANDCAT
                            Master island RACS catalogue.
    -c COMPCAT, --compcat COMPCAT
                            Master component RACS catalogue.
    -v, --verbose         Verbose output [False].
    -l, --load            Load catalogue into database [False].
    -f, --field           Load field table into database [False].
    -e EPOCH, --epoch EPOCH
                            RACS epoch to load [0].
For extra information you can refer to the API:

* :py:mod:`spiceracs.init_database`