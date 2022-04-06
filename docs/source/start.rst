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

Then run the initialisation script: ::

    (spice) $ spice_init -h
    usage: init_database.py [-h] [--username USERNAME] [--password PASSWORD] [-i ISLANDCAT] 
    [-c COMPCAT] [-v] [-l] [--field] host

        
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
    --username USERNAME   Username of mongodb.
    --password PASSWORD   Password of mongodb.
    -i ISLANDCAT, --islandcat ISLANDCAT
                            Master island RACS catalogue.
    -c COMPCAT, --compcat COMPCAT
                            Master component RACS catalogue.
    -v, --verbose         Verbose output [False].
    -l, --load            Load catalogue into database [False].
    --field               Load field table into database [False].

For extra information you can refer to the API:

* :py:mod:`spiceracs.init_database`