:py:mod:`spiceracs.init_database`
=================================

.. py:module:: spiceracs.init_database

.. autoapi-nested-parse::

   Create the SPICE-RACS database

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   spiceracs.init_database.beam_database
   spiceracs.init_database.cat2beams
   spiceracs.init_database.cli
   spiceracs.init_database.field_database
   spiceracs.init_database.get_beams
   spiceracs.init_database.get_catalogue
   spiceracs.init_database.main
   spiceracs.init_database.ndix_unique
   spiceracs.init_database.source2beams
   spiceracs.init_database.source_database



.. py:function:: beam_database(islandcat, host, username=None, password=None, verbose=True)

   
















   ..
       !! processed by numpydoc !!

.. py:function:: cat2beams(mastercat: astropy.table.Table, database: astropy.table.Table, max_sep=1, verbose=True) -> Tuple[numpy.ndarray, numpy.ndarray, astropy.coordinates.Angle]

   
   Find the separations between sources in the master catalogue and the RACS beams

   :Parameters: * **mastercat** (*Table*) -- Master catalogue table.
                * **database** (*Table*) -- RACS database table.
                * **max_sep** (*int, optional*) -- Maxium source separation in degrees. Defaults to 1.
                * **verbose** (*bool, optional*) -- Verbose output. Defaults to True.

   :returns: Output of astropy.coordinates.search_around_sky
   :rtype: Tuple[np.ndarray, np.ndarray, Angle]















   ..
       !! processed by numpydoc !!

.. py:function:: cli()

   
   Command-line interface
















   ..
       !! processed by numpydoc !!

.. py:function:: field_database(host, username, password, verbose=True)

   
















   ..
       !! processed by numpydoc !!

.. py:function:: get_beams(mastercat, database, verbose=True)

   
















   ..
       !! processed by numpydoc !!

.. py:function:: get_catalogue(verbose=True)

   
















   ..
       !! processed by numpydoc !!

.. py:function:: main(args, verbose=True)

   
   Main script

   :Parameters: **args {[type]} -- commandline args**















   ..
       !! processed by numpydoc !!

.. py:function:: ndix_unique(x: numpy.ndarray) -> Tuple[numpy.ndarray, List[numpy.ndarray]]

   
   Find the N-dimensional array of indices of the unique values in x
   From https://stackoverflow.com/questions/54734545/indices-of-unique-values-in-array

   :Parameters: **x** (*np.ndarray*) -- Array of values.

   :returns:     - 1D-array of sorted unique values
                 - Array of arrays. Each array contains the indices where a
                 given value in x is found
   :rtype: Tuple[np.ndarray, np.ndarray]















   ..
       !! processed by numpydoc !!

.. py:function:: source2beams(ra: float, dec: float, database: astropy.table.Table, max_sep=1) -> astropy.table.Table

   
   Find RACS beams that contain a given source position

   :Parameters: * **ra** (*float*) -- RA of source in degrees.
                * **dec** (*float*) -- DEC of source in degrees.
                * **database** (*dict*) -- RACS database table.
                * **max_sep** (*int, optional*) -- Maximum seperation of source to beam centre in degrees. Defaults to 1.

   :returns: Subset of RACS databsae table containing beams that contain the source.
   :rtype: Table















   ..
       !! processed by numpydoc !!

.. py:function:: source_database(islandcat: astropy.table.Table, compcat: astropy.table.Table, host: str, username: str = None, password: str = None, verbose=True)

   
   Insert sources into the database

   Following https://medium.com/analytics-vidhya/how-to-upload-a-pandas-dataframe-to-mongodb-ffa18c0953c1

   :Parameters: * **islandcat** (*Table*) -- Island catalogue table.
                * **compcat** (*Table*) -- Component catalogue table.
                * **host** (*str*) -- MongoDB host IP.
                * **username** (*str, optional*) -- Mongo username. Defaults to None.
                * **password** (*str, optional*) -- Mongo host. Defaults to None.
                * **verbose** (*bool, optional*) -- Verbose output. Defaults to True.















   ..
       !! processed by numpydoc !!

