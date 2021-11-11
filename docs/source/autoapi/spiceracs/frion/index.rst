:py:mod:`spiceracs.frion`
=========================

.. py:module:: spiceracs.frion

.. autoapi-nested-parse::

   Correct for the ionosphere in parallel

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   spiceracs.frion.cli
   spiceracs.frion.correct_worker
   spiceracs.frion.main
   spiceracs.frion.predict_worker



.. py:function:: cli()

   
   Command-line interface
















   ..
       !! processed by numpydoc !!

.. py:function:: correct_worker(beam: Dict, outdir: str, field: str, predict_file: str, island_id: str) -> pymongo.UpdateOne

   
   Apply FRion corrections to a single island

   :Parameters: * **beam** (*Dict*) -- MongoDB beam document
                * **outdir** (*str*) -- Output directory
                * **field** (*str*) -- RACS field name
                * **predict_file** (*str*) -- FRion prediction file
                * **island_id** (*str*) -- RACS island ID

   :returns: Pymongo update query
   :rtype: pymongo.UpdateOne















   ..
       !! processed by numpydoc !!

.. py:function:: main(field: str, outdir: str, host: str, client: dask.distributed.Client, username: str = None, password: str = None, database=False, verbose=True)

   
   Main script

   :Parameters: * **field** (*str*) -- RACS field name
                * **outdir** (*str*) -- Output directory
                * **host** (*str*) -- MongoDB host IP address
                * **client** (*Client*) -- Dask distributed client
                * **username** (*str, optional*) -- Mongo username. Defaults to None.
                * **password** (*str, optional*) -- Mongo passwrod. Defaults to None.
                * **database** (*bool, optional*) -- Update database. Defaults to False.
                * **verbose** (*bool, optional*) -- Verbose output. Defaults to True.















   ..
       !! processed by numpydoc !!

.. py:function:: predict_worker(island: Dict, field: str, beam: Dict, start_time: astropy.time.Time, end_time: astropy.time.Time, freq: numpy.ndarray, cutdir: str, plotdir: str) -> str

   
   Make FRion prediction for a single island

   :Parameters: * **island** (*Dict*) -- Pymongo island document
                * **field** (*str*) -- RACS field name
                * **beam** (*Dict*) -- Pymongo beam document
                * **start_time** (*Time*) -- Start time of the observation
                * **end_time** (*Time*) -- End time of the observation
                * **freq** (*np.ndarray*) -- Array of frequencies with units
                * **cutdir** (*str*) -- Cutout directory
                * **plotdir** (*str*) -- Plot directory

   :returns: Prediction file name
   :rtype: str















   ..
       !! processed by numpydoc !!

