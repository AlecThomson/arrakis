:py:mod:`spiceracs.rmclean_oncuts`
==================================

.. py:module:: spiceracs.rmclean_oncuts

.. autoapi-nested-parse::

   Run RM-synthesis on cutouts in parallel

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   spiceracs.rmclean_oncuts.cli
   spiceracs.rmclean_oncuts.main
   spiceracs.rmclean_oncuts.rmclean1d
   spiceracs.rmclean_oncuts.rmclean3d



.. py:function:: cli()

   
   Command-line interface
















   ..
       !! processed by numpydoc !!

.. py:function:: main(field: str, outdir: str, host: str, client: dask.distributed.Client, username: str = None, password: str = None, dimension='1d', verbose=True, database=False, savePlots=True, validate=False, limit: int = None, cutoff: float = -3, maxIter=10000, gain=0.1, window=False, showPlots=False, rm_verbose=False)

   
   Main script

   :Parameters: * **field** (*str*) -- RACS field name.
                * **outdir** (*str*) -- Output directory.
                * **host** (*str*) -- MongoDB host IP.
                * **client** (*Client*) -- Dask client.
                * **username** (*str, optional*) -- Mongo username. Defaults to None.
                * **password** (*str, optional*) -- Mongo password. Defaults to None.
                * **dimension** (*str, optional*) -- Which dimension to run RM-CLEAN. Defaults to "1d".
                * **verbose** (*bool, optional*) -- Verbose output. Defaults to True.
                * **database** (*bool, optional*) -- Update database. Defaults to False.
                * **savePlots** (*bool, optional*) -- Save plots. Defaults to True.
                * **validate** (*bool, optional*) -- Run validation. Defaults to False.
                * **limit** (*int, optional*) -- Limit number of sources processed. Defaults to None.
                * **cutoff** (*float, optional*) -- CLEAN cutoff (in sigma). Defaults to -3.
                * **maxIter** (*int, optional*) -- Max CLEAN iterations. Defaults to 10000.
                * **gain** (*float, optional*) -- Clean gain. Defaults to 0.1.
                * **showPlots** (*bool, optional*) -- Show interactive plots. Defaults to False.
                * **rm_verbose** (*bool, optional*) -- Verbose output from RM-CLEAN. Defaults to False.















   ..
       !! processed by numpydoc !!

.. py:function:: rmclean1d(comp: dict, outdir: str, cutoff: float = -3, maxIter=10000, gain=0.1, showPlots=False, savePlots=False, rm_verbose=True, window=False) -> pymongo.UpdateOne

   
   1D RM-CLEAN

   :Parameters: * **comp** (*dict*) -- Mongo entry for component.
                * **outdir** (*str*) -- Output directory.
                * **cutoff** (*float, optional*) -- CLEAN cutouff (in sigma). Defaults to -3.
                * **maxIter** (*int, optional*) -- Maximum CLEAN interation. Defaults to 10000.
                * **gain** (*float, optional*) -- CLEAN gain. Defaults to 0.1.
                * **showPlots** (*bool, optional*) -- Show CLEAN plots. Defaults to False.
                * **savePlots** (*bool, optional*) -- Save CLEAN plots. Defaults to False.
                * **rm_verbose** (*bool, optional*) -- Verbose RM-CLEAN. Defaults to True.

   :returns: MongoDB update query.
   :rtype: pymongo.UpdateOne















   ..
       !! processed by numpydoc !!

.. py:function:: rmclean3d(island: dict, outdir: str, cutoff: float = -3, maxIter=10000, gain=0.1, rm_verbose=False) -> pymongo.UpdateOne

   
   Run RM-CLEAN on 3D cube

   :Parameters: * **island** (*dict*) -- MongoDB island entry.
                * **outdir** (*str*) -- Output directory.
                * **cutoff** (*float, optional*) -- CLEAN cutoff (in sigma). Defaults to -3.
                * **maxIter** (*int, optional*) -- Max CLEAN iterations. Defaults to 10000.
                * **gain** (*float, optional*) -- CLEAN gain. Defaults to 0.1.
                * **rm_verbose** (*bool, optional*) -- Verbose output. Defaults to False.

   :returns: MongoDB update query.
   :rtype: pymongo.UpdateOne















   ..
       !! processed by numpydoc !!

