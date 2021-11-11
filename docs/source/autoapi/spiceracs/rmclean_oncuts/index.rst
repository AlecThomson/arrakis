:py:mod:`spiceracs.rmclean_oncuts`
==================================

.. py:module:: spiceracs.rmclean_oncuts

.. autoapi-nested-parse::

   
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

.. py:function:: main(field, outdir, host, client, username=None, password=None, dimension='1d', verbose=True, database=False, savePlots=True, validate=False, limit=None, cutoff=-3, maxIter=10000, gain=0.1, showPlots=False, rm_verbose=False)

   
   Main script

   :Parameters: * **field** (*str*) -- RACS field
                * **outdir** (*str*) -- Work directory (contains 'cutouts' as subdir)
                * **host** (*str*) -- MongoDB host
                * **client** (*Client*) -- Dask client
                * **dimension** (*str, optional*) -- RM-CLEAN dimension. Defaults to '1d'.
                * **verbose** (*bool, optional*) -- Verbose output. Defaults to True.
                * **database** (*bool, optional*) -- Update MongoDB. Defaults to False.
                * **validate** (*bool, optional*) -- Run on Stokes I. Defaults to False.
                * **limit** (*int, optional*) -- Limit number of sources to CLEAN. Defaults to None.
                * **cutoff** (*float, optional*) -- CLEAN cutof. Defaults to -3.
                * **maxIter** (*int, optional*) -- CLEAN max iterations. Defaults to 10000.
                * **gain** (*float, optional*) -- CLEAN gain. Defaults to 0.1.
                * **showPlots** (*bool, optional*) -- Show CLEAN plots. Defaults to False.
                * **rm_verbose** (*bool, optional*) -- Verbose RM-CLEAN. Defaults to False.















   ..
       !! processed by numpydoc !!

.. py:function:: rmclean1d(comp, outdir, cutoff=-3, maxIter=10000, gain=0.1, showPlots=False, savePlots=False, rm_verbose=True)

   
   1D RM-CLEAN

   :Parameters: * **comp_id** (*str*) -- RACS component ID
                * **host** (*str*) -- MongoDB host
                * **field** (*str*) -- RACS field
                * **cutoff** (*int, optional*) -- CLEAN cutoff. Defaults to -3.
                * **maxIter** (*int, optional*) -- CLEAN max iterations. Defaults to 10000.
                * **gain** (*float, optional*) -- CLEAN gain. Defaults to 0.1.
                * **showPlots** (*bool, optional*) -- Show plots. Defaults to False.
                * **savePlots** (*bool, optional*) -- Save plots. Defaults to False.
                * **database** (*bool, optional*) -- Update MongoDB. Defaults to False.
                * **rm_verbose** (*bool, optional*) -- Verbose RM-CLEAN. Defaults to True.















   ..
       !! processed by numpydoc !!

.. py:function:: rmclean3d(island, outdir, cutoff=-3, maxIter=10000, gain=0.1, rm_verbose=False)

   
   3D RM-CLEAN

   :Parameters: * **island_id** (*str*) -- RACS Island ID
                * **host** (*str*) -- MongoDB host
                * **field** (*str*) -- RACS field
                * **cutoff** (*int, optional*) -- CLEAN cutoff. Defaults to -3.
                * **maxIter** (*int, optional*) -- CLEAN max iterations. Defaults to 10000.
                * **gain** (*float, optional*) -- CLEAN gain. Defaults to 0.1.
                * **rm_verbose** (*bool, optional*) -- Verbose RM-CLEAN. Defaults to False.















   ..
       !! processed by numpydoc !!

