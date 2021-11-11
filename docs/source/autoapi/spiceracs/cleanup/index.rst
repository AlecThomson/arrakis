:py:mod:`spiceracs.cleanup`
===========================

.. py:module:: spiceracs.cleanup

.. autoapi-nested-parse::

   DANGER ZONE: Purge directories of un-needed FITS files.

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   spiceracs.cleanup.cleanup
   spiceracs.cleanup.cli
   spiceracs.cleanup.main



.. py:function:: cleanup(workdir: str, stoke: str) -> None

   
   Clean up beam images

   :Parameters: * **workdir** (*str*) -- Directory containing images
                * **stoke** (*str*) -- Stokes parameter















   ..
       !! processed by numpydoc !!

.. py:function:: cli()

   
   Command-line interface
















   ..
       !! processed by numpydoc !!

.. py:function:: main(datadir: str, client: dask.distributed.Client, stokeslist: List[str] = None, verbose=True) -> None

   
   Clean up beam images

   :Parameters: * **datadir** (*str*) -- Directory with sub dir 'cutouts'
                * **client** (*Client*) -- Dask Client
                * **stokeslist** (*List[str], optional*) -- List of Stokes parameters to purge. Defaults to None.
                * **verbose** (*bool, optional*) -- Verbose output. Defaults to True.















   ..
       !! processed by numpydoc !!

