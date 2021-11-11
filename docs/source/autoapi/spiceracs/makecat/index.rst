:py:mod:`spiceracs.makecat`
===========================

.. py:module:: spiceracs.makecat

.. autoapi-nested-parse::

   Make a SPICE-RACS catalogue

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   spiceracs.makecat.cli
   spiceracs.makecat.main



.. py:function:: cli()

   
   Command-line interface
















   ..
       !! processed by numpydoc !!

.. py:function:: main(field: str, host: str, username: str = None, password: str = None, verbose=True, outfile: str = None, cat_format: str = None) -> None

   
   Main

   :Parameters: * **field** (*str*) -- RACS field name
                * **host** (*str*) -- MongoDB host IP
                * **username** (*str, optional*) -- Mongo username. Defaults to None.
                * **password** (*str, optional*) -- Mongo password. Defaults to None.
                * **verbose** (*bool, optional*) -- Verbose output. Defaults to True.
                * **outfile** (*str, optional*) -- Output file name. Defaults to None.
                * **cat_format** (*str, optional*) -- Type of catalogue .e.g. fits. Defaults to None.















   ..
       !! processed by numpydoc !!

