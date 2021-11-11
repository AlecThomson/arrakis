:py:mod:`spiceracs.linmos`
==========================

.. py:module:: spiceracs.linmos

.. autoapi-nested-parse::

   Run LINMOS on cutouts in parallel

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   spiceracs.linmos.cli
   spiceracs.linmos.gen_seps
   spiceracs.linmos.genparset
   spiceracs.linmos.get_yanda
   spiceracs.linmos.linmos
   spiceracs.linmos.main



.. py:function:: cli()

   
   Command-line interface
















   ..
       !! processed by numpydoc !!

.. py:function:: gen_seps(field: str) -> astropy.table.Table

   
   Get separation table for a given RACS field

   :Parameters: **field** (*str*) -- RACS field name.

   :returns: Table of separation for each beam.
   :rtype: Table















   ..
       !! processed by numpydoc !!

.. py:function:: genparset(field: str, src_name: str, beams: dict, stoke: str, datadir: str, septab: astropy.table.Table, holofile: str) -> str

   
   Generate parset for LINMOS

   :Parameters: * **field** (*str*) -- RACS field name.
                * **src_name** (*str*) -- RACE source name.
                * **beams** (*dict*) -- Mongo entry for RACS beams.
                * **stoke** (*str*) -- Stokes parameter.
                * **datadir** (*str*) -- Data directory.
                * **septab** (*Table*) -- Table of separations.
                * **holofile** (*str*) -- Full path to holography file.

   :raises Exception: If no files are found.

   :returns: Path to parset file.
   :rtype: str















   ..
       !! processed by numpydoc !!

.. py:function:: get_yanda(version='1.3.0') -> str

   
   Pull yandasoft image from dockerhub.

   :Parameters: **version** (*str, optional*) -- Yandasoft version. Defaults to "1.3.0".

   :returns: Path to yandasoft image.
   :rtype: str















   ..
       !! processed by numpydoc !!

.. py:function:: linmos(parset: str, fieldname: str, image: str, verbose=False) -> pymongo.UpdateOne

   
   Run linmos

   :Parameters: * **parset** (*str*) -- Path to parset file.
                * **fieldname** (*str*) -- Name of RACS field.
                * **image** (*str*) -- Name of Yandasoft image.
                * **verbose** (*bool, optional*) -- Verbose output. Defaults to False.

   :raises Exception: If LINMOS fails.
   :raises Exception: LINMOS output not found.

   :returns: Mongo update object.
   :rtype: pymongo.UpdateOne















   ..
       !! processed by numpydoc !!

.. py:function:: main(field: str, datadir: str, client: dask.distributed.Client, host: str, holofile: str, username: str = None, password: str = None, yanda='1.3.0', stokeslist: List[str] = None, verbose=True) -> None

   
   Main script

   :Parameters: * **field** (*str*) -- RACS field name.
                * **datadir** (*str*) -- Data directory.
                * **client** (*Client*) -- Dask client.
                * **host** (*str*) -- MongoDB host IP.
                * **holofile** (*str*) -- Path to primary beam file.
                * **username** (*str, optional*) -- Mongo username. Defaults to None.
                * **password** (*str, optional*) -- Mongo password. Defaults to None.
                * **yanda** (*str, optional*) -- Yandasoft version. Defaults to "1.3.0".
                * **stokeslist** (*List[str], optional*) -- Stokes parameters to process. Defaults to None.
                * **verbose** (*bool, optional*) -- Verbose output. Defaults to True.















   ..
       !! processed by numpydoc !!

