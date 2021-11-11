:py:mod:`spiceracs.cutout`
==========================

.. py:module:: spiceracs.cutout

.. autoapi-nested-parse::

   Produce cutouts from RACS cubes

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   spiceracs.cutout.cli
   spiceracs.cutout.cutout
   spiceracs.cutout.cutout_islands
   spiceracs.cutout.find_comps
   spiceracs.cutout.get_args
   spiceracs.cutout.main
   spiceracs.cutout.unpack



Attributes
~~~~~~~~~~

.. autoapisummary::

   spiceracs.cutout.auto_download


.. py:data:: auto_download
   :annotation: = False

   
















   ..
       !! processed by numpydoc !!

.. py:function:: cli() -> None

   
   Command-line interface
















   ..
       !! processed by numpydoc !!

.. py:function:: cutout(image: str, src_name: str, beam: int, ra_hi: float, ra_lo: float, dec_hi: float, dec_lo: float, outdir: str, stoke: str, field: str, pad=3, verbose=False, dryrun=False) -> pymongo.UpdateOne

   
   Perform a cutout.

   :Parameters: * **image** (*str*) -- Name of the image file
                * **src_name** (*str*) -- Name of the RACS source
                * **beam** (*int*) -- Beam number
                * **ra_hi** (*float*) -- Upper RA bound
                * **ra_lo** (*float*) -- Lower RA bound
                * **dec_hi** (*float*) -- Upper DEC bound
                * **dec_lo** (*float*) -- Lower DEC bound
                * **outdir** (*str*) -- Output directgory
                * **stoke** (*str*) -- Stokes parameter
                * **field** (*str*) -- RACS field name
                * **pad** (*int, optional*) -- Number of beamwidths to pad. Defaults to 3.
                * **verbose** (*bool, optional*) -- Verbose output. Defaults to False.
                * **dryrun** (*bool, optional*) -- Don't save FITS files. Defaults to False.

   :returns: Update query for MongoDB
   :rtype: pymongo.UpdateOne















   ..
       !! processed by numpydoc !!

.. py:function:: cutout_islands(field: str, directory: str, host: str, client: dask.distributed.Client, username: str = None, password: str = None, verbose=True, pad=3, stokeslist: List[str] = None, verbose_worker=False, dryrun=True) -> None

   
   Perform cutouts of RACS islands in parallel.

   :Parameters: * **field** (*str*) -- RACS field name.
                * **directory** (*str*) -- Directory to store cutouts.
                * **host** (*str*) -- MongoDB host.
                * **client** (*Client*) -- Dask client.
                * **username** (*str, optional*) -- Mongo username. Defaults to None.
                * **password** (*str, optional*) -- Mongo password. Defaults to None.
                * **verbose** (*bool, optional*) -- Verbose output. Defaults to True.
                * **pad** (*int, optional*) -- Number of beamwidths to pad cutouts. Defaults to 3.
                * **stokeslist** (*List[str], optional*) -- Stokes parameters to cutout. Defaults to None.
                * **verbose_worker** (*bool, optional*) -- Worker function outout. Defaults to False.
                * **dryrun** (*bool, optional*) -- Do everything except write FITS files. Defaults to True.















   ..
       !! processed by numpydoc !!

.. py:function:: find_comps(island_id: str, comp_col: pymongo.collection.Collection) -> List[Dict]

   
   Find components for a given island

   :Parameters: * **island_id** (*str*) -- RACS island ID
                * **comp_col** (*pymongo.collection.Collection*) -- Component collection

   :returns: List of mongo entries for RACS components in island
   :rtype: List[Dict]















   ..
       !! processed by numpydoc !!

.. py:function:: get_args(island: Dict, comps: List[Dict], beam: Dict, island_id: str, outdir: str, field: str, datadir: str, stokeslist: List[str], verbose=True) -> List[Dict]

   
   Get arguments for cutout function

   :Parameters: * **island** (*str*) -- Mongo entry for RACS island
                * **comps** (*List[Dict]*) -- List of mongo entries for RACS components in island
                * **beam** (*Dict*) -- Mongo entry for the RACS beam
                * **island_id** (*str*) -- RACS island ID
                * **outdir** (*str*) -- Output directory
                * **field** (*str*) -- RACS field name
                * **datadir** (*str*) -- Input directory
                * **stokeslist** (*List[str]*) -- List of Stokes parameters to process
                * **verbose** (*bool, optional*) -- Verbose output. Defaults to True.

   :raises e: Exception
   :raises Exception: Problems with coordinates

   :returns: List of cutout arguments for cutout function
   :rtype: List[Dict]















   ..
       !! processed by numpydoc !!

.. py:function:: main(args: argparse.Namespace, verbose=True) -> None

   
   Main script

   :Parameters: * **args** (*argparse.Namespace*) -- Command-line args
                * **verbose** (*bool, optional*) -- Verbose output. Defaults to True.















   ..
       !! processed by numpydoc !!

.. py:function:: unpack(list_sq: List[List[Dict]]) -> List[Dict]

   
   Unpack list of lists

   :Parameters: **list_sq** (*List[List[Dict]]*) -- List of lists of dicts

   :returns: List of dicts
   :rtype: List[Dict]















   ..
       !! processed by numpydoc !!

