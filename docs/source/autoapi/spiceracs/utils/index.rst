:py:mod:`spiceracs.utils`
=========================

.. py:module:: spiceracs.utils

.. autoapi-nested-parse::

   
   ..
       !! processed by numpydoc !!


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   spiceracs.utils.MyEncoder
   spiceracs.utils.TqdmProgressBar



Functions
~~~~~~~~~

.. autoapisummary::

   spiceracs.utils._samefile
   spiceracs.utils.coord_to_string
   spiceracs.utils.copyfile
   spiceracs.utils.copyfileobj
   spiceracs.utils.cpu_to_use
   spiceracs.utils.deg_to_dms
   spiceracs.utils.deg_to_hms
   spiceracs.utils.delayed_to_da
   spiceracs.utils.fix_header
   spiceracs.utils.get_db
   spiceracs.utils.get_field_db
   spiceracs.utils.getdata
   spiceracs.utils.getfreq
   spiceracs.utils.gettable
   spiceracs.utils.head2dict
   spiceracs.utils.port_forward
   spiceracs.utils.test_db
   spiceracs.utils.tqdm_dask
   spiceracs.utils.try_mkdir
   spiceracs.utils.try_symlink
   spiceracs.utils.yes_or_no



Attributes
~~~~~~~~~~

.. autoapisummary::

   spiceracs.utils.print


.. py:data:: print
   

   
















   ..
       !! processed by numpydoc !!

.. py:exception:: Error

   Bases: :py:obj:`OSError`

   
   Base class for I/O related errors.
















   ..
       !! processed by numpydoc !!

.. py:exception:: ExecError

   Bases: :py:obj:`OSError`

   
   Raised when a command could not be executed
















   ..
       !! processed by numpydoc !!

.. py:exception:: ReadError

   Bases: :py:obj:`OSError`

   
   Raised when an archive cannot be read
















   ..
       !! processed by numpydoc !!

.. py:exception:: RegistryError

   Bases: :py:obj:`Exception`

   
   Raised when a registry operation with the archiving
   and unpacking registeries fails
















   ..
       !! processed by numpydoc !!

.. py:exception:: SameFileError

   Bases: :py:obj:`Error`

   
   Raised when source and destination are the same file.
















   ..
       !! processed by numpydoc !!

.. py:exception:: SpecialFileError

   Bases: :py:obj:`OSError`

   
   Raised when trying to do a kind of operation (e.g. copying) which is
   not supported on a special file (e.g. a named pipe)
















   ..
       !! processed by numpydoc !!

.. py:class:: MyEncoder(*, skipkeys=False, ensure_ascii=True, check_circular=True, allow_nan=True, sort_keys=False, indent=None, separators=None, default=None)

   Bases: :py:obj:`json.JSONEncoder`

   
   Cutom JSON encorder.

   Parses the data stored in source_dict to JSON without
   errors.















   ..
       !! processed by numpydoc !!
   .. py:method:: default(self, obj)

      
      Implement this method in a subclass such that it returns
      a serializable object for ``o``, or calls the base implementation
      (to raise a ``TypeError``).

      For example, to support arbitrary iterators, you could
      implement default like this::

          def default(self, o):
              try:
                  iterable = iter(o)
              except TypeError:
                  pass
              else:
                  return list(iterable)
              # Let the base class default method raise the TypeError
              return JSONEncoder.default(self, o)















      ..
          !! processed by numpydoc !!


.. py:class:: TqdmProgressBar(keys, scheduler=None, interval='100ms', loop=None, complete=True, start=True, **tqdm_kwargs)

   Bases: :py:obj:`distributed.diagnostics.progressbar.ProgressBar`

   
   Tqdm for Dask
















   ..
       !! processed by numpydoc !!
   .. py:method:: _draw_bar(self, remaining, all, **kwargs)

      
















      ..
          !! processed by numpydoc !!

   .. py:method:: _draw_stop(self, **kwargs)

      
















      ..
          !! processed by numpydoc !!


.. py:function:: _samefile(src, dst)

   
















   ..
       !! processed by numpydoc !!

.. py:function:: coord_to_string(coord: astropy.coordinates.SkyCoord) -> Tuple[str, str]

   
   Convert coordinate to string without astropy

   :Parameters: **coord** (*SkyCoord*) -- Coordinate

   :returns: Tuple of RA string, Dec string
   :rtype: Tuple[str,str]















   ..
       !! processed by numpydoc !!

.. py:function:: copyfile(src, dst, *, follow_symlinks=True, verbose=True)

   
   Copy data from src to dst.

   If follow_symlinks is not set and src is a symbolic link, a new
   symlink will be created instead of copying the file it points to.















   ..
       !! processed by numpydoc !!

.. py:function:: copyfileobj(fsrc, fdst, length=16 * 1024, verbose=True)

   
















   ..
       !! processed by numpydoc !!

.. py:function:: cpu_to_use(max_cpu: int, count: int) -> int

   
   Find number of cpus to use.

   Find the right number of cpus to use when dividing up a task, such
   that there are no remainders.

   :Parameters: * **max_cpu** (*int*) -- Maximum number of cores to use for a process.
                * **count** (*int*) -- Number of tasks.

   :returns: Maximum number of cores to be used that divides into the number















   ..
       !! processed by numpydoc !!

.. py:function:: deg_to_dms(deg: float) -> astropy.coordinates.angles.dms_tuple

   
   Convert degree to hms without astropy.

   :Parameters: **deg** (*float*) -- Decimal degrees

   :returns: DMS, like coord.dec.dms
   :rtype: hms_tuple















   ..
       !! processed by numpydoc !!

.. py:function:: deg_to_hms(deg: float) -> astropy.coordinates.angles.hms_tuple

   
   Convert degree to hms without astropy.

   :Parameters: **deg** (*float*) -- Decimal degrees

   :returns: HMS, like coord.ra.hms
   :rtype: hms_tuple















   ..
       !! processed by numpydoc !!

.. py:function:: delayed_to_da(list_of_delayed: List[dask.delayed], chunk: int = None) -> dask.array.Array

   
   Convert list of delayed arrays to a dask array

   :Parameters: * **list_of_delayed** (*List[delayed]*) -- List of delayed objects
                * **chunk** (*int, optional*) -- Chunksize to use. Defaults to None.

   :returns: Dask array
   :rtype: da.Array















   ..
       !! processed by numpydoc !!

.. py:function:: fix_header(cutout_header: astropy.io.fits.Header, original_header: astropy.io.fits.Header) -> astropy.io.fits.Header

   
   Make cutout header the same as original header

   :Parameters: * **cutout_header** (*fits.Header*) -- Cutout header
                * **original_header** (*fits.Header*) -- Original header

   :returns: Fixed header
   :rtype: fits.Header















   ..
       !! processed by numpydoc !!

.. py:function:: get_db(host: str, username: str = None, password: str = None) -> Tuple[pymongo.collection.Collection, pymongo.collection.Collection, pymongo.collection.Collection]

   
   Get MongoDBs

   :Parameters: * **host** (*str*) -- Mongo host IP.
                * **username** (*str, optional*) -- Username. Defaults to None.
                * **password** (*str, optional*) -- Password. Defaults to None.

   :returns: beams_col, island_col, comp_col
   :rtype: Tuple[pymongo.Collection, pymongo.Collection, pymongo.Collection]















   ..
       !! processed by numpydoc !!

.. py:function:: get_field_db(host: str, username=None, password=None) -> pymongo.collection.Collection

   
   Get MongoDBs

   :Parameters: * **host** (*str*) -- Mongo host IP.
                * **username** (*str, optional*) -- Username. Defaults to None.
                * **password** (*str, optional*) -- Password. Defaults to None.

   :returns: beams_col, island_col, comp_col
   :rtype: pymongo.Collection















   ..
       !! processed by numpydoc !!

.. py:function:: getdata(cubedir='./', tabledir='./', mapdata=None, verbose=True)

   
   Get the spectral and source-finding data.

   :Parameters: * **cubedir** -- Directory containing data cubes in FITS format.
                * **tabledir** -- Directory containing Selavy results.
                * **mapdata** -- 2D FITS image which corresponds to Selavy table.

   Kwargs:
       verbose (bool): Whether to print messages.

   :returns:

             Dictionary of necessary astropy tables and
                 Spectral cubes.
   :rtype: datadict (dict)















   ..
       !! processed by numpydoc !!

.. py:function:: getfreq(cube: str, outdir: str = None, filename: str = None)

   
   Get list of frequencies from FITS data.

   Gets the frequency list from a given cube. Can optionally save
   frequency list to disk.

   :Parameters: **cube** (*str*) -- File to get spectral axis from.

   Kwargs:
       outdir (str): Where to save the output file. If not given, data
           will not be saved to disk.

       filename (str): Name of frequency list file. Requires 'outdir'
           to also be specified.

       verbose (bool): Whether to print messages.

   :returns: Frequencies of each channel in the input cube.
   :rtype: freq (list)















   ..
       !! processed by numpydoc !!

.. py:function:: gettable(tabledir: str, keyword: str, verbose=True) -> Tuple[astropy.table.Table, str]

   
   Get a table from a directory given a keyword to glob.

   :Parameters: * **tabledir** (*str*) -- Directory.
                * **keyword** (*str*) -- Keyword to glob for.
                * **verbose** (*bool, optional*) -- Verbose output. Defaults to True.

   :returns: Table and it's file location.
   :rtype: Tuple[Table, str]















   ..
       !! processed by numpydoc !!

.. py:function:: head2dict(h: astropy.io.fits.Header) -> Dict[str, Any]

   
   Convert FITS header to a dict.

   Writes a cutout, as stored in source_dict, to disk. The file location
   should already be specified in source_dict. This format is intended
   for parallel use with pool.map syntax.

   :Parameters: **h** -- An astropy FITS header.

   :returns: The FITS head converted to a dict.
   :rtype: data (dict)















   ..
       !! processed by numpydoc !!

.. py:function:: port_forward(port: int, target: str) -> None

   
   Forward ports to local host

   :Parameters: * **port** (*int*) -- port to forward
                * **target** (*str*) -- Target host















   ..
       !! processed by numpydoc !!

.. py:function:: test_db(host: str, username: str = None, password: str = None, verbose=True) -> None

   
   Test connection to MongoDB

   :Parameters: * **host** (*str*) -- Mongo host IP.
                * **username** (*str, optional*) -- Mongo username. Defaults to None.
                * **password** (*str, optional*) -- Mongo password. Defaults to None.
                * **verbose** (*bool, optional*) -- Verbose output. Defaults to True.

   :raises Exception: If connection fails.















   ..
       !! processed by numpydoc !!

.. py:function:: tqdm_dask(futures: dask.distributed.Future, **kwargs) -> None

   
   Tqdm for Dask futures
















   ..
       !! processed by numpydoc !!

.. py:function:: try_mkdir(dir_path: str, verbose=True)

   
   Create directory if it doesn't exist

   :Parameters: * **dir_path** (*str*) -- Path to directory
                * **verbose** (*bool, optional*) -- Verbose output. Defaults to True.















   ..
       !! processed by numpydoc !!

.. py:function:: try_symlink(src: str, dst: str, verbose=True)

   
   Create symlink if it doesn't exist

   :Parameters: * **src** (*str*) -- Source path
                * **dst** (*str*) -- Destination path
                * **verbose** (*bool, optional*) -- Verbose output. Defaults to True.















   ..
       !! processed by numpydoc !!

.. py:function:: yes_or_no(question: str) -> bool

   
   Ask a yes or no question via input()

   :Parameters: **question** (*str*) -- Question to ask

   :returns: True for yes, False for no
   :rtype: bool















   ..
       !! processed by numpydoc !!

