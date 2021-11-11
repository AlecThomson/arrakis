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
   spiceracs.utils.fix_header
   spiceracs.utils.get_db
   spiceracs.utils.get_field_db
   spiceracs.utils.getdata
   spiceracs.utils.getfreq
   spiceracs.utils.gettable
   spiceracs.utils.head2dict
   spiceracs.utils.port_forward
   spiceracs.utils.test_db
   spiceracs.utils.tmatchn
   spiceracs.utils.tmatchtwo
   spiceracs.utils.tqdm_dask
   spiceracs.utils.try_mkdir
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

.. py:function:: coord_to_string(coord)

   
   Convert coordinate to string without astropy

   :Parameters: **coord** (*SkyCoord*) -- Coordinate

   :returns: Tuple of RA string, Dec string
   :rtype: (str,str)















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

.. py:function:: cpu_to_use(max_cpu, count)

   
   Find number of cpus to use.

   Find the right number of cpus to use when dividing up a task, such
   that there are no remainders.

   :Parameters: * **max_cpu** (*int*) -- Maximum number of cores to use for a process.
                * **count** (*float*) -- Number of tasks.

   :returns: Maximum number of cores to be used that divides into the number















   ..
       !! processed by numpydoc !!

.. py:function:: deg_to_dms(deg)

   
   Convert degree to hms without astropy.

   :Parameters: **deg** (*float*) -- Decimal degrees

   :returns: DMS, like coord.dec.dms
   :rtype: hms_tuple















   ..
       !! processed by numpydoc !!

.. py:function:: deg_to_hms(deg)

   
   Convert degree to hms without astropy.

   :Parameters: **deg** (*float*) -- Decimal degrees

   :returns: HMS, like coord.ra.hms
   :rtype: hms_tuple















   ..
       !! processed by numpydoc !!

.. py:function:: fix_header(cutout_header, original_header)

   
















   ..
       !! processed by numpydoc !!

.. py:function:: get_db(host, username=None, password=None)

   
   Get MongoDBs

   :Parameters: * **host** (*str*) -- Mongo host IP.
                * **username** (*str, optional*) -- Username. Defaults to None.
                * **password** (*str, optional*) -- Password. Defaults to None.

   :returns: beams_col, island_col, comp_col
   :rtype: Tuple(Collection)















   ..
       !! processed by numpydoc !!

.. py:function:: get_field_db(host, username=None, password=None)

   
   Get MongoDBs

   :Parameters: * **host** (*str*) -- Mongo host IP.
                * **username** (*str, optional*) -- Username. Defaults to None.
                * **password** (*str, optional*) -- Password. Defaults to None.

   :returns: beams_col, island_col, comp_col
   :rtype: Tuple(Collection)















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

.. py:function:: getfreq(cube, outdir=None, filename=None, verbose=False)

   
   Get list of frequencies from FITS data.

   Gets the frequency list from a given cube. Can optionally save
   frequency list to disk.

   :Parameters: **cube** (*str or SpectralCube*) -- File or cube to get spectral
                axis from. If a file, it will be opened using SpectralCube.

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

.. py:function:: gettable(tabledir, keyword, verbose=True)

   
   Get the spectral and source-finding data.

   :Parameters: * **tabledir** (*str*) -- Directory containing Selavy results.
                * **keyword** (*str*) -- Glob out files containing '*.keyword.*'.

   Kwargs:
       verbose (bool): Whether to print messages.

   :returns:

             Dictionary of necessary astropy tables and
                 Spectral cubes.
   :rtype: datadict (dict)















   ..
       !! processed by numpydoc !!

.. py:function:: head2dict(h)

   
   Convert FITS header to a dict.

   Writes a cutout, as stored in source_dict, to disk. The file location
   should already be specified in source_dict. This format is intended
   for parallel use with pool.map syntax.

   :Parameters: **h** -- An astropy FITS header.

   :returns: The FITS head converted to a dict.
   :rtype: data (dict)















   ..
       !! processed by numpydoc !!

.. py:function:: port_forward(port, target)

   
   Forward ports to local host

   :Parameters: * **port** (*int*) -- port to forward
                * **target** (*str*) -- Target host















   ..
       !! processed by numpydoc !!

.. py:function:: test_db(host, username=None, password=None, verbose=True)

   
















   ..
       !! processed by numpydoc !!

.. py:function:: tmatchn(nin, inN, valuesN, matcher='sky', params=10, omode='out', out='tmatch.default.xml', verbose=True)

   
   Run STILTS tmatchn
   nin = <count>       (Integer)
       The number of input tables for this task. For each of the input
       tables N there will be associated parameters ifmtN, inN and
       icmdN.

   inN = <tableN>       (StarTable)
       The location of input table #N. This may take one of the
       following forms:
           A filename.
           A URL.
           The special value "-", meaning standard input. In this case
           the input format must be given explicitly using the ifmtN
           parameter. Note that not all formats can be streamed in this
           way.
           A system command line with either a "<" character at the
           start, or a "|" character at the end ("<syscmd" or
           "syscmd|"). This executes the given pipeline and reads from
           its standard output. This will probably only work on
           unix-like systems.

   valuesN = <expr-list>       (String[])
       Defines the values from table N which are used to determine
       whether a match has occurred. These will typically be coordinate
       values such as RA and Dec and perhaps some per-row error values
       as well, though exactly what values are required is determined
       by the kind of match as determined by matcher. Depending on the
       kind of match, the number and type of the values required will
       be different. Multiple values should be separated by whitespace;
       if whitespace occurs within a single value it must be 'quoted'
       or "quoted". Elements of the expression list are commonly just
       column names, but may be algebraic expressions calculated from
       zero or more columns as explained in Section 10.

   matcher = <matcher-name>       (MatchEngine)
       Defines the nature of the matching that will be performed.
       Depending on the name supplied, this may be positional matching
       using celestial or Cartesian coordinates, exact matching on the
       value of a string column, or other things. A list and
       explanation of the available matching algorithms is given in
       Section 7.1. The value supplied for this parameter determines
       the meanings of the values required by the params, values* and
       tuning parameter(s).
       [Default: sky]

   params = <match-params>       (String[])
       Determines the parameters of this match. This is typically one
       or more tolerances such as error radii. It may contain zero or
       more values; the values that are required depend on the match
       type selected by the matcher parameter. If it contains multiple
       values, they must be separated by spaces; values which contain a
       space can be 'quoted' or "quoted".

   omode = out|meta|stats|count|cgi|discard|topcat|samp|plastic|tosql|gui
           (ProcessingMode)
       The mode in which the result table will be output. The default
       mode is out, which means that the result will be written as a
       new table to disk or elsewhere, as determined by the out and
       ofmt parameters. However, there are other possibilities, which
       correspond to uses to which a table can be put other than
       outputting it, such as displaying metadata, calculating
       statistics, or populating a table in an SQL database. For some
       values of this parameter, additional parameters (<mode-args>)
       are required to determine the exact behaviour.
       [Default: out]
   out = <out-table>       (TableConsumer)
       The location of the output table. This is usually a filename to
       write to. If it is equal to the special value "-" (the default)
       the output table will be written to standard output.
       This parameter must only be given if omode has its default value
       of "out".
       [Default: -]















   ..
       !! processed by numpydoc !!

.. py:function:: tmatchtwo(inN, valuesN, matcher='sky', params=10, omode='out', out='tmatch.default.xml', join='1or2', verbose=True)

   
   inN = <tableN>       (StarTable)
       The location of input table #N. This may take one of the
       following forms:
           A filename.
           A URL.
           The special value "-", meaning standard input. In this case
           the input format must be given explicitly using the ifmtN
           parameter. Note that not all formats can be streamed in this
           way.
           A system command line with either a "<" character at the
           start, or a "|" character at the end ("<syscmd" or
           "syscmd|"). This executes the given pipeline and reads from
           its standard output. This will probably only work on
           unix-like systems.

   valuesN = <expr-list>       (String[])
       Defines the values from table N which are used to determine
       whether a match has occurred. These will typically be coordinate
       values such as RA and Dec and perhaps some per-row error values
       as well, though exactly what values are required is determined
       by the kind of match as determined by matcher. Depending on the
       kind of match, the number and type of the values required will
       be different. Multiple values should be separated by whitespace;
       if whitespace occurs within a single value it must be 'quoted'
       or "quoted". Elements of the expression list are commonly just
       column names, but may be algebraic expressions calculated from
       zero or more columns as explained in Section 10.

   matcher = <matcher-name>       (MatchEngine)
       Defines the nature of the matching that will be performed.
       Depending on the name supplied, this may be positional matching
       using celestial or Cartesian coordinates, exact matching on the
       value of a string column, or other things. A list and
       explanation of the available matching algorithms is given in
       Section 7.1. The value supplied for this parameter determines
       the meanings of the values required by the params, values* and
       tuning parameter(s).
       [Default: sky]

   params = <match-params>       (String[])
       Determines the parameters of this match. This is typically one
       or more tolerances such as error radii. It may contain zero or
       more values; the values that are required depend on the match
       type selected by the matcher parameter. If it contains multiple
       alues, they must be separated by spaces; values which contain a
       space can be 'quoted' or "quoted".

   omode = out|meta|stats|count|cgi|discard|topcat|samp|plastic
               |tosql|gui       (ProcessingMode)
       The mode in which the result table will be output. The default
       mode is out, which means that the result will be written as a
       new table to disk or elsewhere, as determined by the out and
       ofmt parameters. However, there are other possibilities, which
       correspond to uses to which a table can be put other than
       outputting it, such as displaying metadata, calculating
       statistics, or populating a table in an SQL database. For some
       values of this parameter, additional parameters (<mode-args>)
       are required to determine the exact behaviour.
       [Default: out]

   out = <out-table>       (TableConsumer)
       The location of the output table. This is usually a filename to
       write to. If it is equal to the special value "-" (the default)
       the output table will be written to standard output.
       This parameter must only be given if omode has its default
       value of "out".
       [Default: -]

   join = 1and2|1or2|all1|all2|1not2|2not1|1xor2       (JoinType)
       Determines which rows are included in the output table. The
       matching algorithm determines which of the rows from the first
       table correspond to which rows from the second. This parameter
       determines what to do with that information. Perhaps the most
       obvious thing is to write out a table containing only rows which
       correspond to a row in both of the two input tables. However,
       you may also want to see the unmatched rows from one or both
       input tables, or rows present in one table but unmatched in the
       other, or other possibilities. The options are:
           1and2: An output row for each row represented in both input
               tables (INNER JOIN)
           1or2: An output row for each row represented in either or
               both of the input tables (FULL OUTER JOIN)
           all1: An output row for each matched or unmatched row in
               table 1 (LEFT OUTER JOIN)
           all2: An output row for each matched or unmatched row in
               table 2 (RIGHT OUTER JOIN)
           1not2: An output row only for rows which appear in the first
               table but are not matched in the second table
           2not1: An output row only for rows which appear in the
               second table but are not matched in the first table
           1xor2: An output row only for rows represented in one of the
               input tables but not the other one
       [Default: 1and2]















   ..
       !! processed by numpydoc !!

.. py:function:: tqdm_dask(futures, **kwargs)

   
















   ..
       !! processed by numpydoc !!

.. py:function:: try_mkdir(dir_path, verbose=True)

   
















   ..
       !! processed by numpydoc !!

.. py:function:: yes_or_no(question)

   
















   ..
       !! processed by numpydoc !!

