:py:mod:`casda_prepare`
=======================

.. py:module:: casda_prepare

.. autoapi-nested-parse::

   Prepare files for CASDA upload

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   casda_prepare.cli
   casda_prepare.convert_pdf
   casda_prepare.convert_spectra
   casda_prepare.find_cubes
   casda_prepare.find_plots
   casda_prepare.find_spectra
   casda_prepare.main
   casda_prepare.make_polspec
   casda_prepare.make_thumbnail
   casda_prepare.update_cube



.. py:function:: cli()

   
   Command line interface
















   ..
       !! processed by numpydoc !!

.. py:function:: convert_pdf(pdf_file: str, plots_dir: str, spec_dir: str) -> None

   
   Convert a PDF to a PNG

   :Parameters: **pdf_file** (*str*) -- PDF file to convert















   ..
       !! processed by numpydoc !!

.. py:function:: convert_spectra(spectrum: str, polcat: astropy.table.Table, spec_dir: str = '.') -> Tuple[astropy.units.Quantity, astropy.units.Quantity, astropy.units.Quantity, numpy.ndarray]

   
   Convert a ascii spectrum to FITS

   :Parameters: * **spectrum** (*str*) -- Name of ASCII spectrum file
                * **spec_dir** (*str, optional*) -- Directory to save FITS spectrum. Defaults to '.'.















   ..
       !! processed by numpydoc !!

.. py:function:: find_cubes(data_dir: str = '.') -> list

   
   Find cubelets in a directory

   :Parameters: **data_dir** (*str, optional*) -- Data containg cutouts directory. Defaults to ".".

   :returns: List of cubelets
   :rtype: list















   ..
       !! processed by numpydoc !!

.. py:function:: find_plots(data_dir: str = '.') -> list

   
   Find plots in a directory

   :Parameters: **data_dir** (*str, optional*) -- Data containg cutouts directory. Defaults to ".".

   :returns: List of plots
   :rtype: list















   ..
       !! processed by numpydoc !!

.. py:function:: find_spectra(data_dir: str = '.') -> list

   
   Find spectra in from cutouts directory

   :Parameters: **data_dir** (*str, optional*) -- Directory containing cutouts directory. Defaults to ".".

   :returns: List of spectra in ascii format
   :rtype: list















   ..
       !! processed by numpydoc !!

.. py:function:: main(polcatf: str, client: dask.distributed.Client, data_dir: str = '.', do_update_cubes: bool = False, do_convert_spectra: bool = False, do_convert_plots: bool = False, verbose: bool = False, test: bool = False)

   
   Main function
















   ..
       !! processed by numpydoc !!

.. py:function:: make_polspec(casda_dir: str, polcat: astropy.table.Table, freqs: numpy.ndarray, data: numpy.ndarray, noises: numpy.ndarray, gauss_ids: numpy.ndarray) -> None

   
   Make a PolSpectra table

   :Parameters: * **casda_dir** (*str*) -- CASDA directory
                * **polcat** (*Table*) -- Polarisation catalogue
                * **freqs** (*np.ndarray*) -- Array of frequency arrays
                * **data** (*np.ndarray*) -- Array of data arrays
                * **noises** (*np.ndarray*) -- Array of noise arrays
                * **gauss_ids** (*np.ndarray*) -- Array of Gaussian IDs















   ..
       !! processed by numpydoc !!

.. py:function:: make_thumbnail(cube_f: str, cube_dir: str)

   
















   ..
       !! processed by numpydoc !!

.. py:function:: update_cube(cube: str, cube_dir: str) -> None

   
   Update cube headers and symlink to CASDA area

   :Parameters: * **cube** (*str*) -- Cubelet path
                * **cube_dir** (*str*) -- CASDA cublet directory















   ..
       !! processed by numpydoc !!

