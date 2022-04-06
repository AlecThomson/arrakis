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

   spiceracs.makecat.add_metadata
   spiceracs.makecat.cli
   spiceracs.makecat.cuts_and_flags
   spiceracs.makecat.get_alpha
   spiceracs.makecat.get_fit_func
   spiceracs.makecat.get_integration_time
   spiceracs.makecat.is_leakage
   spiceracs.makecat.lognorm_from_percentiles
   spiceracs.makecat.main
   spiceracs.makecat.sigma_add_fix



.. py:function:: add_metadata(vo_table)

   
   Add metadata to VO Table for CASDA

   :Parameters: **vo_table** (*vot*) -- VO Table object

   :returns: VO Table object with metadata
   :rtype: vot















   ..
       !! processed by numpydoc !!

.. py:function:: cli()

   
   Command-line interface
















   ..
       !! processed by numpydoc !!

.. py:function:: cuts_and_flags(cat)

   
   Cut out bad sources, and add flag columns

   A flag of 'True' means the source is bad.

   :Parameters: **cat** (*rmt*) -- Catalogue to cut and flag















   ..
       !! processed by numpydoc !!

.. py:function:: get_alpha(cat)

   
















   ..
       !! processed by numpydoc !!

.. py:function:: get_fit_func(tab, nbins=21, offset=0.002)

   
   Fit an envelope to define leakage sources

   :Parameters: * **tab** (*Table*) -- Catalogue to fit
                * **nbins** (*int, optional*) -- Number of bins along seperation axis. Defaults to 21.

   :returns: 3rd order polynomial fit.
   :rtype: np.polynomial.Polynomial.fit















   ..
       !! processed by numpydoc !!

.. py:function:: get_integration_time(cat, field_col)

   
















   ..
       !! processed by numpydoc !!

.. py:function:: is_leakage(frac, sep, fit)

   
   Determine if a source is leakage

   :Parameters: * **frac** (*float*) -- Polarised fraction
                * **sep** (*float*) -- Separation from tile centre
                * **fit** (*function*) -- Fitting function

   :returns: True if source is leakage
   :rtype: bool















   ..
       !! processed by numpydoc !!

.. py:function:: lognorm_from_percentiles(x1, p1, x2, p2)

   
   Return a log-normal distribuion X parametrized by:

   P(X < p1) = x1
   P(X < p2) = x2















   ..
       !! processed by numpydoc !!

.. py:function:: main(field: str, host: str, username: str = None, password: str = None, verbose=True, outfile: str = None) -> None

   
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

.. py:function:: sigma_add_fix(tab)

   
















   ..
       !! processed by numpydoc !!

