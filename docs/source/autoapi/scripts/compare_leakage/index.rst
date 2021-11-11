:py:mod:`compare_leakage`
=========================

.. py:module:: compare_leakage

.. autoapi-nested-parse::

   The interpolation works as follows:
   Take pixels offsets x,y from reference pixel in input image, multiply by
   axis increments to get offx and offy.

   Then compute offset = arcsin(offx^2+offy^2) and angle=atan2(offx,offy),
   which should be the angular offset on the sky of the pixel position.

   For the leakage image the inverse is used.
   Take the offset and angle and turn them into pixel positions on the leakage map:

   x = sin(offset)*cos(angle)/incx + refx
   y = sin(offset)*sin(angle)/incy + refy

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   compare_leakage.cli
   compare_leakage.interpolate
   compare_leakage.main
   compare_leakage.make_plot



.. py:function:: cli()

   
   Command-line interface
















   ..
       !! processed by numpydoc !!

.. py:function:: interpolate(field, comp, beams, cutdir, septab, holofile, verbose=True)

   
















   ..
       !! processed by numpydoc !!

.. py:function:: main(field, datadir, client, host, holofile, username=None, password=None, verbose=True, snr_cut=None)

   
















   ..
       !! processed by numpydoc !!

.. py:function:: make_plot(data, comp, imfile)

   
















   ..
       !! processed by numpydoc !!

