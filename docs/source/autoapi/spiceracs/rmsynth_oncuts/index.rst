:py:mod:`spiceracs.rmsynth_oncuts`
==================================

.. py:module:: spiceracs.rmsynth_oncuts

.. autoapi-nested-parse::

   
   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   spiceracs.rmsynth_oncuts.cli
   spiceracs.rmsynth_oncuts.estimate_noise_annulus
   spiceracs.rmsynth_oncuts.main
   spiceracs.rmsynth_oncuts.rms_1d
   spiceracs.rmsynth_oncuts.rmsynthoncut1d
   spiceracs.rmsynth_oncuts.rmsynthoncut3d
   spiceracs.rmsynth_oncuts.rmsynthoncut_i



.. py:function:: cli()

   
   Command-line interface
















   ..
       !! processed by numpydoc !!

.. py:function:: estimate_noise_annulus(x_center, y_center, cube)

   
   Noise estimation for annulus taken around point source. Annulus has fixed
   inner radius of 10 and outer radius of 31. Function makes an annulus shaped
   mask, then for each source applies the mask at each frequency and takes the
   standard deviation.

   ​Inputs: Array of sets of pixel coordinates (y-position,x-position) for
   sources, Stokes cube (assumes 4 axes), array of flagged channels (can be an
   empty array), number of frequency channels.

   ​Output: 2D array of standard deviation values with shape (length of
   coordinate array, number of unflagged frequency channels).















   ..
       !! processed by numpydoc !!

.. py:function:: main(field, outdir, host, client, username=None, password=None, dimension='1d', verbose=True, database=False, validate=False, limit=None, savePlots=False, weightType='variance', fitRMSF=True, phiMax_radm2=None, dPhi_radm2=None, nSamples=5, polyOrd=3, noStokesI=False, showPlots=False, not_RMSF=False, rm_verbose=False, debug=False, fit_function='log', tt0=None, tt1=None, ion=False)

   
















   ..
       !! processed by numpydoc !!

.. py:function:: rms_1d(data)

   
   Compute RMS from bounding pixels
















   ..
       !! processed by numpydoc !!

.. py:function:: rmsynthoncut1d(comp, beam, outdir, freq, field, polyOrd=3, phiMax_radm2=None, dPhi_radm2=None, nSamples=5, weightType='variance', fitRMSF=True, noStokesI=False, showPlots=False, savePlots=False, debug=False, rm_verbose=False, fit_function='log', tt0=None, tt1=None, ion=False)

   
   1D RM synthesis

   :Parameters: * **comp_id** (*str*) -- RACS component ID
                * **outdir** (*str*) -- Output directory
                * **freq** (*list*) -- Frequencies in Hz
                * **host** (*str*) -- MongoDB host
                * **field** (*str*) -- RACS field
                * **database** (*bool, optional*) -- Update MongoDB. Defaults to False.
                * **polyOrd** (*int, optional*) -- Order of fit to I. Defaults to 3.
                * **phiMax_radm2** (*float, optional*) -- Max FD. Defaults to None.
                * **dPhi_radm2** (*float, optional*) -- Delta FD. Defaults to None.
                * **nSamples** (*int, optional*) -- Samples across RMSF. Defaults to 5.
                * **weightType** (*str, optional*) -- Weight type. Defaults to 'variance'.
                * **fitRMSF** (*bool, optional*) -- Fit RMSF. Defaults to False.
                * **noStokesI** (*bool, optional*) -- Ignore Stokes I. Defaults to False.
                * **showPlots** (*bool, optional*) -- Show plots. Defaults to False.
                * **savePlots** (*bool, optional*) -- Save plots. Defaults to False.
                * **debug** (*bool, optional*) -- Turn on debug plots. Defaults to False.
                * **rm_verbose** (*bool, optional*) -- Verbose RMsynth. Defaults to False.















   ..
       !! processed by numpydoc !!

.. py:function:: rmsynthoncut3d(island_id, beam, outdir, freq, field, phiMax_radm2=None, dPhi_radm2=None, nSamples=5, weightType='variance', fitRMSF=True, not_RMSF=False, rm_verbose=False, ion=False)

   
   3D RM-synthesis

   :Parameters: * **island_id** (*str*) -- RACS Island ID
                * **freq** (*list*) -- Frequencies in Hz
                * **host** (*str*) -- Host of MongoDB
                * **field** (*str*) -- RACS field ID
                * **database** (*bool, optional*) -- Update MongoDB. Defaults to False.
                * **phiMax_radm2** (*float, optional*) -- Max Faraday depth. Defaults to None.
                * **dPhi_radm2** (*float, optional*) -- Faraday dpeth channel width. Defaults to None.
                * **nSamples** (*int, optional*) -- Samples acorss RMSF. Defaults to 5.
                * **weightType** (*str, optional*) -- Weighting type. Defaults to 'variance'.
                * **fitRMSF** (*bool, optional*) -- Fit RMSF. Defaults to False.
                * **not_RMSF** (*bool, optional*) -- Skip calculation of RMSF. Defaults to False.
                * **rm_verbose** (*bool, optional*) -- Verbose RMsynth. Defaults to False.















   ..
       !! processed by numpydoc !!

.. py:function:: rmsynthoncut_i(comp_id, outdir, freq, host, field, username=None, password=None, nSamples=5, phiMax_radm2=None, verbose=False, rm_verbose=False)

   
   RMsynth on Stokes I

   :Parameters: * **comp_id** (*str*) -- RACS component ID
                * **freq** (*list*) -- Frequencies in Hz
                * **host** (*str*) -- MongoDB host
                * **field** (*str*) -- RACS field
                * **nSamples** (*[type]*) -- Samples across the RMSF
                * **phiMax_radm2** (*float*) -- Max FD
                * **verbose** (*bool, optional*) -- Verbose output Defaults to False.
                * **rm_verbose** (*bool, optional*) -- Verbose RMsynth. Defaults to False.















   ..
       !! processed by numpydoc !!

