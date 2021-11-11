:py:mod:`spiceracs.processSPICE`
================================

.. py:module:: spiceracs.processSPICE

.. autoapi-nested-parse::

   SPICE-RACS pipeline script

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   spiceracs.processSPICE.cat_task
   spiceracs.processSPICE.cleanup_task
   spiceracs.processSPICE.cli
   spiceracs.processSPICE.cut_task
   spiceracs.processSPICE.frion_task
   spiceracs.processSPICE.linmos_task
   spiceracs.processSPICE.main
   spiceracs.processSPICE.rmclean_task
   spiceracs.processSPICE.rmsynth_task



.. py:function:: cat_task(skip: bool, **kwargs) -> prefect.Task

   
   Catalogue task

   Kwargs passed to makecat.main

   :Parameters: **skip** (*bool*) -- Whether to skip this task

   :raises signals.SUCCESS: If task is skipped

   :returns: Runs makecat.main
   :rtype: Task















   ..
       !! processed by numpydoc !!

.. py:function:: cleanup_task(skip: bool, **kwargs) -> prefect.Task

   
   Cleanup task

   Kwargs passed to cleanup.main

   :Parameters: **skip** (*bool*) -- Whether to skip this task

   :raises signals.SUCCESS: If task is skipped

   :returns: Runs cleanup.main
   :rtype: Task















   ..
       !! processed by numpydoc !!

.. py:function:: cli()

   
   Command-line interface
















   ..
       !! processed by numpydoc !!

.. py:function:: cut_task(skip: bool, **kwargs) -> prefect.Task

   
   Cutout task

   Kwargs passed to cutout.cutout_islands

   :Parameters: **skip** (*bool*) -- Whether to skip this task

   :raises signals.SUCCESS: If task is skipped

   :returns: Runs cutout.cutout_islands
   :rtype: Task















   ..
       !! processed by numpydoc !!

.. py:function:: frion_task(skip: bool, **kwargs) -> prefect.Task

   
   FRion task

   Kwargs passed to frion.main

   :Parameters: **skip** (*bool*) -- Whether to skip this task

   :raises signals.SUCCESS: If task is skipped

   :returns: Runs frion.main
   :rtype: Task















   ..
       !! processed by numpydoc !!

.. py:function:: linmos_task(skip: bool, **kwargs) -> prefect.Task

   
   LINOS task

   Kwargs passed to linmos.main

   :Parameters: **skip** (*bool*) -- Whether to skip this task

   :raises signals.SUCCESS: If task is skipped

   :returns: Runs linmos.main
   :rtype: Task















   ..
       !! processed by numpydoc !!

.. py:function:: main(args: configargparse.Namespace) -> None

   
   Main script

   :Parameters: **args** (*configargparse.Namespace*) -- Command line arguments.















   ..
       !! processed by numpydoc !!

.. py:function:: rmclean_task(skip: bool, **kwargs) -> prefect.Task

   
   RM-CLEAN task

   Kwargs passed to rmclean_oncuts.main

   :Parameters: **skip** (*bool*) -- Whether to skip this task

   :raises signals.SUCCESS: If task is skipped

   :returns: Runs rmclean_oncuts.main
   :rtype: Task















   ..
       !! processed by numpydoc !!

.. py:function:: rmsynth_task(skip: bool, **kwargs) -> prefect.Task

   
   RM synth task

   Kwargs passed to rmsynth_oncuts.main

   :Parameters: **skip** (*bool*) -- Whether to skip this task

   :raises signals.SUCCESS: If task is skipped

   :returns: Runs rmsynth_oncuts.main
   :rtype: Task















   ..
       !! processed by numpydoc !!

