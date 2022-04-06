:py:mod:`spiceracs.process_region`
==================================

.. py:module:: spiceracs.process_region

.. autoapi-nested-parse::

   SPICE-RACS multi-field pipeline

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   spiceracs.process_region.cli
   spiceracs.process_region.main
   spiceracs.process_region.merge_task



.. py:function:: cli()

   
   Command-line interface
















   ..
       !! processed by numpydoc !!

.. py:function:: main(args: configargparse.Namespace) -> None

   
   Main script

   :Parameters: **args** (*configargparse.Namespace*) -- Command line arguments.















   ..
       !! processed by numpydoc !!

.. py:function:: merge_task(skip: bool, **kwargs) -> prefect.Task

   
   Cutout task

   Kwargs passed to merge_fields.main

   :Parameters: **skip** (*bool*) -- Whether to skip this task

   :raises signals.SKIP: If task is skipped

   :returns: Runs merge_fields.main
   :rtype: Task















   ..
       !! processed by numpydoc !!

