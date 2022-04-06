:py:mod:`spiceracs.merge_fields`
================================

.. py:module:: spiceracs.merge_fields

.. autoapi-nested-parse::

   Merge multiple RACS fields

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   spiceracs.merge_fields.cli
   spiceracs.merge_fields.copy_singleton
   spiceracs.merge_fields.copy_singletons
   spiceracs.merge_fields.genparset
   spiceracs.merge_fields.main
   spiceracs.merge_fields.make_short_name
   spiceracs.merge_fields.merge_multiple_field
   spiceracs.merge_fields.merge_multiple_fields



.. py:function:: cli()

   
   Command-line interface
















   ..
       !! processed by numpydoc !!

.. py:function:: copy_singleton(beam: dict, vals: dict, merge_name: str, field_dir: str, data_dir: str) -> pymongo.UpdateOne

   
















   ..
       !! processed by numpydoc !!

.. py:function:: copy_singletons(field_dict: Dict[str, str], data_dir: str, beams_col: pymongo.collection.Collection, merge_name: str) -> list

   
















   ..
       !! processed by numpydoc !!

.. py:function:: genparset(old_ims: list, stokes: str, new_dir: str) -> str

   
















   ..
       !! processed by numpydoc !!

.. py:function:: main(fields: List[str], field_dirs: List[str], merge_name: str, output_dir: str, client: dask.distributed.Client, host: str, username: str = None, password: str = None, yanda='1.3.0', verbose: bool = True) -> str

   
















   ..
       !! processed by numpydoc !!

.. py:function:: make_short_name(name: str) -> str

   
















   ..
       !! processed by numpydoc !!

.. py:function:: merge_multiple_field(beam: dict, field_dict: dict, merge_name: str, data_dir: str, image: str) -> list

   
















   ..
       !! processed by numpydoc !!

.. py:function:: merge_multiple_fields(field_dict: Dict[str, str], data_dir: str, beams_col: pymongo.collection.Collection, merge_name: str, image: str) -> list

   
















   ..
       !! processed by numpydoc !!

