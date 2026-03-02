Simple LSM configuration
========================

The "simple LSM" example shows how to couple a DALES simulation to the
land surface model (LSM) using spatially uniform soil and skin
properties.

What this example demonstrates
------------------------------

- Adding an ``LSMModule`` on top of a basic atmosphere and grid.
- Specifying soil temperature and moisture profiles that are uniform in
  space.
- Applying a simple land-use modification over the whole domain.

Core setup
----------

The Python structure is

.. literalinclude: ../../tests/sim_builders/test_lsm.py
   :pyobject: simple_LSM_case

Use cases
---------

- Testing that LSM input files are produced correctly for a small,
  uniform domain.
- Providing a starting point for more detailed, spatially varying soil
  and land-use specifications.
  TODO