Open boundary conditions: overview
==================================

Several examples in the test suite exercise open boundary conditions.
This overview page explains the main ideas and how to approach such
configurations.

Concepts illustrated
--------------------

- Generating NetCDF files with atmospheric and surface fields that can
  be used as open boundary input.
- Inspecting those files (e.g. with ``ncview``) to verify their
  contents.

Typical workflow
----------------

A simple smoke-test style setup (see tests/openbc_test_input and
related tests) typically:

- Creates a small NetCDF file for the atmosphere and surface using
  helper functions.
- Opens those files with a visualisation tool such as ``ncview``.
- Later examples (see the dedicated open-BC examples) show how to wire
  these profiles directly into a ``dales_simulation`` using dedicated
  open-boundary helpers.

For fully integrated, AtmosphereProfiles-based open boundaries, see the
:doc:`openbc_atmosphere_profiles` example.
  TODO