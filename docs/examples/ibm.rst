Immersed boundary method (IBM) example
======================================

The IBM example shows how to represent obstacles (e.g. buildings)
within the domain using the immersed boundary method.

What this example demonstrates
------------------------------

- Adding an ``IBMModule`` to a standard atmospheric and LSM setup.
- Using ``IBMModification`` objects to describe simple geometries such
  as a rectangular building.

Configuration sketch
--------------------

The IBM-specific part of the configuration looks like

.. code-block:: python
  from modular_dales.IBM.IBM import IBMModule, IBMModification

  ibm = IBMModule()
  # Background (no obstacle) region
  ibm += IBMModification("all", height=0, params={})

  # A simple rectangular obstacle in index space
  ibm += IBMModification(
      "rectangle_idx",
      height=20,
      params={"minx": 4, "maxx": 6, "miny": 4, "maxy": 4},
  )

  sim += ibm

IBM is typically combined with a grid, atmosphere, LSM and radiation
setup similar to the basic or LSM cases. Once everything is attached,
you run ``sim.sim_preprocessing_pipeline()`` to generate the IBM input
files.

When to use this pattern
------------------------

- Idealised building-resolving experiments.
- Testing how flow responds to simple obstacles before moving to more
  complex geometry import workflows.
  TODO