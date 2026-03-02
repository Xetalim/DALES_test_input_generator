Point-source emission example
=============================

This example shows how to add point-source emissions (e.g. stacks) for
tracer species such as CO2 or SO2.

What this example demonstrates
------------------------------

- Defining tracers via ``EmissionTracer``.
- Placing multiple ``EmissionPointSource`` objects in the horizontal
grid.
- Combining emission with a simple atmospheric and surface setup.

Core emission configuration
---------------------------

The emission-related part of the example follows

.. code-block:: python

  from modular_dales.Emission.emission import (
      EmissionModule,
      EmissionPointSource,
      EmissionTracer,
  )

  emis = EmissionModule()

  emis += EmissionTracer(
      name="co2",
      long_name="Carbon Dioxide (CO2)",
      unit="ppm",
      molar_mass=44.009,
      lemis=True,
  )
  emis += EmissionTracer(
      name="so2",
      long_name="Sulfur Dioxide (SO2)",
      unit="ppb",
      molar_mass=64.066,
      lemis=True,
  )

  emis += EmissionPointSource(
      tracer_name="co2",
      x_idx=4,
      y_idx=4,
      height=10.0,
      temperature=293.0,
      volume=1.0,
      emission=10.0,
      stack_exit_area=1.0,
  )
  emis += EmissionPointSource(
      tracer_name="so2",
      x_idx=2,
      y_idx=2,
      height=15.0,
      temperature=293.0,
      volume=1.0,
      emission=3.0,
      stack_exit_area=1.0,
  )

  sim += emis

You can attach this to a small grid, basic atmosphere and a constant
surface module, then run the preprocessing pipeline to produce tracer
emission input files.

Use cases
---------

- Idealised plume dispersion studies with a few stacks.
- Testing that emitted tracers and their units are correctly integrated
  into the model input.
  TODO