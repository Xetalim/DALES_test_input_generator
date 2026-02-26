Longer time-dependent parameter example
======================================

This variant illustrates a more compact two-time-level setup where
``TimeDependentScalar`` objects are used inside ``AtmosphericProfile``
parameter dictionaries directly.

What this example demonstrates
------------------------------

- Embedding ``TimeDependentScalar`` values for both ``surf_val`` and
  ``ddz`` in atmosphere profiles.
- Combining those with a short list of timesteps in
  ``TimedependentModule``.

Sketch
------

The key elements are

.. code-block:: python

  from modular_dales.modular.time_dependent import TimedependentModule
  from modular_dales.modular.time_dependent_scalars import TimeDependentScalar

  atmo = AtmosphereModule()
  atmo += AtmosphericProfile(
      variable=wa,
      shape="lin",
      params={
          "surf_val": TimeDependentScalar(times=[0.0, 60.0], values=[0.0, 0.02]),
          "ddz": 0.0,
      },
  )
  atmo += AtmosphericProfile(
      variable=ug,
      shape="lin",
      params={
          "surf_val": TimeDependentScalar(times=[0.0, 60.0], values=[5.0, 8.0]),
          "ddz": TimeDependentScalar(times=[0.0, 60.0], values=[1e-3, 2e-3]),
      },
  )

  sim += TimedependentModule(timesteps=[0.0, 60.0])
  sim += ConstantSurfaceTemperatureModule(
      thls=TimeDependentScalar(times=[0.0, 60.0], values=[293.15, 294.15]),
      z0mav=0.0001,
      z0hav=0.0001,
      ps=100000,
      albedoav=0.22,
  )

  sim += atmo

When to use this pattern
------------------------

- Short experiments with just a few forcing times, where you want to
  keep the configuration concise.
- As a template for extending time-dependent profiles to more
  variables.
  TODO