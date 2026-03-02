SLURB with LCZ-based land use
=============================

This example combines SLURB with an LCZ (Local Climate Zone) based
surface description, allowing you to represent heterogeneous urban and
rural patterns.

What this example demonstrates
------------------------------

- Using ``FromLCZ`` within the LSM to derive land-use properties from an
  LCZ classification.
- Applying a SLURB modification on top of LCZ-based surface fields.
- Working with a slightly larger domain and a projected grid.

Key steps
---------

The surface part of the configuration uses

.. code-block:: python

  from modular_dales.Surface.LSM.LSM import LSMModule, LandUseModification, FromLCZ
  from modular_dales.Surface.LSM.SLuRB.slurb import SLURBModule, SLURBModification
  from modular_dales.Surface.LSM.modular_temps_moisture import (
      UniformSkinTemperature,
      UniformSoilMoisture,
      UniformSoilTemperature,
  )

  lsm = LSMModule(...)
  lsm += UniformSkinTemperature(291)
  lsm += UniformSoilTemperature([288, 288, 288, 288])
  lsm += UniformSoilMoisture([
      0.36867549,
      0.25300502,
      0.14997292,
      0.16459982,
  ])
  lsm += LandUseModification(geometry="all", type="grs", params={})
  lsm += FromLCZ()
  sim += lsm

  slurb = SLURBModule(deep_soil_temperature=293)
  slurb += SLURBModification(
      geometry="all",
      vars=[{"varname": "albedo_av", "value": 10}],
      params={},
  )
  sim += slurb

When to use this pattern
------------------------

- Urban or peri-urban case studies where LCZ data is available.
- Sensitivity experiments where you adjust SLURB variables on top of
  LCZ-based defaults.
  TODO