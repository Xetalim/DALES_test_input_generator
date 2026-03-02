LSM with HARMONIE soil initialisation
=====================================

This example configures the land-surface model using soil temperature
and moisture profiles taken from an external HARMONIE NetCDF file.

Code setup
---------------

The setup uses

.. code-block:: python

  from modular_dales.Surface.LSM.LSM import LSMModule, LandUseModification
  from modular_dales.Surface.LSM.modular_temps_moisture import (
      UniformSkinTemperature,
      SoilTemperatureMoistureFromHarmonie,
  )

  harmonie_soil_levels = [
      0.01, 0.04, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0,
      1.5, 2.0, 3.0, 5.0, 8.0, 12.0,
  ]

  lsm = LSMModule(
      ps=100000,
      z0mav=0.0001,
      z0hav=0.0001,
      iinterp_t=1,
      iinterp_theta=1,
      dz_soil=[1.89, 0.72, 0.21, 0.07],
      albedoav=0.22,
  )
  lsm += UniformSkinTemperature(305)
  lsm += SoilTemperatureMoistureFromHarmonie(
      "path_to_soil/file.nc",
      harmonie_soil_valid_time="2023-08-20T00:00:00",
      harmonie_soil_height_levels=harmonie_soil_levels,
  )
  lsm += LandUseModification(geometry="all", type="grs", params={})

  sim += lsm