YAML round-trip example (sim_2)
================================

The "sim_2" example mirrors a slightly more complex configuration that
is often used to exercise YAML serialization and deserialization.

What this example demonstrates
------------------------------

- A small 3D grid with basic land-surface coupling.
- A simple radiation setup.
- A configuration that is convenient to serialize to YAML and reload.

Key building blocks
-------------------

The case uses

.. code-block:: python

  from modular_dales import dales_simulation
  from modular_dales.Configuration.defaultnamelist import DefaultNamelistModule
  from modular_dales.Configuration.run_and_time import TimeModule
  from modular_dales.Geometry.GridDales import GridDales
  from modular_dales.Atmosphere import AtmosphereModule, AtmosphericProfile, InterpolatedProfile
  from modular_dales.Radiation.radiation import RadiationModule
  from modular_dales.Surface.LSM.LSM import LSMModule, LandUseModification
  from modular_dales.Surface.LSM.modular_temps_moisture import (
      UniformSkinTemperature,
      UniformSoilMoisture,
      UniformSoilTemperature,
  )
  from modular_dales.vars import ua, va, thetal, qt, wa, tke

  sim = dales_simulation("complicated", machine_conf)

  sim += DefaultNamelistModule()
  sim += GridDales(
      itot=16,
      jtot=16,
      kmax=20,
      xsize=160.0,
      ysize=160.0,
      kmax_soil=4,
      xlat=52.25,
      xlon=5.45,
      x0=0,
      y0=0,
      alpha=1.0,
      dz0=5.0,
  )

  atmo = AtmosphereModule()
  atmo += AtmosphericProfile(variable=ua, shape="lin", params=dict(surf_val=3, ddz=1e-3))
  atmo += AtmosphericProfile(variable=va, shape="lin", params=dict(surf_val=0, ddz=0))
  atmo += AtmosphericProfile(variable=thetal, shape="lin", params=dict(surf_val=293.15, ddz=1e-2))
  atmo += AtmosphericProfile(variable=qt, shape="lin", params=dict(surf_val=0, ddz=0))
  atmo += AtmosphericProfile(variable=wa, shape="lin", params=dict(surf_val=0, ddz=0))
  atmo += InterpolatedProfile(variable=tke, z=[0, 4000, 5000], points=[1, 1e-8, 1e-8])
  sim += atmo

  sim += TimeModule(xday=0, xtime=12.0, xyear=2023, runtime=60)
  sim += RadiationModule(iradiation=4)

  lsm = LSMModule(
      ps=100000,
      z0mav=0.0001,
      z0hav=0.0001,
      iinterp_t=1,
      iinterp_theta=1,
      dz_soil=[1.89, 0.72, 0.21, 0.07],
      albedoav=0.22,
  )
  lsm += UniformSkinTemperature(293)
  lsm += UniformSoilTemperature([293, 293, 293, 293])
  lsm += UniformSoilMoisture([0.3, 0.3, 0.3, 0.3])
  lsm += LandUseModification(geometry="all", type="grs", params={})
  sim += lsm

After constructing the simulation you can::

  yaml_text = sim.save_sim_to_yaml()
  # ... store yaml_text, edit it, or load it on another machine ...

  from modular_dales.modular import dales_simulation as dsim
  sim2 = dsim.load_sim_from_yaml(yaml_text, machine_conf=machine_conf)
  sim2.sim_preprocessing_pipeline()

When to use this pattern
------------------------

- When you want a compact, but realistic, case to test YAML
  round-tripping.
- When you need an example of how land-surface and radiation settings
  appear in the serialized YAML.
  TODO