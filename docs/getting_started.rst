Getting started
===============

This page gives a minimal overview of how to build and process a DALES
case using the modular_dales public API. The examples mirror the patterns
used in the tests under tests/sim_builders/.

Basic workflow
--------------

The central object is dales_simulation, which coordinates a set of modules
that each configure part of the model (grid, atmosphere, surface, LSM,
timing, etc.). A typical workflow looks like:

1. Create a dales_simulation with a case name and machine configuration.
2. Add configuration modules (e.g. default namelist, grid, atmosphere,
   surface or LSM, radiation, time settings).
3. Run the preprocessing pipeline to generate input files.

In code, this is similar to:


.. code-block:: python

   from modular_dales import dales_simulation
   from modular_dales.Configuration.defaultnamelist import DefaultNamelistModule
   from modular_dales.Configuration.run_and_time import TimeModule
   from modular_dales.Geometry.GridDales import GridDales
   from modular_dales.Atmosphere import (
       AtmosphereModule,
       AtmosphericProfile,
       InterpolatedProfile,
   )
   from modular_dales.Surface.surface import ConstantSurfaceTemperatureModule
   from modular_dales.vars import ua, va, thetal, qt, wa, tke

   # machine_conf is a dict-like configuration, see your own
   # machine_conf.yaml and tests/sim_builders/ for concrete examples.

   sim = dales_simulation("basic_case", machine_conf)

   # Start from a default namelist
   sim += DefaultNamelistModule()

   # Define the grid
   domain_info = GridDales(
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
   sim += domain_info

   # Atmospheric profiles
   atmo = AtmosphereModule()
   atmo += AtmosphericProfile(variable=ua, shape="lin", params=dict(surf_val=3, ddz=1e-3))
   atmo += AtmosphericProfile(variable=va, shape="lin", params=dict(surf_val=0, ddz=0))
   atmo += AtmosphericProfile(variable=thetal, shape="lin", params=dict(surf_val=293.15, ddz=1e-2))
   atmo += AtmosphericProfile(variable=qt, shape="lin", params=dict(surf_val=0, ddz=0))
   atmo += AtmosphericProfile(variable=wa, shape="lin", params=dict(surf_val=0, ddz=0))
   atmo += InterpolatedProfile(variable=tke, z=[0, 4000, 5000], points=[1, 1e-8, 1e-8])
   sim += atmo

   # Simple surface module
   sim += ConstantSurfaceTemperatureModule(
       thls=293.15, z0mav=0.0001, z0hav=0.0001, ps=100000, albedoav=0.22
   )

   # Time settings
   sim += TimeModule(xtime=0.0, xday=1, xyear=2025, runtime=60)

   # Run the preprocessing pipeline to generate DALES input
   sim.sim_preprocessing_pipeline()

This will create an output directory (based on the case name and
machine configuration) containing a namelist, job file, and all input
files needed by DALES for this configuration.

More advanced cases
-------------------

For more complex setups, additional modules can be added in the same
way, for example the land surface model (LSM) and radiation:


.. code-block:: python

   from modular_dales.Radiation.radiation import RadiationModule
   from modular_dales.Surface.LSM.LSM import LSMModule, LandUseModification
   from modular_dales.Surface.LSM.modular_temps_moisture import (
       UniformSkinTemperature,
       UniformSoilMoisture,
       UniformSoilTemperature,
   )

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

   sim += RadiationModule(iradiation=4)

   sim.sim_preprocessing_pipeline()

This pattern is used in the tests in
`tests/sim_builders/test_lsm.py` and
`tests/sim_builders/test_sim2.py`, which you can consult for
runnable, end-to-end examples.

YAML round-tripping
-------------------

The dales_simulation class also supports serializing the module
configuration to YAML and loading it back, which is useful for
storing and sharing cases. See the docstring of
modular_dales.modular.dales_simulation.dales_simulation and the
round-trip tests under tests/sim_builders/ for details.
