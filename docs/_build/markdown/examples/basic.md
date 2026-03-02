# Basic simulation

The “basic” example is a minimal, yet non-trivial, DALES case that
illustrates how to combine a default namelist, a small 3D grid, a set
of simple atmospheric profiles, and a constant surface module.

## What this example demonstrates

- Creating a `dales_simulation` with a short runtime suitable for quick
  experiments.
- Using `DefaultNamelistModule` as a starting point.
- Defining a small grid with `GridDales`.
- Building up an `AtmosphereModule` from simple linear profiles.
- Adding a constant surface temperature via
  `ConstantSurfaceTemperatureModule`.

## Typical usage pattern

A basic configuration follows the pattern below (compare to the
implementation in tests/sim_builders/test_basic.py)

```python
from modular_dales.modular import dales_simulation
from modular_dales.Configuration.defaultnamelist import DefaultNamelistModule
from modular_dales.Configuration.run_and_time import TimeModule
from modular_dales.Geometry.GridDales import GridDales
from modular_dales.Atmosphere import AtmosphereModule, AtmosphericProfile, InterpolatedProfile
from modular_dales.Surface.surface import ConstantSurfaceTemperatureModule
from modular_dales.vars import ua, va, thetal, qt, wa, tke

sim = dales_simulation("basic", machine_conf)

sim += DefaultNamelistModule()

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

atmo = AtmosphereModule()
atmo += AtmosphericProfile(variable=ua, shape="lin", params=dict(surf_val=3, ddz=1e-3))
atmo += AtmosphericProfile(variable=va, shape="lin", params=dict(surf_val=0, ddz=0))
atmo += AtmosphericProfile(variable=thetal, shape="lin", params=dict(surf_val=293.15, ddz=1e-2))
atmo += AtmosphericProfile(variable=qt, shape="lin", params=dict(surf_val=0, ddz=0))
atmo += AtmosphericProfile(variable=wa, shape="lin", params=dict(surf_val=0, ddz=0))
atmo += InterpolatedProfile(variable=tke, z=[0, 4000, 5000], points=[1, 1e-8, 1e-8])
sim += atmo

sim += ConstantSurfaceTemperatureModule(
    thls=293.15, z0mav=0.0001, z0hav=0.0001, ps=100000, albedoav=0.22
)

sim += TimeModule(xtime=0.0, xday=1, xyear=2025, runtime=60)

sim.sim_preprocessing_pipeline()
```

## When to start from this case

Use this example as a template when you want to:

- Verify your installation and machine configuration.
- Prototype new modules or small changes to the grid or atmosphere.
- Create very lightweight input for debugging downstream DALES runs.
  TODO
