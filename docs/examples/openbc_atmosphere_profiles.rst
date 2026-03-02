Open boundary conditions from AtmosphereProfiles
===============================================

This example shows how to configure open boundary conditions driven by
AtmosphereProfiles and an external ``AtmosphereModule``.

What this example demonstrates
------------------------------

- Constructing a small simulation domain that uses open boundaries.
- Creating an external atmosphere description specifically for the
  boundaries.
- Nesting that description into an open-boundary helper.

Core structure
--------------

The high-level pattern is:


.. code-block:: python

    from modular_dales import (
        dales_simulation,
        do_openboundary,
        AtmosphereModule,
        AtmosphericProfile,
        Nest_in_AtmosphereProfiles,
    )
    from modular_dales.Atmosphere import TimedAtmosphereProfile
    from modular_dales.Configuration.defaultnamelist import DefaultNamelistModule
    from modular_dales.Configuration.run_and_time import TimeModule
    from modular_dales.Configuration.output_modules import EasyOutputModule
    from modular_dales.Geometry.GridDales import GridDales
    from modular_dales.Surface.surface import ConstantSurfaceTemperatureModule
    from modular_dales.vars import ua, va, w, thetal, qt, tke, ug, vg

    sim = dales_simulation("openbc_case", machine_conf)
    sim += DefaultNamelistModule()
    grid = GridDales(...)
    sim += grid

    # Internal atmosphere used inside the domain
    atmo = AtmosphereModule()
    atmo += AtmosphericProfile(variable=ua, shape="lin", params=dict(surf_val=3.0, ddz=1e-3))
    ...
    sim += atmo

    # External atmosphere used at the open boundaries
    atmo_external = AtmosphereModule()
    atmo_external += AtmosphericProfile(variable=ug, shape="lin", params=dict(surf_val=3.0, ddz=1e-3))
    atmo_external += AtmosphericProfile(variable=vg, shape="lin", params=dict(surf_val=0.0, ddz=0))
    ...

    openbc = do_openboundary(
        time0="2025-01-01T00:00:00",
        start="2025-01-01T00:00:00",
        end="2025-01-01T00:10:00",
        e12=0.1,
        dxint=grid.xsize / grid.itot,
        dyint=grid.ysize / grid.jtot,
    )

    openbc += Nest_in_AtmosphereProfiles(
        atmosphere_module=atmo_external,
        noise_boundaries=["south", "west"],
        noise_std=0.1,
        noise_seed=42,
        noise_variables=["ug", "vg"],
    )
    sim += openbc

    sim += ConstantSurfaceTemperatureModule(...)
    sim += TimeModule(xtime=0.0, xday=1, xyear=2025, runtime=3600)
    sim += EasyOutputModule(output_interval=10)

    sim.sim_preprocessing_pipeline()

Time-dependent variant
----------------------

A closely related configuration uses ``TimedAtmosphereProfile`` entries
for the external geostrophic wind profiles, so that the open boundaries
contain multiple time slices (e.g. at t=0 and t=3600 seconds). The rest
of the pattern is the same.

Use cases
---------

- Idealised inflow/outflow experiments with controlled large-scale
  forcing.
- Sensitivity studies on how boundary noise or time dependence affect
  the interior solution.
  TODO