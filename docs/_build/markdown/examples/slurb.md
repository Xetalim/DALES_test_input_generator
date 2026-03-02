# SLURB surface configuration

The “SLURB” example demonstrates how to add a SLURB surface module on
Top of a land-surface configuration. SLURB provides additional control
Over sub-grid surface properties, such as modified albedo fields.

## What this example demonstrates

- Attaching a `SLURBModule` to an existing `LSMModule` setup.
- Modifying surface-related fields (e.g. albedo) over a chosen
  geometry.

## Basic pattern

Starting from a small LSM-enabled case, the SLURB part looks like

```python
from modular_dales.Surface.LSM.SLuRB.slurb import SLURBModule, SLURBModification

slurb = SLURBModule(deep_soil_temperature=293)
slurb += SLURBModification(
    geometry="all",
    vars=[{"varname": "albedo_av", "value": 10}],
    params={},
)
sim += slurb
```

This is layered on top of the usual grid, atmosphere, LSM and radiation
modules, and then processed via `sim.sim_preprocessing_pipeline()`.

## When to use this pattern

- When you want a simple, domain-wide perturbation to land-surface
  parameters.
- As a starting point for more complex SLURB setups with spatially
  varying geometries or multiple variable modifications.
  TODO
