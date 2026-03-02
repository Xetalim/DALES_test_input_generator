# Time-dependent scalar parameters in profiles

This example demonstrates how to use `TimeDependentScalar` inside the
parameter dictionaries of `AtmosphericProfile` and surface modules.

## What this example demonstrates

- Providing a dense list of timesteps to `TimedependentModule`.
- Using `TimeDependentScalar` for surface temperature that evolves
  smoothly over time.

## Pattern

The time-dependent configuration is

```python
from modular_dales.modular.time_dependent import TimedependentModule
from modular_dales.modular.time_dependent_scalars import TimeDependentScalar

timesteps = np.append(np.insert(np.linspace(0, 300, num=5), 1, 0.1), 720)
sim += TimedependentModule(timesteps=timesteps)

atmo = AtmosphereModule()
# basic, time-independent background profiles
...

sim += ConstantSurfaceTemperatureModule(
    thls=TimeDependentScalar(
        times=timesteps,
        values=np.append(
            np.insert(
                np.linspace(293.15, 295.15, num=len(timesteps) - 2),
                0,
                293.15,
            ),
            295.15,
        ),
    ),
    z0mav=0.0001,
    z0hav=0.0001,
    ps=100000,
    albedoav=0.22,
)
```

This leads to a forcing file with multiple time levels, where the
surface potential temperature follows a prescribed curve.

## Use cases

- Diurnal or multi-hour experiments where the surface signal is
  prescribed rather than computed interactively.
- Tests where you need fine-grained control over the time axis of the
  forcing file.
  TODO
