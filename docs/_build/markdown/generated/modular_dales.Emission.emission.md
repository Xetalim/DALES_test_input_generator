# modular_dales.Emission.emission

Emission module for tracer and emission handling.

This module mirrors the IBM module pattern:

- Emission tracers and point sources are represented as dataclasses that can be
  : added to the module using the `+=` operator.
- The module then constructs the internal emissions data structures and writes
  : the required NetCDF files during `prepare_calculation` / `write_files`.

### Classes

| [`EmissionModule`](#modular_dales.Emission.emission.EmissionModule)(sim, tracers, point_sources, ...)   | Emission simulation module for tracer and emission handling.   |
|---------------------------------------------------------------------------------------------------------|----------------------------------------------------------------|
| [`EmissionPointSource`](#modular_dales.Emission.emission.EmissionPointSource)(tracer_name, x_idx, ...)  | Single emission point source tied to a tracer.                 |
| [`EmissionTracer`](#modular_dales.Emission.emission.EmissionTracer)(name, long_name, unit, molar_mass)  | Single tracer configuration.                                   |

### *class* modular_dales.Emission.emission.EmissionModule(sim: simulation_module | [None](https://docs.python.org/3/library/constants.html#None) = None, tracers: [List](https://docs.python.org/3/library/typing.html#typing.List)[[EmissionTracer](#modular_dales.Emission.emission.EmissionTracer)] = <factory>, point_sources: [List](https://docs.python.org/3/library/typing.html#typing.List)[[EmissionPointSource](#modular_dales.Emission.emission.EmissionPointSource)] = <factory>, l_emission: [bool](https://docs.python.org/3/library/functions.html#bool) = True, l_points: [bool](https://docs.python.org/3/library/functions.html#bool) = False, explicit_plume_rise: [bool](https://docs.python.org/3/library/functions.html#bool) = False, emisnames: [List](https://docs.python.org/3/library/typing.html#typing.List)[[str](https://docs.python.org/3/library/stdtypes.html#str)] = <factory>, nemis: [int](https://docs.python.org/3/library/functions.html#int) = 0)

Bases: `simulation_module`

Emission simulation module for tracer and emission handling.

Tracers and point sources can be added directly to this module, e.g.:

```default
emis = EmissionModule()
emis += EmissionTracer(
    name="co2",
    long_name="Carbon Dioxide (CO2)",
    unit="ppm",
    molar_mass=44.009,
    lemis=True,
)
emis += EmissionPointSource(
    tracer_name="co2",
    x_idx=8,
    y_idx=8,
    height=10,
    temperature=293.0,
    volume=1.0,
    emission=10.0,
    stack_exit_area=1.0,
)
```

* **Parameters:**
  **sim** – Parent simulation instance

## Namelist parameters

*Section* `NAMEMISSION`:
: - `emisnames` (field `emisnames`) (required)
  - `explicit_plume_rise` (field `explicit_plume_rise`) (required)
  - `l_emission` (field `l_emission`) (required)
  - `l_points` (field `l_points`) (required)
  - `nemis` (field `nemis`) (required)

#### check_settings()

Check emission settings validity.

If this module configures emissions, ensure the runtime knows the
`emissions` folder is required so the job script can link/copy it.

#### emisnames *: [List](https://docs.python.org/3/library/typing.html#typing.List)[[str](https://docs.python.org/3/library/stdtypes.html#str)]*

#### emissions_instance *: [emissions](modular_dales.Emission.create_emis.md#modular_dales.Emission.create_emis.emissions) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### explicit_plume_rise *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

#### l_emission *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= True*

#### l_points *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

#### nemis *: [int](https://docs.python.org/3/library/functions.html#int)* *= 0*

#### point_sources *: [List](https://docs.python.org/3/library/typing.html#typing.List)[[EmissionPointSource](#modular_dales.Emission.emission.EmissionPointSource)]*

#### prepare_calculation()

Set up emission tracers and update namelist.

If no tracers or point sources have been added, emissions are disabled
and this method is a no-op.

#### sim *: simulation_module | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### tracers *: [List](https://docs.python.org/3/library/typing.html#typing.List)[[EmissionTracer](#modular_dales.Emission.emission.EmissionTracer)]*

#### write_files()

Write emission files to disk.

### *class* modular_dales.Emission.emission.EmissionPointSource(tracer_name: [str](https://docs.python.org/3/library/stdtypes.html#str), x_idx: [int](https://docs.python.org/3/library/functions.html#int), y_idx: [int](https://docs.python.org/3/library/functions.html#int), height: [float](https://docs.python.org/3/library/functions.html#float), temperature: [float](https://docs.python.org/3/library/functions.html#float), volume: [float](https://docs.python.org/3/library/functions.html#float), emission: [float](https://docs.python.org/3/library/functions.html#float), stack_exit_area: [float](https://docs.python.org/3/library/functions.html#float))

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Single emission point source tied to a tracer.

`tracer_name` must match the `name` of an `EmissionTracer` that is
added to the same `EmissionModule`.

#### emission *: [float](https://docs.python.org/3/library/functions.html#float)*

#### height *: [float](https://docs.python.org/3/library/functions.html#float)*

#### stack_exit_area *: [float](https://docs.python.org/3/library/functions.html#float)*

#### temperature *: [float](https://docs.python.org/3/library/functions.html#float)*

#### tracer_name *: [str](https://docs.python.org/3/library/stdtypes.html#str)*

#### volume *: [float](https://docs.python.org/3/library/functions.html#float)*

#### x_idx *: [int](https://docs.python.org/3/library/functions.html#int)*

#### y_idx *: [int](https://docs.python.org/3/library/functions.html#int)*

### *class* modular_dales.Emission.emission.EmissionTracer(name: [str](https://docs.python.org/3/library/stdtypes.html#str), long_name: [str](https://docs.python.org/3/library/stdtypes.html#str), unit: [str](https://docs.python.org/3/library/stdtypes.html#str), molar_mass: [float](https://docs.python.org/3/library/functions.html#float), lemis: [bool](https://docs.python.org/3/library/functions.html#bool) = False, lreact: [bool](https://docs.python.org/3/library/functions.html#bool) = False, ldep: [bool](https://docs.python.org/3/library/functions.html#bool) = False, lags: [bool](https://docs.python.org/3/library/functions.html#bool) = False)

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Single tracer configuration.

This closely follows `tracer_info` in `create_emis` so that it can be
converted directly when constructing the internal `emissions` object.

#### lags *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

#### ldep *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

#### lemis *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

#### long_name *: [str](https://docs.python.org/3/library/stdtypes.html#str)*

#### lreact *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

#### molar_mass *: [float](https://docs.python.org/3/library/functions.html#float)*

#### name *: [str](https://docs.python.org/3/library/stdtypes.html#str)*

#### unit *: [str](https://docs.python.org/3/library/stdtypes.html#str)*
