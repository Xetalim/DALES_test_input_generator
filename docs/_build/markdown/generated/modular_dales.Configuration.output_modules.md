# modular_dales.Configuration.output_modules

### Classes

| [`CheckSimulationModule`](#modular_dales.Configuration.output_modules.CheckSimulationModule)([sim, check_interval, ...])   | Module for checking simulation validity and configuration.   |
|----------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------|
| [`EasyOutputModule`](#modular_dales.Configuration.output_modules.EasyOutputModule)([sim, output_interval, ...])            | Output configuration module for DALES simulation.            |

### *class* modular_dales.Configuration.output_modules.CheckSimulationModule(sim: [dales_simulation](modular_dales.dales_simulation.md#modular_dales.dales_simulation) | [None](https://docs.python.org/3/library/constants.html#None) = None, check_interval: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = 60, stop_on_invalid: [bool](https://docs.python.org/3/library/functions.html#bool) = False, check_tendencies: [bool](https://docs.python.org/3/library/functions.html#bool) = False)

Bases: `simulation_module`

Module for checking simulation validity and configuration.

This module provides functionality to monitor and validate a DALES simulation
during runtime. It allows output at an interval of check_interval seconds, at which DALES prints
the currrent simulation time and ETA. It can also check for invalid values, verify tendencies, and
optionally stop the simulation if invalid conditions are detected.

#### check_interval

Time interval (in seconds) between output.
Default is 60 seconds. Mapped to NAMCHECKSIM:tcheck in Fortran namelist.
If 0, DALES prints every timestep.

* **Type:**
  Optional[[int](https://docs.python.org/3/library/functions.html#int)]

#### stop_on_invalid

If True, stops the simulation when invalid values are detected.
Default is False. Mapped to NAMCHECKSIM:lstop in Fortran namelist.

* **Type:**
  [bool](https://docs.python.org/3/library/functions.html#bool)

#### check_tendencies

If True, validates model tendencies during runtime checks.
Default is False. Mapped to NAMCHECKSIM:lchecktend in Fortran namelist.

* **Type:**
  [bool](https://docs.python.org/3/library/functions.html#bool)

## Namelist parameters

*Section* `NAMCHECKSIM`:
: - `lchecktend` (field `check_tendencies`) (required)
  - `lstop` (field `stop_on_invalid`) (required)
  - `tcheck` (field `check_interval`) (required)

#### check_interval *: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None)* *= 60*

#### check_settings()

No additional checks needed.

#### check_tendencies *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

#### do_config()

No configuration needed.

#### prepare_calculation()

Check simulation settings and configuration.

#### sim *: [dales_simulation](modular_dales.dales_simulation.md#modular_dales.dales_simulation) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### stop_on_invalid *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

#### write_files()

No files to write.

### *class* modular_dales.Configuration.output_modules.EasyOutputModule(sim: [dales_simulation](modular_dales.dales_simulation.md#modular_dales.dales_simulation) | [None](https://docs.python.org/3/library/constants.html#None) = None, output_interval: [float](https://docs.python.org/3/library/functions.html#float) | [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)] | [None](https://docs.python.org/3/library/constants.html#None) = 60, enable_output: [bool](https://docs.python.org/3/library/functions.html#bool) = True)

Bases: `simulation_module`

Output configuration module for DALES simulation.
Outputs fielddumps, cape, cross-sections, and statistics at specified intervals.
When added to a simulation, automatically updates the corresponding namelist parameters in the appropriate sections.

* **Parameters:**
  * **sim** – Parent dales_simulation instance
  * **output_interval** – Time interval(s) for output (can be a single value or list for different output types)
  * **enable_output** – Whether to enable output (controls lfielddump, lcape, lcross, lstat in namelist)

## Namelist parameters

*Section* `nambudget`:
: - `dtav` (field `output_interval`) (required)
  - `lbudget` (field `enable_output`) (required)
  - `timeav` (field `output_interval`) (required)

*Section* `namcape`:
: - `dtav` (field `output_interval`) (required)
  - `lcape` (field `enable_output`) (required)

*Section* `namcrosssection`:
: - `dtav` (field `output_interval`) (required)
  - `lcross` (field `enable_output`) (required)

*Section* `namfielddump`:
: - `dtav` (field `output_interval`) (required)
  - `lfielddump` (field `enable_output`) (required)

*Section* `namgenstat`:
: - `dtav` (field `output_interval`) (required)
  - `lstat` (field `enable_output`) (required)
  - `timeav` (field `output_interval`) (required)

*Section* `namlsmcrosssection`:
: - `dtav` (field `output_interval`) (required)
  - `lcross` (field `enable_output`) (required)

*Section* `namtimestat`:
: - `dtav` (field `output_interval`) (required)
  - `ltimestat` (field `enable_output`) (required)

#### check_settings()

Check output settings validity.

#### do_config()

Configure output-related namelist parameters.

#### enable_output *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= True*

#### output_interval *: [float](https://docs.python.org/3/library/functions.html#float) | [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)] | [None](https://docs.python.org/3/library/constants.html#None)* *= 60*

#### prepare_calculation()

No preparation needed.

#### sim *: [dales_simulation](modular_dales.dales_simulation.md#modular_dales.dales_simulation) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### write_files()

Write output files if needed.
