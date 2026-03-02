# modular_dales.Configuration

Public API for configuration modules.

These modules configure namelist sections for a DALES simulation.

### *class* modular_dales.Configuration.DefaultNamelistModule(sim: [dales_simulation](modular_dales.dales_simulation.md#modular_dales.dales_simulation) | [None](https://docs.python.org/3/library/constants.html#None) = None, iexpnr: [int](https://docs.python.org/3/library/functions.html#int) = 1, iinput: [int](https://docs.python.org/3/library/functions.html#int) = 2, ladaptive: [bool](https://docs.python.org/3/library/functions.html#bool) = True)

Bases: `simulation_module`

Module that configures default RUN section.

This module sets up default values for the RUN section only.
Domain configuration is handled by GridModule.
Other configuration modules handle their respective sections.

* **Parameters:**
  **sim** – Parent dales_simulation instance

## Namelist parameters

*Section* `RUN`:
: - `iexpnr` (field `iexpnr`) (required)
  - `iinput` (field `iinput`) (required)
  - `ladaptive` (field `ladaptive`) (required)

#### check_settings()

Check defaults validity.

#### do_config()

Configure RUN section defaults.

#### iexpnr *: [int](https://docs.python.org/3/library/functions.html#int)* *= 1*

#### iinput *: [int](https://docs.python.org/3/library/functions.html#int)* *= 2*

#### ladaptive *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= True*

#### prepare_calculation()

No calculation work needed.

#### sim *: [dales_simulation](modular_dales.dales_simulation.md#modular_dales.dales_simulation) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### write_files()

No files to write.

### *class* modular_dales.Configuration.EasyOutputModule(sim: [dales_simulation](modular_dales.dales_simulation.md#modular_dales.dales_simulation) | [None](https://docs.python.org/3/library/constants.html#None) = None, output_interval: [float](https://docs.python.org/3/library/functions.html#float) | [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)] | [None](https://docs.python.org/3/library/constants.html#None) = 60, enable_output: [bool](https://docs.python.org/3/library/functions.html#bool) = True)

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

### *class* modular_dales.Configuration.RunModule(sim: [dales_simulation](modular_dales.dales_simulation.md#modular_dales.dales_simulation) | [None](https://docs.python.org/3/library/constants.html#None) = None, iexpnr: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, runtime: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, ladaptive: [bool](https://docs.python.org/3/library/functions.html#bool) | [None](https://docs.python.org/3/library/constants.html#None) = None, dtmax: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None)

Bases: `simulation_module`

Run configuration module for DALES simulation.

Handles all run-related namelist parameters in the RUN section.

* **Parameters:**
  * **sim** – Parent dales_simulation instance
  * **iexpnr** – Experiment number
  * **runtime** – Total simulation time in seconds
  * **ladaptive** – Allow adaptive timestepping (highly recommended)
  * **dtmax** – Maximum timestep

## Namelist parameters

*Section* `RUN`:
: - `dtmax` (field `dtmax`)
  - `iexpnr` (field `iexpnr`)
  - `ladaptive` (field `ladaptive`)
  - `runtime` (field `runtime`)

#### check_settings()

Check run settings validity.

#### dtmax *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### iexpnr *: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### ladaptive *: [bool](https://docs.python.org/3/library/functions.html#bool) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### prepare_calculation()

Set up run-related namelist parameters.

#### runtime *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### sim *: [dales_simulation](modular_dales.dales_simulation.md#modular_dales.dales_simulation) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### write_files()

No files to write for run module.

### *class* modular_dales.Configuration.TimeModule(sim: [dales_simulation](modular_dales.dales_simulation.md#modular_dales.dales_simulation) | [None](https://docs.python.org/3/library/constants.html#None) = None, xtime: [float](https://docs.python.org/3/library/functions.html#float) = 0.0, xday: [int](https://docs.python.org/3/library/functions.html#int) = 265, xyear: [int](https://docs.python.org/3/library/functions.html#int) = 2025, enable_datetime: [bool](https://docs.python.org/3/library/functions.html#bool) = True, timezone: [int](https://docs.python.org/3/library/functions.html#int) = 0, startyear: [int](https://docs.python.org/3/library/functions.html#int) = 2025, startmonth: [int](https://docs.python.org/3/library/functions.html#int) = 1, startday: [int](https://docs.python.org/3/library/functions.html#int) = 1, runtime: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, trestart: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None)

Bases: `simulation_module`

Simulation time configuration module.

Handles time-related parameters in the RUN and DOMAIN sections.

* **Parameters:**
  * **sim** – Parent dales_simulation instance
  * **xtime** – Time within the day in seconds (0-86400)
  * **xday** – Day of year (1-365/366)
  * **xyear** – Year
  * **runtime** – Total runtime in seconds
  * **trestart** – Restart file write interval in seconds

## Namelist parameters

*Section* `DOMAIN`:
: - `xday` (field `xday`) (required)
  - `xtime` (field `xtime`) (required)
  - `xyear` (field `xyear`) (required)

*Section* `RUN`:
: - `runtime` (field `runtime`) (required)
  - `trestart` (field `trestart`)

*Section* `namdatetime`:
: - `l_datetime` (field `enable_datetime`) (required)
  - `startday` (field `startday`) (required)
  - `startmonth` (field `startmonth`) (required)
  - `startyear` (field `startyear`) (required)
  - `timezone` (field `timezone`)

#### check_settings()

Check time settings validity.

#### do_config()

Configure time-related namelist parameters.

#### enable_datetime *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= True*

#### prepare_calculation()

No calculation work needed.

#### runtime *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### sim *: [dales_simulation](modular_dales.dales_simulation.md#modular_dales.dales_simulation) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### startday *: [int](https://docs.python.org/3/library/functions.html#int)* *= 1*

#### startmonth *: [int](https://docs.python.org/3/library/functions.html#int)* *= 1*

#### startyear *: [int](https://docs.python.org/3/library/functions.html#int)* *= 2025*

#### timezone *: [int](https://docs.python.org/3/library/functions.html#int)* *= 0*

#### trestart *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### write_files()

No files to write.

#### xday *: [int](https://docs.python.org/3/library/functions.html#int)* *= 265*

#### xtime *: [float](https://docs.python.org/3/library/functions.html#float)* *= 0.0*

#### xyear *: [int](https://docs.python.org/3/library/functions.html#int)* *= 2025*

### Modules

| [`defaultnamelist`](modular_dales.Configuration.defaultnamelist.md#module-modular_dales.Configuration.defaultnamelist)   | Default configuration modules and base surface module.   |
|--------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------|
| [`output_modules`](modular_dales.Configuration.output_modules.md#module-modular_dales.Configuration.output_modules)      |                                                          |
| [`run_and_time`](modular_dales.Configuration.run_and_time.md#module-modular_dales.Configuration.run_and_time)            |                                                          |
