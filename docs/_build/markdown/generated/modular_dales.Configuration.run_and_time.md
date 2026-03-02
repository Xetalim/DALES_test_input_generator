# modular_dales.Configuration.run_and_time

### Classes

| [`RunModule`](#modular_dales.Configuration.run_and_time.RunModule)([sim, iexpnr, runtime, ladaptive, ...])   | Run configuration module for DALES simulation.   |
|--------------------------------------------------------------------------------------------------------------|--------------------------------------------------|
| [`TimeModule`](#modular_dales.Configuration.run_and_time.TimeModule)([sim, xtime, xday, xyear, ...])         | Simulation time configuration module.            |

### *class* modular_dales.Configuration.run_and_time.RunModule(sim: [dales_simulation](modular_dales.dales_simulation.md#modular_dales.dales_simulation) | [None](https://docs.python.org/3/library/constants.html#None) = None, iexpnr: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, runtime: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, ladaptive: [bool](https://docs.python.org/3/library/functions.html#bool) | [None](https://docs.python.org/3/library/constants.html#None) = None, dtmax: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None)

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

### *class* modular_dales.Configuration.run_and_time.TimeModule(sim: [dales_simulation](modular_dales.dales_simulation.md#modular_dales.dales_simulation) | [None](https://docs.python.org/3/library/constants.html#None) = None, xtime: [float](https://docs.python.org/3/library/functions.html#float) = 0.0, xday: [int](https://docs.python.org/3/library/functions.html#int) = 265, xyear: [int](https://docs.python.org/3/library/functions.html#int) = 2025, enable_datetime: [bool](https://docs.python.org/3/library/functions.html#bool) = True, timezone: [int](https://docs.python.org/3/library/functions.html#int) = 0, startyear: [int](https://docs.python.org/3/library/functions.html#int) = 2025, startmonth: [int](https://docs.python.org/3/library/functions.html#int) = 1, startday: [int](https://docs.python.org/3/library/functions.html#int) = 1, runtime: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, trestart: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None)

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
