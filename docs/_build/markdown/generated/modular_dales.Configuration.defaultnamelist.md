# modular_dales.Configuration.defaultnamelist

Default configuration modules and base surface module.

### Classes

| [`DefaultNamelistModule`](#modular_dales.Configuration.defaultnamelist.DefaultNamelistModule)([sim, iexpnr, iinput, ...])   | Module that configures default RUN section.   |
|-----------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------|

### *class* modular_dales.Configuration.defaultnamelist.DefaultNamelistModule(sim: [dales_simulation](modular_dales.dales_simulation.md#modular_dales.dales_simulation) | [None](https://docs.python.org/3/library/constants.html#None) = None, iexpnr: [int](https://docs.python.org/3/library/functions.html#int) = 1, iinput: [int](https://docs.python.org/3/library/functions.html#int) = 2, ladaptive: [bool](https://docs.python.org/3/library/functions.html#bool) = True)

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
