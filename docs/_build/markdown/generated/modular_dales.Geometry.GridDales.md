# modular_dales.Geometry.GridDales

### *class* modular_dales.Geometry.GridDales(sim: simulation_module | [None](https://docs.python.org/3/library/constants.html#None) = None, itot: [int](https://docs.python.org/3/library/functions.html#int) = None, jtot: [int](https://docs.python.org/3/library/functions.html#int) = None, kmax: [int](https://docs.python.org/3/library/functions.html#int) = None, xsize: [float](https://docs.python.org/3/library/functions.html#float) = None, ysize: [float](https://docs.python.org/3/library/functions.html#float) = None, kmax_soil: [int](https://docs.python.org/3/library/functions.html#int) = None, xlat: [float](https://docs.python.org/3/library/functions.html#float) = None, xlon: [float](https://docs.python.org/3/library/functions.html#float) = None, x0: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, y0: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, alpha: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, dz0: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, proj4: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None)

Bases: `simulation_module`

Grid configuration module for DALES simulation.

Handles all grid-related namelist parameters in the DOMAIN section.
Initialized with a complete domain_info dict containing all grid parameters.
This module should be added after DefaultNamelistModule and before other modules,
as it configures the grid that other modules depend on.

* **Parameters:**
  * **domain_info** – Dictionary containing domain configuration with keys:
    - itot, jtot, kmax: Grid dimensions
    - xsize, ysize: Domain physical dimensions
    - kmax_soil: Number of soil levels
    - xlat, xlon: Latitude/longitude
    - x0, y0: Grid origin
    - alpha: Vertical stretch factor
    - dz0: Initial vertical spacing
    - proj4 (optional): Projection string
  * **sim** – Parent dales_simulation instance

## Namelist parameters

*Section* `DOMAIN`:
: - `itot` (field `itot`) (required)
  - `jtot` (field `jtot`) (required)
  - `kmax` (field `kmax`) (required)
  - `kmax_soil` (field `kmax_soil`) (required)
  - `x0` (field `x0`)
  - `xlat` (field `xlat`)
  - `xlon` (field `xlon`)
  - `xsize` (field `xsize`) (required)
  - `y0` (field `y0`)
  - `ysize` (field `ysize`) (required)

#### \_\_init_\_(sim: simulation_module | [None](https://docs.python.org/3/library/constants.html#None) = None, itot: [int](https://docs.python.org/3/library/functions.html#int) = None, jtot: [int](https://docs.python.org/3/library/functions.html#int) = None, kmax: [int](https://docs.python.org/3/library/functions.html#int) = None, xsize: [float](https://docs.python.org/3/library/functions.html#float) = None, ysize: [float](https://docs.python.org/3/library/functions.html#float) = None, kmax_soil: [int](https://docs.python.org/3/library/functions.html#int) = None, xlat: [float](https://docs.python.org/3/library/functions.html#float) = None, xlon: [float](https://docs.python.org/3/library/functions.html#float) = None, x0: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, y0: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, alpha: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, dz0: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, proj4: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None) → [None](https://docs.python.org/3/library/constants.html#None)

Initialize module with reference to parent simulation.

* **Parameters:**
  **sim** – Parent dales_simulation instance providing config, grid, paths, etc.
  Can be None if module will be added to simulation later.

### Methods

| `__init__`([sim, itot, jtot, kmax, xsize, ...])   | Initialize module with reference to parent simulation.                      |
|---------------------------------------------------|-----------------------------------------------------------------------------|
| `apply_namelist_from_fields`()                    | Apply dataclass field values to the namelist using field metadata.          |
| `as_dic`()                                        |                                                                             |
| `as_openbc`()                                     |                                                                             |
| `check_settings`()                                | Check grid settings validity.                                               |
| `do_config`()                                     | Configure DOMAIN section from domain_info and create GridDales.             |
| `get_timedep_atmosphere_forcings`()               | Optional hook: provide time-dependent forcing series for Atmosphere output. |
| `module_exists`(module_name)                      | Check if a module with the given name exists in the parent simulation.      |
| `prepare_calculation`()                           | No additional preparation needed.                                           |
| `retrieve_module`(module_name)                    | Helper to retrieve another module from the parent simulation by name.       |
| `set_nml_section`(section, key, value[, ...])     | Helper to set a single namelist section/key from within module methods.     |
| `write_files`()                                   | No files to write for grid module.                                          |

### Attributes

| `alpha`                |                                                       |
|------------------------|-------------------------------------------------------|
| `dz0`                  |                                                       |
| `exp_id`               | Access experiment ID from parent simulation.          |
| `grid`                 | Access grid from parent simulation.                   |
| `itot`                 |                                                       |
| `jtot`                 |                                                       |
| `kmax`                 |                                                       |
| `kmax_soil`            |                                                       |
| `nml`                  | Access namelist from parent simulation.               |
| `nml_docs`             | Access namelist documentation from parent simulation. |
| `output_path`          | Access output_path from parent simulation.            |
| `proj4`                |                                                       |
| `required_folder_list` | Access required_folder_list from parent simulation.   |
| `sim`                  |                                                       |
| `x0`                   |                                                       |
| `xlat`                 |                                                       |
| `xlon`                 |                                                       |
| `xsize`                |                                                       |
| `y0`                   |                                                       |
| `ysize`                |                                                       |

#### alpha *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### as_dic()

#### as_openbc()

#### check_settings()

Check grid settings validity.

#### do_config()

Configure DOMAIN section from domain_info and create GridDales.

#### dz0 *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### itot *: [int](https://docs.python.org/3/library/functions.html#int)* *= None*

#### jtot *: [int](https://docs.python.org/3/library/functions.html#int)* *= None*

#### kmax *: [int](https://docs.python.org/3/library/functions.html#int)* *= None*

#### kmax_soil *: [int](https://docs.python.org/3/library/functions.html#int)* *= None*

#### prepare_calculation()

No additional preparation needed.

#### proj4 *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### sim *: simulation_module | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### write_files()

No files to write for grid module.

#### x0 *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### xlat *: [float](https://docs.python.org/3/library/functions.html#float)* *= None*

#### xlon *: [float](https://docs.python.org/3/library/functions.html#float)* *= None*

#### xsize *: [float](https://docs.python.org/3/library/functions.html#float)* *= None*

#### y0 *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### ysize *: [float](https://docs.python.org/3/library/functions.html#float)* *= None*
