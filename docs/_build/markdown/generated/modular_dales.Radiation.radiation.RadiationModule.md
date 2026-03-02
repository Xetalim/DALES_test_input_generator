# modular_dales.Radiation.radiation.RadiationModule

### *class* modular_dales.Radiation.radiation.RadiationModule(sim: simulation_module | [None](https://docs.python.org/3/library/constants.html#None) = None, iradiation: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, ssa: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, laero: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, ide: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, lCnstZenith: [bool](https://docs.python.org/3/library/functions.html#bool) | [None](https://docs.python.org/3/library/constants.html#None) = None, ioverlap: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, inflglw: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, iceflglw: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, liqflglw: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, inflgsw: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, iceflgsw: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, liqflgsw: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, iyear: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, ocean: [bool](https://docs.python.org/3/library/functions.html#bool) | [None](https://docs.python.org/3/library/constants.html#None) = None, nbatch: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, usepade: [bool](https://docs.python.org/3/library/functions.html#bool) | [None](https://docs.python.org/3/library/constants.html#None) = None, doclearsky: [bool](https://docs.python.org/3/library/functions.html#bool) | [None](https://docs.python.org/3/library/constants.html#None) = None, timerad: [int](https://docs.python.org/3/library/functions.html#int) = 60)

Bases: `simulation_module`

Radiation simulation module.

Possible values for iradiation:
: 0: no radiation
  1: full radiation
  2: parameterized radiation
  3: simple surface radiation for land surface model
  4: RRTMG radiation
  5: RTE-RRTMGP radiation
  10: user specified radiation

* **Parameters:**
  * **sim** – Parent dales_simulation instance
  * **iradiation** – Radiation type
  * **ssa** – Representative single scattering albedo (0 <= x <= 1)
  * **ide** – Scalar field used as aerosols if laero set to .true.
  * **laero** – .true. for aerosols, .false. for clouds
  * **lCnstZenith** – Switch to apply a fixed solar zenith angle
  * **ioverlap** – Flag for cloud overlap method
  * **inflglw** – Flag for RRTMG longwave input
  * **iceflglw** – Flag for ice particle specification in longwave
  * **liqflglw** – Flag for effect of liquid water in longwave
  * **inflgsw** – Flag for RRTMG shortwave input
  * **iceflgsw** – Flag for ice particle specification in shortwave
  * **liqflgsw** – Flag for effect of liquid water in shortwave
  * **iyear** – Year of the simulation
  * **ocean** – Switch to calculate radiation over ocean
  * **nbatch** – Number of batch of vertical columns sent to RTE-RRTMGP routines
  * **usepade** – Use Pade coefficients for cloud optical properties instead of lookup tables
  * **doclearsky** – Use clear sky radiation in the calculation

## Namelist parameters

*Section* `NAMDE`:
: - `ide` (field `ide`)
  - `laero` (field `laero`)
  - `ssa` (field `ssa`)

*Section* `NAMRADIATION`:
: - `iceflglw` (field `iceflglw`)
  - `iceflgsw` (field `iceflgsw`)
  - `inflglw` (field `inflglw`)
  - `inflgsw` (field `inflgsw`)
  - `ioverlap` (field `ioverlap`)
  - `iyear` (field `iyear`)
  - `lCnstZenith` (field `lCnstZenith`)
  - `liqflglw` (field `liqflglw`)
  - `liqflgsw` (field `liqflgsw`)
  - `ocean` (field `ocean`)

*Section* `NAMRTERRTMGP`:
: - `doclearsky` (field `doclearsky`)
  - `nbatch` (field `nbatch`)
  - `usepade` (field `usepade`)

*Section* `PHYSICS`:
: - `IRADIATION` (field `iradiation`)
  - `timerad` (field `timerad`) (required)

#### \_\_init_\_(sim: simulation_module | [None](https://docs.python.org/3/library/constants.html#None) = None, iradiation: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, ssa: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, laero: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, ide: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, lCnstZenith: [bool](https://docs.python.org/3/library/functions.html#bool) | [None](https://docs.python.org/3/library/constants.html#None) = None, ioverlap: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, inflglw: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, iceflglw: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, liqflglw: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, inflgsw: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, iceflgsw: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, liqflgsw: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, iyear: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, ocean: [bool](https://docs.python.org/3/library/functions.html#bool) | [None](https://docs.python.org/3/library/constants.html#None) = None, nbatch: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, usepade: [bool](https://docs.python.org/3/library/functions.html#bool) | [None](https://docs.python.org/3/library/constants.html#None) = None, doclearsky: [bool](https://docs.python.org/3/library/functions.html#bool) | [None](https://docs.python.org/3/library/constants.html#None) = None, timerad: [int](https://docs.python.org/3/library/functions.html#int) = 60) → [None](https://docs.python.org/3/library/constants.html#None)

Initialize module with reference to parent simulation.

* **Parameters:**
  **sim** – Parent dales_simulation instance providing config, grid, paths, etc.
  Can be None if module will be added to simulation later.

### Methods

| [`__init__`](#modular_dales.Radiation.radiation.RadiationModule.__init__)([sim, iradiation, ssa, laero, ide, ...])   | Initialize module with reference to parent simulation.                      |
|----------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------|
| `apply_namelist_from_fields`()                                                                                       | Apply dataclass field values to the namelist using field metadata.          |
| [`check_settings`](#modular_dales.Radiation.radiation.RadiationModule.check_settings)()                              | Validate constant fluxes settings.                                          |
| [`do_config`](#modular_dales.Radiation.radiation.RadiationModule.do_config)()                                        | Ensure radiation configuration is set.                                      |
| `get_timedep_atmosphere_forcings`()                                                                                  | Optional hook: provide time-dependent forcing series for Atmosphere output. |
| `module_exists`(module_name)                                                                                         | Check if a module with the given name exists in the parent simulation.      |
| [`prepare_calculation`](#modular_dales.Radiation.radiation.RadiationModule.prepare_calculation)()                    | No additional preparation needed.                                           |
| `retrieve_module`(module_name)                                                                                       | Helper to retrieve another module from the parent simulation by name.       |
| `set_nml_section`(section, key, value[, ...])                                                                        | Helper to set a single namelist section/key from within module methods.     |
| [`write_files`](#modular_dales.Radiation.radiation.RadiationModule.write_files)()                                    | Write output files for this module.                                         |

### Attributes

| [`doclearsky`](#modular_dales.Radiation.radiation.RadiationModule.doclearsky)         |                                                       |
|---------------------------------------------------------------------------------------|-------------------------------------------------------|
| `exp_id`                                                                              | Access experiment ID from parent simulation.          |
| `grid`                                                                                | Access grid from parent simulation.                   |
| [`iceflglw`](#modular_dales.Radiation.radiation.RadiationModule.iceflglw)             |                                                       |
| [`iceflgsw`](#modular_dales.Radiation.radiation.RadiationModule.iceflgsw)             |                                                       |
| [`ide`](#modular_dales.Radiation.radiation.RadiationModule.ide)                       |                                                       |
| [`inflglw`](#modular_dales.Radiation.radiation.RadiationModule.inflglw)               |                                                       |
| [`inflgsw`](#modular_dales.Radiation.radiation.RadiationModule.inflgsw)               |                                                       |
| [`ioverlap`](#modular_dales.Radiation.radiation.RadiationModule.ioverlap)             |                                                       |
| [`iradiation`](#modular_dales.Radiation.radiation.RadiationModule.iradiation)         |                                                       |
| [`iyear`](#modular_dales.Radiation.radiation.RadiationModule.iyear)                   |                                                       |
| [`lCnstZenith`](#modular_dales.Radiation.radiation.RadiationModule.lCnstZenith)       |                                                       |
| [`laero`](#modular_dales.Radiation.radiation.RadiationModule.laero)                   |                                                       |
| [`liqflglw`](#modular_dales.Radiation.radiation.RadiationModule.liqflglw)             |                                                       |
| [`liqflgsw`](#modular_dales.Radiation.radiation.RadiationModule.liqflgsw)             |                                                       |
| [`nbatch`](#modular_dales.Radiation.radiation.RadiationModule.nbatch)                 |                                                       |
| `nml`                                                                                 | Access namelist from parent simulation.               |
| `nml_docs`                                                                            | Access namelist documentation from parent simulation. |
| [`ocean`](#modular_dales.Radiation.radiation.RadiationModule.ocean)                   |                                                       |
| `output_path`                                                                         | Access output_path from parent simulation.            |
| `required_folder_list`                                                                | Access required_folder_list from parent simulation.   |
| [`sim`](#modular_dales.Radiation.radiation.RadiationModule.sim)                       |                                                       |
| [`ssa`](#modular_dales.Radiation.radiation.RadiationModule.ssa)                       |                                                       |
| [`surface_module`](#modular_dales.Radiation.radiation.RadiationModule.surface_module) |                                                       |
| [`timerad`](#modular_dales.Radiation.radiation.RadiationModule.timerad)               |                                                       |
| [`usepade`](#modular_dales.Radiation.radiation.RadiationModule.usepade)               |                                                       |

#### check_settings()

Validate constant fluxes settings.

#### do_config()

Ensure radiation configuration is set.

#### doclearsky *: [bool](https://docs.python.org/3/library/functions.html#bool) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### iceflglw *: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### iceflgsw *: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### ide *: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### inflglw *: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### inflgsw *: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### ioverlap *: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### iradiation *: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### iyear *: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### lCnstZenith *: [bool](https://docs.python.org/3/library/functions.html#bool) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### laero *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### liqflglw *: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### liqflgsw *: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### nbatch *: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### ocean *: [bool](https://docs.python.org/3/library/functions.html#bool) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### prepare_calculation()

No additional preparation needed.

#### sim *: simulation_module | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### ssa *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### surface_module *: [SurfaceModule](modular_dales.Surface.surface.md#modular_dales.Surface.surface.SurfaceModule) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### timerad *: [int](https://docs.python.org/3/library/functions.html#int)* *= 60*

#### usepade *: [bool](https://docs.python.org/3/library/functions.html#bool) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### write_files()

Write output files for this module.
