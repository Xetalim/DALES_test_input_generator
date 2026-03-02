# modular_dales.Surface.LSM.LSMModule

### *class* modular_dales.Surface.LSM.LSMModule(sim: simulation_module | [None](https://docs.python.org/3/library/constants.html#None) = None, ps: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, z0mav: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, z0hav: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, albedoav: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, iinterp_t: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, iinterp_theta: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, kmax_soil: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = 4, dz_soil: [List](https://docs.python.org/3/library/typing.html#typing.List)[[float](https://docs.python.org/3/library/functions.html#float)] | [None](https://docs.python.org/3/library/constants.html#None) = None, lheterogeneous: [bool](https://docs.python.org/3/library/functions.html#bool) = True, land_use_modifications: LandUseModifications = <factory>, from_lcz: FromLCZ | [None](https://docs.python.org/3/library/constants.html#None) = None, from_ls2d: [FromLS2D](modular_dales.Atmosphere.ls2d_atmosphere.md#modular_dales.Atmosphere.ls2d_atmosphere.FromLS2D) | [None](https://docs.python.org/3/library/constants.html#None) = None, skin_temperature: UniformSkinTemperature | VaryingSkinTemperature | [None](https://docs.python.org/3/library/constants.html#None) = None, soil_temperature: UniformSoilTemperature | VaryingSoilTemperature | SoilTemperatureMoistureFromHarmonie | [None](https://docs.python.org/3/library/constants.html#None) = None, soil_moisture: UniformSoilMoisture | VaryingSoilMoisture | SoilTemperatureMoistureFromHarmonie | [None](https://docs.python.org/3/library/constants.html#None) = None)

Bases: [`SurfaceModule`](modular_dales.Surface.surface.md#modular_dales.Surface.surface.SurfaceModule)

Land Surface Model simulation module.

When added to a simulation, automatically enables LSM by setting isurf=11
in the namelist (unless overridden by a surface module like ConstantFluxesModule).

* **Parameters:**
  * **sim** – Parent simulation instance
  * **isurf** – Surface scheme selector (11 for LSM)
  * **ps** – Surface pressure (Pa)
  * **z0mav** – Momentum roughness length (m)
  * **z0hav** – Heat roughness length (m)
  * **albedoav** – Albedo (dimensionless)
  * **iinterp_t** – Interpolation method for temperature (1-4)
  * **iinterp_theta** – Interpolation method for potential temperature (1-4)
  * **dz_soil** – List of soil layer thicknesses (m)

## Namelist parameters

*Section* `DOMAIN`:
: - `kmax_soil` (field `kmax_soil`) (required)

*Section* `NAMLSM`:
: - `dz_soil` (field `dz_soil`) (required)
  - `iinterp_t` (field `iinterp_t`) (required)
  - `iinterp_theta` (field `iinterp_theta`) (required)
  - `lheterogeneous` (field `lheterogeneous`)
  - `nlu` (field `nlu`)

*Section* `NAMSURFACE`:
: - `albedoav` (field `albedoav`)
  - `isurf` (field `isurf`)
  - `ps` (field `ps`) (required)
  - `z0hav` (field `z0hav`) (required)
  - `z0mav` (field `z0mav`) (required)

#### \_\_init_\_(sim: simulation_module | [None](https://docs.python.org/3/library/constants.html#None) = None, ps: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, z0mav: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, z0hav: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, albedoav: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, iinterp_t: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, iinterp_theta: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, kmax_soil: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = 4, dz_soil: [List](https://docs.python.org/3/library/typing.html#typing.List)[[float](https://docs.python.org/3/library/functions.html#float)] | [None](https://docs.python.org/3/library/constants.html#None) = None, lheterogeneous: [bool](https://docs.python.org/3/library/functions.html#bool) = True, land_use_modifications: LandUseModifications = <factory>, from_lcz: FromLCZ | [None](https://docs.python.org/3/library/constants.html#None) = None, from_ls2d: [FromLS2D](modular_dales.Atmosphere.ls2d_atmosphere.md#modular_dales.Atmosphere.ls2d_atmosphere.FromLS2D) | [None](https://docs.python.org/3/library/constants.html#None) = None, skin_temperature: UniformSkinTemperature | VaryingSkinTemperature | [None](https://docs.python.org/3/library/constants.html#None) = None, soil_temperature: UniformSoilTemperature | VaryingSoilTemperature | SoilTemperatureMoistureFromHarmonie | [None](https://docs.python.org/3/library/constants.html#None) = None, soil_moisture: UniformSoilMoisture | VaryingSoilMoisture | SoilTemperatureMoistureFromHarmonie | [None](https://docs.python.org/3/library/constants.html#None) = None) → [None](https://docs.python.org/3/library/constants.html#None)

Initialize module with reference to parent simulation.

* **Parameters:**
  **sim** – Parent dales_simulation instance providing config, grid, paths, etc.
  Can be None if module will be added to simulation later.

### Methods

| [`__init__`](#modular_dales.Surface.LSM.LSMModule.__init__)(sim, ps, z0mav, z0hav, albedoav, ...)                   | Initialize module with reference to parent simulation.                      |
|---------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------|
| `apply_namelist_from_fields`()                                                                                      | Apply dataclass field values to the namelist using field metadata.          |
| [`apply_soil_temp_moisture_skin_temp`](#modular_dales.Surface.LSM.LSMModule.apply_soil_temp_moisture_skin_temp)()   |                                                                             |
| [`check_settings`](#modular_dales.Surface.LSM.LSMModule.check_settings)()                                           | Check LSM settings validity.                                                |
| [`do_config`](#modular_dales.Surface.LSM.LSMModule.do_config)()                                                     | Configure namelist and surface configuration for LSM.                       |
| [`exists_soil_temp_moisture_skin_temp`](#modular_dales.Surface.LSM.LSMModule.exists_soil_temp_moisture_skin_temp)() |                                                                             |
| `get_timedep_atmosphere_forcings`()                                                                                 | Optional hook: provide time-dependent forcing series for Atmosphere output. |
| `module_exists`(module_name)                                                                                        | Check if a module with the given name exists in the parent simulation.      |
| [`prepare_calculation`](#modular_dales.Surface.LSM.LSMModule.prepare_calculation)()                                 | Prepare data and setup for calculations.                                    |
| [`prepare_calculations`](#modular_dales.Surface.LSM.LSMModule.prepare_calculations)()                               | Set up LSM with land use types and soil configuration.                      |
| `retrieve_module`(module_name)                                                                                      | Helper to retrieve another module from the parent simulation by name.       |
| `set_nml_section`(section, key, value[, ...])                                                                       | Helper to set a single namelist section/key from within module methods.     |
| [`write_files`](#modular_dales.Surface.LSM.LSMModule.write_files)()                                                 | Write LSM input files and generate plots.                                   |

### Attributes

| [`albedoav`](#modular_dales.Surface.LSM.LSMModule.albedoav)                             |                                                       |
|-----------------------------------------------------------------------------------------|-------------------------------------------------------|
| [`dz_soil`](#modular_dales.Surface.LSM.LSMModule.dz_soil)                               |                                                       |
| `exp_id`                                                                                | Access experiment ID from parent simulation.          |
| [`from_lcz`](#modular_dales.Surface.LSM.LSMModule.from_lcz)                             |                                                       |
| [`from_ls2d`](#modular_dales.Surface.LSM.LSMModule.from_ls2d)                           |                                                       |
| `grid`                                                                                  | Access grid from parent simulation.                   |
| [`iinterp_t`](#modular_dales.Surface.LSM.LSMModule.iinterp_t)                           |                                                       |
| [`iinterp_theta`](#modular_dales.Surface.LSM.LSMModule.iinterp_theta)                   |                                                       |
| [`isurf`](#modular_dales.Surface.LSM.LSMModule.isurf)                                   |                                                       |
| [`kmax_soil`](#modular_dales.Surface.LSM.LSMModule.kmax_soil)                           |                                                       |
| [`lheterogeneous`](#modular_dales.Surface.LSM.LSMModule.lheterogeneous)                 |                                                       |
| [`lsm_writer`](#modular_dales.Surface.LSM.LSMModule.lsm_writer)                         |                                                       |
| [`nlu`](#modular_dales.Surface.LSM.LSMModule.nlu)                                       |                                                       |
| `nml`                                                                                   | Access namelist from parent simulation.               |
| `nml_docs`                                                                              | Access namelist documentation from parent simulation. |
| `output_path`                                                                           | Access output_path from parent simulation.            |
| [`ps`](#modular_dales.Surface.LSM.LSMModule.ps)                                         |                                                       |
| `required_folder_list`                                                                  | Access required_folder_list from parent simulation.   |
| [`sim`](#modular_dales.Surface.LSM.LSMModule.sim)                                       |                                                       |
| [`skin_temperature`](#modular_dales.Surface.LSM.LSMModule.skin_temperature)             |                                                       |
| [`slurb_module`](#modular_dales.Surface.LSM.LSMModule.slurb_module)                     |                                                       |
| [`soil_moisture`](#modular_dales.Surface.LSM.LSMModule.soil_moisture)                   |                                                       |
| [`soil_temperature`](#modular_dales.Surface.LSM.LSMModule.soil_temperature)             |                                                       |
| [`z0hav`](#modular_dales.Surface.LSM.LSMModule.z0hav)                                   |                                                       |
| [`z0mav`](#modular_dales.Surface.LSM.LSMModule.z0mav)                                   |                                                       |
| [`land_use_modifications`](#modular_dales.Surface.LSM.LSMModule.land_use_modifications) |                                                       |

#### albedoav *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### apply_soil_temp_moisture_skin_temp()

#### check_settings()

Check LSM settings validity.

#### do_config()

Configure namelist and surface configuration for LSM.

#### dz_soil *: [List](https://docs.python.org/3/library/typing.html#typing.List)[[float](https://docs.python.org/3/library/functions.html#float)] | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### exists_soil_temp_moisture_skin_temp()

#### from_lcz *: FromLCZ | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### from_ls2d *: [FromLS2D](modular_dales.Atmosphere.ls2d_atmosphere.md#modular_dales.Atmosphere.ls2d_atmosphere.FromLS2D) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### iinterp_t *: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### iinterp_theta *: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### isurf *: [int](https://docs.python.org/3/library/functions.html#int)* *= 11*

#### kmax_soil *: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None)* *= 4*

#### land_use_modifications *: LandUseModifications*

#### lheterogeneous *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= True*

#### lsm_writer *: LSM_output_dales | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### nlu *: [int](https://docs.python.org/3/library/functions.html#int)* *= 0*

#### prepare_calculation()

Prepare data and setup for calculations.

#### prepare_calculations()

Set up LSM with land use types and soil configuration.

#### ps *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### sim *: simulation_module | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### skin_temperature *: UniformSkinTemperature | VaryingSkinTemperature | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### slurb_module *: SLURBModule | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### soil_moisture *: UniformSoilMoisture | VaryingSoilMoisture | SoilTemperatureMoistureFromHarmonie | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### soil_temperature *: UniformSoilTemperature | VaryingSoilTemperature | SoilTemperatureMoistureFromHarmonie | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### write_files()

Write LSM input files and generate plots.

#### z0hav *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### z0mav *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*
