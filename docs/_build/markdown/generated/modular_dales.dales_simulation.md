# modular_dales.dales_simulation

### *class* modular_dales.dales_simulation(case_name: [str](https://docs.python.org/3/library/stdtypes.html#str), machine_conf)

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Main simulation class that manages multiple simulation modules.

Usage:
: sim = dales_simulation(case_name, machine_conf)
  sim += DefaultNamelistModule()
  sim += GridModule()
  sim += LSMModule()
  # … add other modules …
  sim._initialize_modules()  # Initialize grid and validation
  sim.do_config()             # Configure namelist
  sim.check_settings()        # Validate settings
  sim.prepare_all_calculations()  # Prepare calculations
  sim.write_module_files()       # Write output files

#### \_\_init_\_(case_name: [str](https://docs.python.org/3/library/stdtypes.html#str), machine_conf)

Initialize DALES simulation with case name and machine settings.

Only reads configuration at this stage. All work is deferred to
explicit method calls (_initialize_modules, do_config, check_settings,
prepare_all_calculations, write_module_files).

* **Parameters:**
  * **case_name** – Name of the simulation case
  * **machine_conf** – Machine configuration dictionary

### Methods

| [`__init__`](#modular_dales.dales_simulation.__init__)(case_name, machine_conf)                 | Initialize DALES simulation with case name and machine settings.          |
|-------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------|
| [`apply_job_configuration`](#modular_dales.dales_simulation.apply_job_configuration)()          | Apply job configuration and set up file transfers.                        |
| [`check_settings`](#modular_dales.dales_simulation.check_settings)()                            | Check settings and validate configuration for all modules.                |
| [`do_config`](#modular_dales.dales_simulation.do_config)()                                      | Configure all modules.                                                    |
| [`init_output_folder`](#modular_dales.dales_simulation.init_output_folder)()                    | Initialize output folder and paths.                                       |
| [`load_sim_from_yaml`](#modular_dales.dales_simulation.load_sim_from_yaml)(txt[, machine_conf]) | Load a YAML file into a dales_simulation instance.                        |
| [`module_exists`](#modular_dales.dales_simulation.module_exists)(module_name)                   | Check if a module with the given name exists.                             |
| [`prepare_all_calculations`](#modular_dales.dales_simulation.prepare_all_calculations)()        | Prepare all modules for calculations.                                     |
| [`retrieve_module`](#modular_dales.dales_simulation.retrieve_module)(module_name)               | Helper to retrieve another module by name.                                |
| [`save_sim_to_yaml`](#modular_dales.dales_simulation.save_sim_to_yaml)()                        | Save a dales_simulation instance to YAML including case_name and modules. |
| [`setup_module_links`](#modular_dales.dales_simulation.setup_module_links)()                    |                                                                           |
| [`sim_preprocessing_pipeline`](#modular_dales.dales_simulation.sim_preprocessing_pipeline)()    | Run the full modular pipeline to produce input files.                     |
| [`write_module_files`](#modular_dales.dales_simulation.write_module_files)()                    | Call write_files on all registered modules.                               |
| [`write_simulation_files`](#modular_dales.dales_simulation.write_simulation_files)()            | Write namelist and required files to case directory.                      |

#### apply_job_configuration()

Apply job configuration and set up file transfers.

#### check_settings()

Check settings and validate configuration for all modules.

Call this explicitly after do_config() but before prepare_all_calculations().

#### do_config()

Configure all modules.

Call this explicitly after adding all modules but before check_settings().

#### init_output_folder()

Initialize output folder and paths.

This method is called explicitly at the start of the workflow.

#### *static* load_sim_from_yaml(txt, machine_conf=None) → [dales_simulation](#modular_dales.dales_simulation)

Load a YAML file into a dales_simulation instance.

* **Parameters:**
  * **path** – str or Path to YAML file
  * **machine_conf**
* **Returns:**
  populated dales_simulation
* **Return type:**
  [sim](modular_dales.Atmosphere.md#modular_dales.Atmosphere.AtmosphereModule.sim)

#### module_exists(module_name: [str](https://docs.python.org/3/library/stdtypes.html#str) | simulation_module | [type](https://docs.python.org/3/library/functions.html#type)) → [bool](https://docs.python.org/3/library/functions.html#bool)

Check if a module with the given name exists.

#### prepare_all_calculations()

Prepare all modules for calculations.

Call this explicitly after check_settings().

#### retrieve_module(module_name: [str](https://docs.python.org/3/library/stdtypes.html#str) | simulation_module | [type](https://docs.python.org/3/library/functions.html#type)) → simulation_module

Helper to retrieve another module by name.

#### save_sim_to_yaml() → [str](https://docs.python.org/3/library/stdtypes.html#str)

Save a dales_simulation instance to YAML including case_name and modules.

* **Parameters:**
  **sim** – dales_simulation instance

#### setup_module_links()

#### sim_preprocessing_pipeline() → [None](https://docs.python.org/3/library/constants.html#None)

Run the full modular pipeline to produce input files.

This matches the intended workflow:
- init_output_folder
- setup_module_links
- do_config
- check_settings
- prepare_all_calculations
- apply_job_configuration
- write_module_files (per-module outputs)
- write_simulation_files (namelist and required files)

#### write_module_files()

Call write_files on all registered modules.

#### write_simulation_files()

Write namelist and required files to case directory.
