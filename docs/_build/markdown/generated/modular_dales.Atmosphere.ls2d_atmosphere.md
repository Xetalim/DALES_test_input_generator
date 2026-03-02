# modular_dales.Atmosphere.ls2d_atmosphere

### Functions

| [`retry_exponential`](#modular_dales.Atmosphere.ls2d_atmosphere.retry_exponential)(func, \*[, max_attempts, ...])   | Retry a function with exponential backoff.   |
|---------------------------------------------------------------------------------------------------------------------|----------------------------------------------|

### Classes

| [`FromLS2D`](#modular_dales.Atmosphere.ls2d_atmosphere.FromLS2D)()                                                | Marker class to enable LS2D-driven soil/roughness in LSM.   |
|-------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------|
| [`LS2DAtmosphereModule`](#modular_dales.Atmosphere.ls2d_atmosphere.LS2DAtmosphereModule)([sim, central_lat, ...]) | Atmosphere forcing module using LS2D ERA5 processing.       |

### *class* modular_dales.Atmosphere.ls2d_atmosphere.FromLS2D

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Marker class to enable LS2D-driven soil/roughness in LSM.
Also enables injection of LS2D time series into atmosphere.

When added to `LSMModule` or `TimeDependentModule` (`lsm += FromLS2D()`), the
module will, if an [`LS2DAtmosphereModule`](#modular_dales.Atmosphere.ls2d_atmosphere.LS2DAtmosphereModule) is present,
override soil temperature/moisture and soil type index from LS2D
and also set the bulk roughness lengths `z0mav` and `z0hav`
from LS2D time-mean values.

### *class* modular_dales.Atmosphere.ls2d_atmosphere.LS2DAtmosphereModule(sim: simulation_module | [None](https://docs.python.org/3/library/constants.html#None) = None, central_lat: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, central_lon: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, area_size: [float](https://docs.python.org/3/library/functions.html#float) = 1.0, case_name: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'ls2d_case', era5_path: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None, start_date: [datetime](https://docs.python.org/3/library/datetime.html#datetime.datetime) | [None](https://docs.python.org/3/library/constants.html#None) = None, end_date: [datetime](https://docs.python.org/3/library/datetime.html#datetime.datetime) | [None](https://docs.python.org/3/library/constants.html#None) = None, write_log: [bool](https://docs.python.org/3/library/functions.html#bool) = True, data_source: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'CDS', n_av: [int](https://docs.python.org/3/library/functions.html#int) = 0, method: [str](https://docs.python.org/3/library/stdtypes.html#str) = '2nd', init_tke: [float](https://docs.python.org/3/library/functions.html#float) = 0.1)

Bases: `simulation_module`

Atmosphere forcing module using LS2D ERA5 processing.

* Uses LS2D to download and process ERA5 data.
* Calls `era.calculate_forcings` in LS2D.
* Interpolates the resulting forcings onto the DALES
  vertical grid given by `GridDales` (`grid.zt`).
  > * Provides large-scale forcings and base profiles that can be
  >   : injected into `AtmosphereModule` and surface modules.

#### area_size *: [float](https://docs.python.org/3/library/functions.html#float)* *= 1.0*

#### case_name *: [str](https://docs.python.org/3/library/stdtypes.html#str)* *= 'ls2d_case'*

#### central_lat *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### central_lon *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### check_settings()

Basic validation of configuration before running LS2D.

#### data_source *: [str](https://docs.python.org/3/library/stdtypes.html#str)* *= 'CDS'*

#### end_date *: [datetime](https://docs.python.org/3/library/datetime.html#datetime.datetime) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### era5_path *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### init_tke *: [float](https://docs.python.org/3/library/functions.html#float)* *= 0.1*

#### les_input *: [Any](https://docs.python.org/3/library/typing.html#typing.Any)* *= None*

#### method *: [str](https://docs.python.org/3/library/stdtypes.html#str)* *= '2nd'*

#### n_av *: [int](https://docs.python.org/3/library/functions.html#int)* *= 0*

#### prepare_calculation()

Run LS2D, build `les_input` on `GridDales.zt` and cache forcings.

The resulting time-dependent fields are stored internally in a
form that can be used both for writing `forcings.<exp_id>.nc`
and for injecting time series into other modules.

* Surface scalars are stored as 1D arrays with length `nt`.
* Profile variables are stored as 2D arrays of shape
  : `(nz, nt)` where first axis is height and second axis is time.

The time axis is taken from `les_input.time_sec` and cached in
`_times_with_zero`.

#### sim *: simulation_module | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### start_date *: [datetime](https://docs.python.org/3/library/datetime.html#datetime.datetime) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### write_files()

Write LS2D-derived files.

This module no longer writes the `init.<exp_id>.nc` file
itself; base profiles are injected into an
`AtmosphereModule` (when present), which then writes the
init file. The radiation background is
written to `backrad.inp.<exp_id>.nc`.

#### write_log *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= True*

### modular_dales.Atmosphere.ls2d_atmosphere.retry_exponential(func: [Callable](https://docs.python.org/3/library/typing.html#typing.Callable), , max_attempts: [int](https://docs.python.org/3/library/functions.html#int) = 5, max_total_time: [float](https://docs.python.org/3/library/functions.html#float) = 60.0, base_delay: [float](https://docs.python.org/3/library/functions.html#float) = 1.0, max_delay: [float](https://docs.python.org/3/library/functions.html#float) = 30.0)

Retry a function with exponential backoff.

* **Parameters:**
  * **func** – Function to execute
  * **max_attempts** – Maximum number of attempts
  * **max_total_time** – Maximum total retry time (seconds)
  * **base_delay** – Initial delay in seconds
  * **max_delay** – Maximum delay between retries
