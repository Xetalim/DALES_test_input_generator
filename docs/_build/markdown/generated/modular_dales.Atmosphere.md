# modular_dales.Atmosphere

### *class* modular_dales.Atmosphere.AtmosphereModule(sim: simulation_module | [None](https://docs.python.org/3/library/constants.html#None) = None, shaped_profiles: [List](https://docs.python.org/3/library/typing.html#typing.List)[[AtmosphericProfile](modular_dales.Atmosphere.atmosphere.md#modular_dales.Atmosphere.atmosphere.AtmosphericProfile)] = <factory>, interpolated_profiles: [List](https://docs.python.org/3/library/typing.html#typing.List)[[InterpolatedProfile](modular_dales.Atmosphere.atmosphere.md#modular_dales.Atmosphere.atmosphere.InterpolatedProfile)] = <factory>, timed_profiles: [List](https://docs.python.org/3/library/typing.html#typing.List)[[TimedAtmosphereProfile](modular_dales.Atmosphere.atmosphere.md#modular_dales.Atmosphere.atmosphere.TimedAtmosphereProfile)] = <factory>)

Bases: `simulation_module`

Atmosphere module for profile setup.

* **Parameters:**
  * **sim** – Parent simulation instance
  * **shaped_profiles** – Mapping of variable name to shaped profile configuration
  * **interpolated_profiles** – Mapping of variable name to interpolated profile configuration

#### check_settings()

Check atmosphere settings validity.

#### collected_base_profiles *: [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[VariableDefinition](modular_dales.vars.md#modular_dales.vars.VariableDefinition), [AtmosphericProfile](modular_dales.Atmosphere.atmosphere.md#modular_dales.Atmosphere.atmosphere.AtmosphericProfile) | [InterpolatedProfile](modular_dales.Atmosphere.atmosphere.md#modular_dales.Atmosphere.atmosphere.InterpolatedProfile)]*

#### collected_timed_profiles_by_time *: [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[float](https://docs.python.org/3/library/functions.html#float), [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[VariableDefinition](modular_dales.vars.md#modular_dales.vars.VariableDefinition), [AtmosphericProfile](modular_dales.Atmosphere.atmosphere.md#modular_dales.Atmosphere.atmosphere.AtmosphericProfile) | [InterpolatedProfile](modular_dales.Atmosphere.atmosphere.md#modular_dales.Atmosphere.atmosphere.InterpolatedProfile)]]*

#### interpolated_profiles *: [List](https://docs.python.org/3/library/typing.html#typing.List)[[InterpolatedProfile](modular_dales.Atmosphere.atmosphere.md#modular_dales.Atmosphere.atmosphere.InterpolatedProfile)]*

#### prepare_calculation()

Prepare atmosphere profiles and update internal state.

#### profile_writer *: AtmosphereProfileWriter*

#### shaped_profiles *: [List](https://docs.python.org/3/library/typing.html#typing.List)[[AtmosphericProfile](modular_dales.Atmosphere.atmosphere.md#modular_dales.Atmosphere.atmosphere.AtmosphericProfile)]*

#### sim *: simulation_module | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### timed_profiles *: [List](https://docs.python.org/3/library/typing.html#typing.List)[[TimedAtmosphereProfile](modular_dales.Atmosphere.atmosphere.md#modular_dales.Atmosphere.atmosphere.TimedAtmosphereProfile)]*

#### variables *: [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[VariableDefinition](modular_dales.vars.md#modular_dales.vars.VariableDefinition), [AtmosphereVariable](modular_dales.Atmosphere.atmosphere.md#modular_dales.Atmosphere.atmosphere.AtmosphereVariable)]*

#### write_files()

Write atmosphere profile files and generate plots.

### *class* modular_dales.Atmosphere.AtmosphereVariable(definition: [VariableDefinition](modular_dales.vars.md#modular_dales.vars.VariableDefinition), is_time_dependent: [bool](https://docs.python.org/3/library/functions.html#bool) = False, values: ndarray | [None](https://docs.python.org/3/library/constants.html#None) = None)

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Represents an atmospheric variable with its definition and associated values.

#### definition

The variable definition containing metadata and specifications.

* **Type:**
  [modular_dales.vars.VariableDefinition](modular_dales.vars.md#modular_dales.vars.VariableDefinition)

#### is_time_dependent

Whether the variable varies over time. Defaults to False.

* **Type:**
  [bool](https://docs.python.org/3/library/functions.html#bool)

#### values

Optional numpy array containing the variable’s values.

* **Type:**
  numpy.ndarray | None

#### definition *: [VariableDefinition](modular_dales.vars.md#modular_dales.vars.VariableDefinition)*

#### is_time_dependent *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

#### values *: ndarray | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

### *class* modular_dales.Atmosphere.AtmosphericProfile(variable: [VariableDefinition](modular_dales.vars.md#modular_dales.vars.VariableDefinition), shape: [str](https://docs.python.org/3/library/stdtypes.html#str), params: [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), [float](https://docs.python.org/3/library/functions.html#float) | TimeDependentScalar])

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Single atmospheric profile configuration.

#### evaluate(z: ndarray) → ndarray

#### params *: [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), [float](https://docs.python.org/3/library/functions.html#float) | TimeDependentScalar]*

Shape-specific parameters

#### shape *: [str](https://docs.python.org/3/library/stdtypes.html#str)*

Shape function (lin, exp, linmlsurf, expsinw, etc.)

#### timed_forcing_series(z: ndarray) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[float](https://docs.python.org/3/library/functions.html#float), ndarray] | [None](https://docs.python.org/3/library/constants.html#None)

#### variable *: [VariableDefinition](modular_dales.vars.md#modular_dales.vars.VariableDefinition)*

Variable definition (see modular_dales.Atmosphere.vars).

### *class* modular_dales.Atmosphere.FromLS2D

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Marker class to enable LS2D-driven soil/roughness in LSM.
Also enables injection of LS2D time series into atmosphere.

When added to `LSMModule` or `TimeDependentModule` (`lsm += FromLS2D()`), the
module will, if an [`LS2DAtmosphereModule`](#modular_dales.Atmosphere.LS2DAtmosphereModule) is present,
override soil temperature/moisture and soil type index from LS2D
and also set the bulk roughness lengths `z0mav` and `z0hav`
from LS2D time-mean values.

### *class* modular_dales.Atmosphere.InterpolatedProfile(variable: [VariableDefinition](modular_dales.vars.md#modular_dales.vars.VariableDefinition), z: [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)], points: [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float) | TimeDependentScalar], fill_value: [float](https://docs.python.org/3/library/functions.html#float) | [str](https://docs.python.org/3/library/stdtypes.html#str) = 'extrapolate')

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Single atmospheric profile configuration.

#### evaluate(z: ndarray) → ndarray

#### fill_value *: [float](https://docs.python.org/3/library/functions.html#float) | [str](https://docs.python.org/3/library/stdtypes.html#str)* *= 'extrapolate'*

Value to use for beyond the provided z range, or ‘extrapolate’ to extrapolate.

#### points *: [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float) | TimeDependentScalar]*

Values at the specified z coordinates

#### timed_forcing_series(z: ndarray) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[float](https://docs.python.org/3/library/functions.html#float), ndarray] | [None](https://docs.python.org/3/library/constants.html#None)

#### variable *: [VariableDefinition](modular_dales.vars.md#modular_dales.vars.VariableDefinition)*

Variable definition (see modular_dales.Atmosphere.vars).

#### z *: [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)]*

list of z coordinates for interpolation

### *class* modular_dales.Atmosphere.LS2DAtmosphereModule(sim: simulation_module | [None](https://docs.python.org/3/library/constants.html#None) = None, central_lat: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, central_lon: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, area_size: [float](https://docs.python.org/3/library/functions.html#float) = 1.0, case_name: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'ls2d_case', era5_path: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None, start_date: [datetime](https://docs.python.org/3/library/datetime.html#datetime.datetime) | [None](https://docs.python.org/3/library/constants.html#None) = None, end_date: [datetime](https://docs.python.org/3/library/datetime.html#datetime.datetime) | [None](https://docs.python.org/3/library/constants.html#None) = None, write_log: [bool](https://docs.python.org/3/library/functions.html#bool) = True, data_source: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'CDS', n_av: [int](https://docs.python.org/3/library/functions.html#int) = 0, method: [str](https://docs.python.org/3/library/stdtypes.html#str) = '2nd', init_tke: [float](https://docs.python.org/3/library/functions.html#float) = 0.1)

Bases: `simulation_module`

Atmosphere forcing module using LS2D ERA5 processing.

* Uses LS2D to download and process ERA5 data.
* Calls `era.calculate_forcings` in LS2D.
* Interpolates the resulting forcings onto the DALES
  vertical grid given by `GridDales` (`grid.zt`).
  > * Provides large-scale forcings and base profiles that can be
  >   : injected into [`AtmosphereModule`](#modular_dales.Atmosphere.AtmosphereModule) and surface modules.

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
[`AtmosphereModule`](#modular_dales.Atmosphere.AtmosphereModule) (when present), which then writes the
init file. The radiation background is
written to `backrad.inp.<exp_id>.nc`.

#### write_log *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= True*

### *class* modular_dales.Atmosphere.TimedAtmosphereProfile(time: [float](https://docs.python.org/3/library/functions.html#float), profile: [AtmosphericProfile](modular_dales.Atmosphere.atmosphere.md#modular_dales.Atmosphere.atmosphere.AtmosphericProfile) | [InterpolatedProfile](modular_dales.Atmosphere.atmosphere.md#modular_dales.Atmosphere.atmosphere.InterpolatedProfile))

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Profile assignment at a specific simulation time (seconds).

#### profile *: [AtmosphericProfile](modular_dales.Atmosphere.atmosphere.md#modular_dales.Atmosphere.atmosphere.AtmosphericProfile) | [InterpolatedProfile](modular_dales.Atmosphere.atmosphere.md#modular_dales.Atmosphere.atmosphere.InterpolatedProfile)*

#### time *: [float](https://docs.python.org/3/library/functions.html#float)*

### modular_dales.Atmosphere.exp(z, surf_val, lambda_val, z_ml=500)

### modular_dales.Atmosphere.expsinw(z, surf_val, H, amp, Hp)

### modular_dales.Atmosphere.lin(z, surf_val, ddz)

### modular_dales.Atmosphere.linmlsurf(z, lapse_rate, surf_val, offset_val=1.25, z_ml=500)

### Modules

| [`atmosphere`](modular_dales.Atmosphere.atmosphere.md#module-modular_dales.Atmosphere.atmosphere)                |    |
|------------------------------------------------------------------------------------------------------------------|----|
| [`ls2d_atmosphere`](modular_dales.Atmosphere.ls2d_atmosphere.md#module-modular_dales.Atmosphere.ls2d_atmosphere) |    |
| [`shapes`](modular_dales.Atmosphere.shapes.md#module-modular_dales.Atmosphere.shapes)                            |    |
