# modular_dales.Atmosphere.atmosphere

### Functions

| [`evaluate_profile_map`](#modular_dales.Atmosphere.atmosphere.evaluate_profile_map)(z, profiles_by_var)   |    |
|-----------------------------------------------------------------------------------------------------------|----|

### Classes

| [`AtmosphereModule`](#modular_dales.Atmosphere.atmosphere.AtmosphereModule)(sim, shaped_profiles, ...)        | Atmosphere module for profile setup.                                          |
|---------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------|
| [`AtmosphereVariable`](#modular_dales.Atmosphere.atmosphere.AtmosphereVariable)(definition[, ...])            | Represents an atmospheric variable with its definition and associated values. |
| [`AtmosphericProfile`](#modular_dales.Atmosphere.atmosphere.AtmosphericProfile)(variable, shape, params)      | Single atmospheric profile configuration.                                     |
| [`InterpolatedProfile`](#modular_dales.Atmosphere.atmosphere.InterpolatedProfile)(variable, z, points[, ...]) | Single atmospheric profile configuration.                                     |
| [`TimedAtmosphereProfile`](#modular_dales.Atmosphere.atmosphere.TimedAtmosphereProfile)(time, profile)        | Profile assignment at a specific simulation time (seconds).                   |

### *class* modular_dales.Atmosphere.atmosphere.AtmosphereModule(sim: simulation_module | [None](https://docs.python.org/3/library/constants.html#None) = None, shaped_profiles: [List](https://docs.python.org/3/library/typing.html#typing.List)[[AtmosphericProfile](#modular_dales.Atmosphere.atmosphere.AtmosphericProfile)] = <factory>, interpolated_profiles: [List](https://docs.python.org/3/library/typing.html#typing.List)[[InterpolatedProfile](#modular_dales.Atmosphere.atmosphere.InterpolatedProfile)] = <factory>, timed_profiles: [List](https://docs.python.org/3/library/typing.html#typing.List)[[TimedAtmosphereProfile](#modular_dales.Atmosphere.atmosphere.TimedAtmosphereProfile)] = <factory>)

Bases: `simulation_module`

Atmosphere module for profile setup.

* **Parameters:**
  * **sim** – Parent simulation instance
  * **shaped_profiles** – Mapping of variable name to shaped profile configuration
  * **interpolated_profiles** – Mapping of variable name to interpolated profile configuration

#### check_settings()

Check atmosphere settings validity.

#### collected_base_profiles *: [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[VariableDefinition](modular_dales.vars.md#modular_dales.vars.VariableDefinition), [AtmosphericProfile](#modular_dales.Atmosphere.atmosphere.AtmosphericProfile) | [InterpolatedProfile](#modular_dales.Atmosphere.atmosphere.InterpolatedProfile)]*

#### collected_timed_profiles_by_time *: [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[float](https://docs.python.org/3/library/functions.html#float), [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[VariableDefinition](modular_dales.vars.md#modular_dales.vars.VariableDefinition), [AtmosphericProfile](#modular_dales.Atmosphere.atmosphere.AtmosphericProfile) | [InterpolatedProfile](#modular_dales.Atmosphere.atmosphere.InterpolatedProfile)]]*

#### interpolated_profiles *: [List](https://docs.python.org/3/library/typing.html#typing.List)[[InterpolatedProfile](#modular_dales.Atmosphere.atmosphere.InterpolatedProfile)]*

#### prepare_calculation()

Prepare atmosphere profiles and update internal state.

#### profile_writer *: AtmosphereProfileWriter*

#### shaped_profiles *: [List](https://docs.python.org/3/library/typing.html#typing.List)[[AtmosphericProfile](#modular_dales.Atmosphere.atmosphere.AtmosphericProfile)]*

#### sim *: simulation_module | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### timed_profiles *: [List](https://docs.python.org/3/library/typing.html#typing.List)[[TimedAtmosphereProfile](#modular_dales.Atmosphere.atmosphere.TimedAtmosphereProfile)]*

#### variables *: [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[VariableDefinition](modular_dales.vars.md#modular_dales.vars.VariableDefinition), [AtmosphereVariable](#modular_dales.Atmosphere.atmosphere.AtmosphereVariable)]*

#### write_files()

Write atmosphere profile files and generate plots.

### *class* modular_dales.Atmosphere.atmosphere.AtmosphereVariable(definition: [VariableDefinition](modular_dales.vars.md#modular_dales.vars.VariableDefinition), is_time_dependent: [bool](https://docs.python.org/3/library/functions.html#bool) = False, values: ndarray | [None](https://docs.python.org/3/library/constants.html#None) = None)

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

### *class* modular_dales.Atmosphere.atmosphere.AtmosphericProfile(variable: [VariableDefinition](modular_dales.vars.md#modular_dales.vars.VariableDefinition), shape: [str](https://docs.python.org/3/library/stdtypes.html#str), params: [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), [float](https://docs.python.org/3/library/functions.html#float) | TimeDependentScalar])

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

### *class* modular_dales.Atmosphere.atmosphere.InterpolatedProfile(variable: [VariableDefinition](modular_dales.vars.md#modular_dales.vars.VariableDefinition), z: [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)], points: [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float) | TimeDependentScalar], fill_value: [float](https://docs.python.org/3/library/functions.html#float) | [str](https://docs.python.org/3/library/stdtypes.html#str) = 'extrapolate')

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

### *class* modular_dales.Atmosphere.atmosphere.TimedAtmosphereProfile(time: [float](https://docs.python.org/3/library/functions.html#float), profile: [AtmosphericProfile](#modular_dales.Atmosphere.atmosphere.AtmosphericProfile) | [InterpolatedProfile](#modular_dales.Atmosphere.atmosphere.InterpolatedProfile))

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Profile assignment at a specific simulation time (seconds).

#### profile *: [AtmosphericProfile](#modular_dales.Atmosphere.atmosphere.AtmosphericProfile) | [InterpolatedProfile](#modular_dales.Atmosphere.atmosphere.InterpolatedProfile)*

#### time *: [float](https://docs.python.org/3/library/functions.html#float)*

### modular_dales.Atmosphere.atmosphere.evaluate_profile_map(z: ndarray, profiles_by_var: [Mapping](https://docs.python.org/3/library/typing.html#typing.Mapping)[[VariableDefinition](modular_dales.vars.md#modular_dales.vars.VariableDefinition), [AtmosphericProfile](#modular_dales.Atmosphere.atmosphere.AtmosphericProfile) | [InterpolatedProfile](#modular_dales.Atmosphere.atmosphere.InterpolatedProfile)])
