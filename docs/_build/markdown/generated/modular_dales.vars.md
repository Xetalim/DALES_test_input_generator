# modular_dales.vars

### Classes

| [`VariableDefinition`](#modular_dales.vars.VariableDefinition)(name, long_name, unit[, ...])   | Definition of an atmospheric variable.   |
|------------------------------------------------------------------------------------------------|------------------------------------------|

### *class* modular_dales.vars.VariableDefinition(name: [str](https://docs.python.org/3/library/stdtypes.html#str), long_name: [str](https://docs.python.org/3/library/stdtypes.html#str), unit: [str](https://docs.python.org/3/library/stdtypes.html#str), can_be_time_dependent: [bool](https://docs.python.org/3/library/functions.html#bool) = False, must_only_be_time_dependent: [bool](https://docs.python.org/3/library/functions.html#bool) = False, is_profile: [bool](https://docs.python.org/3/library/functions.html#bool) = True, time_dependent_name: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None)

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Definition of an atmospheric variable.

This is the user-facing object you should import and pass to
AtmosphericProfile / InterpolatedProfile instead of raw strings.

#### can_be_time_dependent *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

#### is_profile *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= True*

#### long_name *: [str](https://docs.python.org/3/library/stdtypes.html#str)*

#### must_only_be_time_dependent *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

#### name *: [str](https://docs.python.org/3/library/stdtypes.html#str)*

#### time_dependent_name *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### unit *: [str](https://docs.python.org/3/library/stdtypes.html#str)*
