# modular_dales.Emission.create_emis

### Functions

| [`get_emis_sim_hours`](#modular_dales.Emission.create_emis.get_emis_sim_hours)(nml)       |    |
|-------------------------------------------------------------------------------------------|----|
| [`parse_tracer_type`](#modular_dales.Emission.create_emis.parse_tracer_type)(tracer_conf) |    |

### Classes

| [`PointSource`](#modular_dales.Emission.create_emis.PointSource)(x_idx, y_idx, height, ...)       |    |
|---------------------------------------------------------------------------------------------------|----|
| [`Tracer`](#modular_dales.Emission.create_emis.Tracer)(info, profile, pointsources)               |    |
| [`TracerInfo`](#modular_dales.Emission.create_emis.TracerInfo)(name, long_name, unit, molar_mass) |    |
| [`emissions`](#modular_dales.Emission.create_emis.emissions)(output_path, grid)                   |    |

### *class* modular_dales.Emission.create_emis.PointSource(x_idx: [int](https://docs.python.org/3/library/functions.html#int), y_idx: [int](https://docs.python.org/3/library/functions.html#int), height: [float](https://docs.python.org/3/library/functions.html#float), temperature: [float](https://docs.python.org/3/library/functions.html#float), volume: [float](https://docs.python.org/3/library/functions.html#float), emission: [float](https://docs.python.org/3/library/functions.html#float), stack_exit_area: [float](https://docs.python.org/3/library/functions.html#float))

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

#### emission *: [float](https://docs.python.org/3/library/functions.html#float)*

#### height *: [float](https://docs.python.org/3/library/functions.html#float)*

#### stack_exit_area *: [float](https://docs.python.org/3/library/functions.html#float)*

#### temperature *: [float](https://docs.python.org/3/library/functions.html#float)*

#### volume *: [float](https://docs.python.org/3/library/functions.html#float)*

#### x_idx *: [int](https://docs.python.org/3/library/functions.html#int)*

#### y_idx *: [int](https://docs.python.org/3/library/functions.html#int)*

### *class* modular_dales.Emission.create_emis.Tracer(info: [modular_dales.Emission.create_emis.TracerInfo](#modular_dales.Emission.create_emis.TracerInfo), profile: [list](https://docs.python.org/3/library/stdtypes.html#list) | numpy.ndarray, pointsources: List[[modular_dales.Emission.create_emis.PointSource](#modular_dales.Emission.create_emis.PointSource)])

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

#### info *: [TracerInfo](#modular_dales.Emission.create_emis.TracerInfo)*

#### pointsources *: [List](https://docs.python.org/3/library/typing.html#typing.List)[[PointSource](#modular_dales.Emission.create_emis.PointSource)]*

#### profile *: [list](https://docs.python.org/3/library/stdtypes.html#list) | ndarray*

### *class* modular_dales.Emission.create_emis.TracerInfo(name: [str](https://docs.python.org/3/library/stdtypes.html#str), long_name: [str](https://docs.python.org/3/library/stdtypes.html#str), unit: [str](https://docs.python.org/3/library/stdtypes.html#str), molar_mass: [float](https://docs.python.org/3/library/functions.html#float), lemis: [bool](https://docs.python.org/3/library/functions.html#bool) = False, lreact: [bool](https://docs.python.org/3/library/functions.html#bool) = False, ldep: [bool](https://docs.python.org/3/library/functions.html#bool) = False, lags: [bool](https://docs.python.org/3/library/functions.html#bool) = False, profile: numpy.ndarray | [None](https://docs.python.org/3/library/constants.html#None) = None)

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

#### lags *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

#### ldep *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

#### lemis *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

#### long_name *: [str](https://docs.python.org/3/library/stdtypes.html#str)*

#### lreact *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

#### molar_mass *: [float](https://docs.python.org/3/library/functions.html#float)*

#### name *: [str](https://docs.python.org/3/library/stdtypes.html#str)*

#### profile *: ndarray | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### unit *: [str](https://docs.python.org/3/library/stdtypes.html#str)*

### *class* modular_dales.Emission.create_emis.emissions(output_path, grid: [GridDales](modular_dales.Geometry.GridDales.md#modular_dales.Geometry.GridDales))

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

#### add_pointsource(tracername, source_info: [PointSource](#modular_dales.Emission.create_emis.PointSource))

#### add_pre_existing_tracer(tracername)

#### add_tracer(info: [TracerInfo](#modular_dales.Emission.create_emis.TracerInfo))

#### emissionmap_to_netcdf(datetime_str, tracername)

#### pointsources_to_netcdf(datetime_str, tracername)

#### set_profile(tracername, profile)

#### set_tracer_property(tracername, property_name, property_value)

#### tracerinfo_to_netcdf(file, tracername)

#### write_hourly_emission_data(datetime_str)

#### write_tracers_file()

### modular_dales.Emission.create_emis.get_emis_sim_hours(nml)

### modular_dales.Emission.create_emis.parse_tracer_type(tracer_conf) → [TracerInfo](#modular_dales.Emission.create_emis.TracerInfo)
