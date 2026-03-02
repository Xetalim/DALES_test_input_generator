# modular_dales.LBC.openbc

Open Boundary Conditions module for lateral boundary forcing.

### Classes

| [`Nest_in_AtmosphereProfiles`](#modular_dales.LBC.openbc.Nest_in_AtmosphereProfiles)(variable_mapping, ...)   | Nest DALES in a horizontally homogeneous atmosphere from profiles.   |
|---------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------|
| [`Nest_in_Dales`](#modular_dales.LBC.openbc.Nest_in_Dales)([outpath_coarse, ...])                             | Nest DALES inside another DALES simulation                           |
| [`Nest_in_Harmonie`](#modular_dales.LBC.openbc.Nest_in_Harmonie)([ml_glob, sfc_glob])                         | Nest a DALES simulation inside a HARMONIE model environment.         |
| [`do_openboundary`](#modular_dales.LBC.openbc.do_openboundary)(sim, nest_in_harmonie, ...)                    | Namelist parameters                                                  |

### *class* modular_dales.LBC.openbc.Nest_in_AtmosphereProfiles(variable_mapping: Dict[str, str]=<factory>, noise_std: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, noise_seed: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, noise_boundaries: [List](https://docs.python.org/3/library/typing.html#typing.List)[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [None](https://docs.python.org/3/library/constants.html#None) = None, noise_variables: [List](https://docs.python.org/3/library/typing.html#typing.List)[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [None](https://docs.python.org/3/library/constants.html#None) = None, atmosphere_module_name: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None, atmosphere_module: [AtmosphereModule](modular_dales.Atmosphere.atmosphere.md#modular_dales.Atmosphere.atmosphere.AtmosphereModule) | [None](https://docs.python.org/3/library/constants.html#None) = None)

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Nest DALES in a horizontally homogeneous atmosphere from profiles.

This mode does not require a separate coarse DALES or HARMONIE run.
Instead, it reuses the vertical profiles configured in the
`AtmosphereModule` of the *same* simulation to construct
open boundary fields that are uniform in x and y.

By default, the following mapping is used between open-boundary
variables and atmospheric profile variables (see `vars.ATMO_VARS_BY_NAME`):

> - u   <- ug
> - v   <- vg
> - w   <- w
> - thl <- thetal
> - qt  <- qt
> - e12 <- tke

You can override this mapping via `variable_mapping` if needed.

> This mode does *not* use any existing `AtmosphereModule` that is
> part of the `dales_simulation` modules. Instead, you must
> provide an explicit `AtmosphereModule` instance via
> `atmosphere_module`. This instance does not need to be added to
> the simulation’s module list, but it must have a valid `sim` with
> a grid attached.

#### atmosphere_module *: [AtmosphereModule](modular_dales.Atmosphere.atmosphere.md#modular_dales.Atmosphere.atmosphere.AtmosphereModule) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### atmosphere_module_name *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### noise_boundaries *: [List](https://docs.python.org/3/library/typing.html#typing.List)[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### noise_seed *: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### noise_std *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### noise_variables *: [List](https://docs.python.org/3/library/typing.html#typing.List)[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### variable_mapping *: [Dict](https://docs.python.org/3/library/typing.html#typing.Dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), [str](https://docs.python.org/3/library/stdtypes.html#str)]*

### *class* modular_dales.LBC.openbc.Nest_in_Dales(outpath_coarse: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None, outpath_coarse_old: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None, inpath_coarse: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None, inpath: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None)

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Nest DALES inside another DALES simulation

* **Parameters:**
  * **outpath_coarse** – Optional[str]
    Path to the output directory of the coarse parent DALES simulation
    used to source boundary fields for all time steps except t=0.
  * **outpath_coarse_old** – Optional[str]
    Path to the output directory of the previous DALES simulation,
    used to source initial boundary fields at t=0
    from the last time step in the output of the previous simulation
  * **inpath_coarse** – Optional[str]
    Path to the input directory of the previous DALES simulation,
    used to source initial boundary fields at t=0
    from the initfields.nc input of the previous simulation
  * **inpath** – Optional[str]
    Path to the input directory of the current DALES simulation,
    used to source initial fields for the current simulation

#### inpath *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### inpath_coarse *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### outpath_coarse *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### outpath_coarse_old *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

### *class* modular_dales.LBC.openbc.Nest_in_Harmonie(ml_glob: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None, sfc_glob: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None)

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Nest a DALES simulation inside a HARMONIE model environment.

This class configures and manages the coupling between a DALES
simulation and an enclosing HARMONIE model, handling the exchange
of necessary meteorological and surface fields.

* **Parameters:**
  * **ml_glob** – Optional[str]
    Glob pattern or path template for HARMONIE model-level fields
    (e.g., 3D atmospheric state) used to drive or constrain DALES.
  * **sfc_glob** – Optional[str]
    Glob pattern or path template for HARMONIE surface fields
    (e.g., surface fluxes, skin temperature) required for the nesting.

#### ml_glob *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### sfc_glob *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

### *class* modular_dales.LBC.openbc.do_openboundary(sim: ForwardRef('simulation_module') | [None](https://docs.python.org/3/library/constants.html#None) = None, nest_in_harmonie: [modular_dales.LBC.openbc.Nest_in_Harmonie](#modular_dales.LBC.openbc.Nest_in_Harmonie) | [None](https://docs.python.org/3/library/constants.html#None) = None, nest_in_dales: [modular_dales.LBC.openbc.Nest_in_Dales](#modular_dales.LBC.openbc.Nest_in_Dales) | [None](https://docs.python.org/3/library/constants.html#None) = None, nest_in_atmosphere: [modular_dales.LBC.openbc.Nest_in_AtmosphereProfiles](#modular_dales.LBC.openbc.Nest_in_AtmosphereProfiles) | [None](https://docs.python.org/3/library/constants.html#None) = None, e12: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, tracernames: List[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [None](https://docs.python.org/3/library/constants.html#None) = <factory>, tchunk: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, start: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None, time0: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None, end: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None, lopenbc: [bool](https://docs.python.org/3/library/functions.html#bool) = True, linithetero: [bool](https://docs.python.org/3/library/functions.html#bool) = True, lper: [bool](https://docs.python.org/3/library/functions.html#bool) = False, lbuoytop: [bool](https://docs.python.org/3/library/functions.html#bool) = True, dxint: [int](https://docs.python.org/3/library/functions.html#int) = 100, dyint: [int](https://docs.python.org/3/library/functions.html#int) = 100, dzint: [int](https://docs.python.org/3/library/functions.html#int) = -1, dxturb: [int](https://docs.python.org/3/library/functions.html#int) = -1, dyturb: [int](https://docs.python.org/3/library/functions.html#int) = -1, taum: [int](https://docs.python.org/3/library/functions.html#int) = 0, tauh: [int](https://docs.python.org/3/library/functions.html#int) = 20, pbc: [int](https://docs.python.org/3/library/functions.html#int) = 3, lsynturb: [bool](https://docs.python.org/3/library/functions.html#bool) = False, iturb: [int](https://docs.python.org/3/library/functions.html#int) = 0, tau: [int](https://docs.python.org/3/library/functions.html#int) = 180, lambda_: [int](https://docs.python.org/3/library/functions.html#int) = 1, nmodes: [int](https://docs.python.org/3/library/functions.html#int) = 100, lambdas: [int](https://docs.python.org/3/library/functions.html#int) = 1, lambdas_x: [int](https://docs.python.org/3/library/functions.html#int) = 1, lambdas_y: [int](https://docs.python.org/3/library/functions.html#int) = 1, lambdas_z: [int](https://docs.python.org/3/library/functions.html#int) = 1, igrw_damp: [int](https://docs.python.org/3/library/functions.html#int) = 0, lconstexner: [bool](https://docs.python.org/3/library/functions.html#bool) = True)

Bases: `simulation_module`

## Namelist parameters

*Section* `OPENBC`:
: - `dxint` (field `dxint`)
  - `dxturb` (field `dxturb`)
  - `dyint` (field `dyint`)
  - `dyturb` (field `dyturb`)
  - `dzint` (field `dzint`)
  - `iturb` (field `iturb`)
  - `lambda` (field `lambda_`)
  - `lambdas` (field `lambdas`)
  - `lambdas_x` (field `lambdas_x`)
  - `lambdas_y` (field `lambdas_y`)
  - `lambdas_z` (field `lambdas_z`)
  - `lbuoytop` (field `lbuoytop`)
  - `linithetero` (field `linithetero`)
  - `lopenbc` (field `lopenbc`)
  - `lper` (field `lper`)
  - `lsynturb` (field `lsynturb`)
  - `nmodes` (field `nmodes`)
  - `pbc` (field `pbc`)
  - `tau` (field `tau`)
  - `tauh` (field `tauh`)
  - `taum` (field `taum`)

*Section* `PHYSICS`:
: - `igrw_damp` (field `igrw_damp`)

*Section* `thermodynamics`:
: - `lconstexner` (field `lconstexner`)

#### boundaries *: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### dxint *: [int](https://docs.python.org/3/library/functions.html#int)* *= 100*

#### dxturb *: [int](https://docs.python.org/3/library/functions.html#int)* *= -1*

#### dyint *: [int](https://docs.python.org/3/library/functions.html#int)* *= 100*

#### dyturb *: [int](https://docs.python.org/3/library/functions.html#int)* *= -1*

#### dzint *: [int](https://docs.python.org/3/library/functions.html#int)* *= -1*

#### e12 *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### end *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### harmonieprepper *: harmoniePrepper | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### igrw_damp *: [int](https://docs.python.org/3/library/functions.html#int)* *= 0*

#### indices *: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### initfields *: [Any](https://docs.python.org/3/library/typing.html#typing.Any) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### iturb *: [int](https://docs.python.org/3/library/functions.html#int)* *= 0*

#### lambda_ *: [int](https://docs.python.org/3/library/functions.html#int)* *= 1*

#### lambdas *: [int](https://docs.python.org/3/library/functions.html#int)* *= 1*

#### lambdas_x *: [int](https://docs.python.org/3/library/functions.html#int)* *= 1*

#### lambdas_y *: [int](https://docs.python.org/3/library/functions.html#int)* *= 1*

#### lambdas_z *: [int](https://docs.python.org/3/library/functions.html#int)* *= 1*

#### lbuoytop *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= True*

#### lconstexner *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= True*

#### linithetero *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= True*

#### lopenbc *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= True*

#### lper *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

#### lsynturb *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= False*

#### nest_in_atmosphere *: [Nest_in_AtmosphereProfiles](#modular_dales.LBC.openbc.Nest_in_AtmosphereProfiles) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### nest_in_dales *: [Nest_in_Dales](#modular_dales.LBC.openbc.Nest_in_Dales) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### nest_in_harmonie *: [Nest_in_Harmonie](#modular_dales.LBC.openbc.Nest_in_Harmonie) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### nmodes *: [int](https://docs.python.org/3/library/functions.html#int)* *= 100*

#### openBCgrid *: GridDalesOpenBC* *= None*

#### pbc *: [int](https://docs.python.org/3/library/functions.html#int)* *= 3*

#### prepare_calculation()

Prepare data and setup for calculations.

#### sim *: simulation_module | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### start *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### tau *: [int](https://docs.python.org/3/library/functions.html#int)* *= 180*

#### tauh *: [int](https://docs.python.org/3/library/functions.html#int)* *= 20*

#### taum *: [int](https://docs.python.org/3/library/functions.html#int)* *= 0*

#### tchunk *: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### time0 *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### tracernames *: [List](https://docs.python.org/3/library/typing.html#typing.List)[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [None](https://docs.python.org/3/library/constants.html#None)*

#### write_files()

Write output files for this module.
