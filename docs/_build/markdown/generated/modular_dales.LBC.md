# modular_dales.LBC

Public API for lateral boundary condition (LBC) helpers.

### *class* modular_dales.LBC.Nest_in_AtmosphereProfiles(variable_mapping: Dict[str, str]=<factory>, noise_std: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, noise_seed: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, noise_boundaries: [List](https://docs.python.org/3/library/typing.html#typing.List)[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [None](https://docs.python.org/3/library/constants.html#None) = None, noise_variables: [List](https://docs.python.org/3/library/typing.html#typing.List)[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [None](https://docs.python.org/3/library/constants.html#None) = None, atmosphere_module_name: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None, atmosphere_module: [AtmosphereModule](modular_dales.Atmosphere.atmosphere.md#modular_dales.Atmosphere.atmosphere.AtmosphereModule) | [None](https://docs.python.org/3/library/constants.html#None) = None)

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

### *class* modular_dales.LBC.Nest_in_Dales(outpath_coarse: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None, outpath_coarse_old: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None, inpath_coarse: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None, inpath: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None)

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

### *class* modular_dales.LBC.NestingIndices(ix_west: [int](https://docs.python.org/3/library/functions.html#int), ix_east: [int](https://docs.python.org/3/library/functions.html#int), iy_south: [int](https://docs.python.org/3/library/functions.html#int), iy_north: [int](https://docs.python.org/3/library/functions.html#int), supergrid_x0: [float](https://docs.python.org/3/library/functions.html#float), supergrid_y0: [float](https://docs.python.org/3/library/functions.html#float), subgrid_x0: [float](https://docs.python.org/3/library/functions.html#float), subgrid_y0: [float](https://docs.python.org/3/library/functions.html#float), supergrid: ForwardRef('GridDales') | [None](https://docs.python.org/3/library/constants.html#None) = None)

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

#### ix_east *: [int](https://docs.python.org/3/library/functions.html#int)*

#### ix_west *: [int](https://docs.python.org/3/library/functions.html#int)*

#### iy_north *: [int](https://docs.python.org/3/library/functions.html#int)*

#### iy_south *: [int](https://docs.python.org/3/library/functions.html#int)*

#### subgrid_x0 *: [float](https://docs.python.org/3/library/functions.html#float)*

#### subgrid_y0 *: [float](https://docs.python.org/3/library/functions.html#float)*

#### supergrid *: [GridDales](modular_dales.Geometry.GridDales.md#modular_dales.Geometry.GridDales) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### supergrid_x0 *: [float](https://docs.python.org/3/library/functions.html#float)*

#### supergrid_y0 *: [float](https://docs.python.org/3/library/functions.html#float)*

### *class* modular_dales.LBC.NestingTopology(sim: simulation_module | [None](https://docs.python.org/3/library/constants.html#None) = None, nestings: [List](https://docs.python.org/3/library/typing.html#typing.List)[[GridDales](modular_dales.Geometry.GridDales.md#modular_dales.Geometry.GridDales)] = <factory>, my_idx: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None)

Bases: `simulation_module`

Nesting topology module to determine grid relationships for nested simulations.

* **Parameters:**
  **sim** – Parent simulation instance

## Namelist parameters

*Section* `namcrosssection`:
: - `lcross` (field `lcross`) (required)

#### indices *: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### lcross *: [bool](https://docs.python.org/3/library/functions.html#bool) | [None](https://docs.python.org/3/library/constants.html#None)* *= True*

#### my_idx *: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### nestings *: [List](https://docs.python.org/3/library/typing.html#typing.List)[[GridDales](modular_dales.Geometry.GridDales.md#modular_dales.Geometry.GridDales)]*

#### openbc_module *: [do_openboundary](modular_dales.LBC.openbc.md#modular_dales.LBC.openbc.do_openboundary) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### prepare_calculation()

Return nesting info for a specific index if available.

#### sim *: simulation_module | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### subnest_subgrid *: GridDalesOpenBC | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### subnest_supergrid *: GridDalesOpenBC | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### supernest_subgrid *: GridDalesOpenBC | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### supernest_supergrid *: [GridDales](modular_dales.Geometry.GridDales.md#modular_dales.Geometry.GridDales) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### write_files()

Write output files for this module.

### *class* modular_dales.LBC.do_openboundary(sim: ForwardRef('simulation_module') | [None](https://docs.python.org/3/library/constants.html#None) = None, nest_in_harmonie: [modular_dales.LBC.openbc.Nest_in_Harmonie](modular_dales.LBC.openbc.md#modular_dales.LBC.openbc.Nest_in_Harmonie) | [None](https://docs.python.org/3/library/constants.html#None) = None, nest_in_dales: [modular_dales.LBC.openbc.Nest_in_Dales](modular_dales.LBC.openbc.md#modular_dales.LBC.openbc.Nest_in_Dales) | [None](https://docs.python.org/3/library/constants.html#None) = None, nest_in_atmosphere: [modular_dales.LBC.openbc.Nest_in_AtmosphereProfiles](modular_dales.LBC.openbc.md#modular_dales.LBC.openbc.Nest_in_AtmosphereProfiles) | [None](https://docs.python.org/3/library/constants.html#None) = None, e12: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, tracernames: List[[str](https://docs.python.org/3/library/stdtypes.html#str)] | [None](https://docs.python.org/3/library/constants.html#None) = <factory>, tchunk: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None, start: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None, time0: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None, end: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None, lopenbc: [bool](https://docs.python.org/3/library/functions.html#bool) = True, linithetero: [bool](https://docs.python.org/3/library/functions.html#bool) = True, lper: [bool](https://docs.python.org/3/library/functions.html#bool) = False, lbuoytop: [bool](https://docs.python.org/3/library/functions.html#bool) = True, dxint: [int](https://docs.python.org/3/library/functions.html#int) = 100, dyint: [int](https://docs.python.org/3/library/functions.html#int) = 100, dzint: [int](https://docs.python.org/3/library/functions.html#int) = -1, dxturb: [int](https://docs.python.org/3/library/functions.html#int) = -1, dyturb: [int](https://docs.python.org/3/library/functions.html#int) = -1, taum: [int](https://docs.python.org/3/library/functions.html#int) = 0, tauh: [int](https://docs.python.org/3/library/functions.html#int) = 20, pbc: [int](https://docs.python.org/3/library/functions.html#int) = 3, lsynturb: [bool](https://docs.python.org/3/library/functions.html#bool) = False, iturb: [int](https://docs.python.org/3/library/functions.html#int) = 0, tau: [int](https://docs.python.org/3/library/functions.html#int) = 180, lambda_: [int](https://docs.python.org/3/library/functions.html#int) = 1, nmodes: [int](https://docs.python.org/3/library/functions.html#int) = 100, lambdas: [int](https://docs.python.org/3/library/functions.html#int) = 1, lambdas_x: [int](https://docs.python.org/3/library/functions.html#int) = 1, lambdas_y: [int](https://docs.python.org/3/library/functions.html#int) = 1, lambdas_z: [int](https://docs.python.org/3/library/functions.html#int) = 1, igrw_damp: [int](https://docs.python.org/3/library/functions.html#int) = 0, lconstexner: [bool](https://docs.python.org/3/library/functions.html#bool) = True)

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

#### nest_in_atmosphere *: [Nest_in_AtmosphereProfiles](modular_dales.LBC.openbc.md#modular_dales.LBC.openbc.Nest_in_AtmosphereProfiles) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### nest_in_dales *: [Nest_in_Dales](modular_dales.LBC.openbc.md#modular_dales.LBC.openbc.Nest_in_Dales) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### nest_in_harmonie *: [Nest_in_Harmonie](modular_dales.LBC.openbc.md#modular_dales.LBC.openbc.Nest_in_Harmonie) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

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

### Modules

| [`Nesting_Topology`](modular_dales.LBC.Nesting_Topology.md#module-modular_dales.LBC.Nesting_Topology)   |                                                               |
|---------------------------------------------------------------------------------------------------------|---------------------------------------------------------------|
| [`boundary_info`](modular_dales.LBC.boundary_info.md#module-modular_dales.LBC.boundary_info)            |                                                               |
| [`nesting_idx`](modular_dales.LBC.nesting_idx.md#module-modular_dales.LBC.nesting_idx)                  |                                                               |
| [`openbc`](modular_dales.LBC.openbc.md#module-modular_dales.LBC.openbc)                                 | Open Boundary Conditions module for lateral boundary forcing. |
