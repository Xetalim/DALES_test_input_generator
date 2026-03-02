# modular_dales.LBC.Nesting_Topology

### Classes

| [`NestingTopology`](#modular_dales.LBC.Nesting_Topology.NestingTopology)(sim, nestings, my_idx)   | Nesting topology module to determine grid relationships for nested simulations.   |
|---------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------|

### *class* modular_dales.LBC.Nesting_Topology.NestingTopology(sim: simulation_module | [None](https://docs.python.org/3/library/constants.html#None) = None, nestings: [List](https://docs.python.org/3/library/typing.html#typing.List)[[GridDales](modular_dales.Geometry.GridDales.md#modular_dales.Geometry.GridDales)] = <factory>, my_idx: [int](https://docs.python.org/3/library/functions.html#int) | [None](https://docs.python.org/3/library/constants.html#None) = None)

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
