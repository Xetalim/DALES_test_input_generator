# modular_dales.IBM.IBM

Immersed Boundary Method (IBM) module for building and geometry representation.

### Classes

| [`FromAHN`](#modular_dales.IBM.IBM.FromAHN)()                                              | AHN (Actueel Hoogtebestand Nederland) approach for IBM.                                                |
|--------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------|
| [`IBMCreatorClass`](#modular_dales.IBM.IBM.IBMCreatorClass)(grid)                          | A class for creating and managing IBM (Immersed Boundary Method) modifications to DALES grid geometry. |
| [`IBMModification`](#modular_dales.IBM.IBM.IBMModification)(geometry, height, params, ...) | Single IBM modification.                                                                               |
| [`IBMModifications`](#modular_dales.IBM.IBM.IBMModifications)(modifications)               | Collection of IBM modifications.                                                                       |
| [`IBMModule`](#modular_dales.IBM.IBM.IBMModule)(sim, ibm_modifications, from_ahn, ...)     | IBM Module for DALES simulation.                                                                       |

### *class* modular_dales.IBM.IBM.FromAHN

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

AHN (Actueel Hoogtebestand Nederland) approach for IBM.
Add this to IBMModule to generate IBM geometry from AHN data.

### *class* modular_dales.IBM.IBM.IBMCreatorClass(grid: [GridDales](modular_dales.Geometry.GridDales.md#modular_dales.Geometry.GridDales))

Bases: `ModifierClass`

A class for creating and managing IBM (Immersed Boundary Method) modifications to DALES grid geometry.

This class handles the creation and output of boundary condition heights for immersed boundary
method simulations. It maintains a 2D array of boundary condition heights that can be modified
based on geometry specifications and exported to NetCDF format.

#### bc_height

A 2D array storing boundary condition heights, initialized with zeros
matching the grid mesh dimensions.

* **Type:**
  np.ndarray

#### do_modification(geometry, modification)

#### output_nc(filename)

### *class* modular_dales.IBM.IBM.IBMModification(geometry: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None, height: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, params: Dict[str, ~typing.Any]=<factory>)

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Single IBM modification. Set Geometry type and parameters to define the modification.

#### geometry *: [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

Geometry type

#### height *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

Building height

#### params *: [Dict](https://docs.python.org/3/library/typing.html#typing.Dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]*

Geometry-specific parameters

### *class* modular_dales.IBM.IBM.IBMModifications(modifications: [List](https://docs.python.org/3/library/typing.html#typing.List)[[IBMModification](#modular_dales.IBM.IBM.IBMModification)] = <factory>)

Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Collection of IBM modifications.

#### apply_to_config(config: [Dict](https://docs.python.org/3/library/typing.html#typing.Dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), [Any](https://docs.python.org/3/library/typing.html#typing.Any)]) → [None](https://docs.python.org/3/library/constants.html#None)

Apply modifications to config dictionary.

#### modifications *: [List](https://docs.python.org/3/library/typing.html#typing.List)[[IBMModification](#modular_dales.IBM.IBM.IBMModification)]*

### *class* modular_dales.IBM.IBM.IBMModule(sim: simulation_module | [None](https://docs.python.org/3/library/constants.html#None) = None, ibm_modifications: [List](https://docs.python.org/3/library/typing.html#typing.List)[[IBMModification](#modular_dales.IBM.IBM.IBMModification)] = <factory>, from_ahn: [FromAHN](#modular_dales.IBM.IBM.FromAHN) | [None](https://docs.python.org/3/library/constants.html#None) = None, ahn_zeroes_buffer: [int](https://docs.python.org/3/library/functions.html#int) = 5, subtract_ahn_mode: [bool](https://docs.python.org/3/library/functions.html#bool) = True, apply_ibm: [bool](https://docs.python.org/3/library/functions.html#bool) = True, ibas_prf: [int](https://docs.python.org/3/library/functions.html#int) = 2, iadv_mom: [int](https://docs.python.org/3/library/functions.html#int) = 2)

Bases: `simulation_module`

IBM Module for DALES simulation.

This module handles Immersed Boundary Method (IBM) configuration and setup for DALES
(Dutch Atmospheric Large-Eddy Simulation) simulations. It manages IBM modifications,
AHN (Actueel Hoogtebestand Nederland) terrain data integration, and generates IBM
geometry files.

#### sim

Reference to parent simulation module.

* **Type:**
  Optional[simulation_module]

#### ibm_modifications

List of IBM geometry modifications to apply.

* **Type:**
  List[[IBMModification](#modular_dales.IBM.IBM.IBMModification)]

#### ibm_generator

IBM geometry creator instance, initialized during preparation.

* **Type:**
  Any

#### from_ahn

AHN terrain data configuration object.

* **Type:**
  Optional[[FromAHN](#modular_dales.IBM.IBM.FromAHN)]

#### ahn_zeroes_buffer

Buffer size for AHN zero elevation handling. Default: 5.

* **Type:**
  [int](https://docs.python.org/3/library/functions.html#int)

#### subtract_ahn_mode

Whether to subtract AHN elevation data. Default: True.

* **Type:**
  [bool](https://docs.python.org/3/library/functions.html#bool)

#### apply_ibm

Enable/disable IBM application. Default: True.

* **Type:**
  [bool](https://docs.python.org/3/library/functions.html#bool)

#### ibas_prf

PRF advection scheme Default: 2.

* **Type:**
  [int](https://docs.python.org/3/library/functions.html#int)

#### iadv_mom

Advection scheme for momentum. Default: 2.

* **Type:**
  [int](https://docs.python.org/3/library/functions.html#int)

#### \_\_post_init_\_()

Initialize module and set module name.

#### \_\_add_\_()

Add IBM modification or AHN configuration to module.

#### \_\_iadd_\_()

In-place addition operator for fluent API.

#### prepare_calculation()

Set up IBM geometry, integrate AHN data if configured, and apply modifications.

#### check_settings()

–

#### write_files()

Export IBM geometry to NetCDF file.

## Namelist parameters

*Section* `DYNAMICS`:
: - `IADV_MOM` (field `iadv_mom`) (required)
  - `IBAS_PRF` (field `ibas_prf`) (required)

*Section* `NAMIBM`:
: - `lapply_ibm` (field `apply_ibm`) (required)

#### ahn_zeroes_buffer *: [int](https://docs.python.org/3/library/functions.html#int)* *= 5*

#### apply_ibm *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= True*

#### check_settings()

Check IBM settings validity.

#### from_ahn *: [FromAHN](#modular_dales.IBM.IBM.FromAHN) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### iadv_mom *: [int](https://docs.python.org/3/library/functions.html#int)* *= 2*

#### ibas_prf *: [int](https://docs.python.org/3/library/functions.html#int)* *= 2*

#### ibm_generator *: [Any](https://docs.python.org/3/library/typing.html#typing.Any)* *= None*

#### ibm_modifications *: [List](https://docs.python.org/3/library/typing.html#typing.List)[[IBMModification](#modular_dales.IBM.IBM.IBMModification)]*

#### prepare_calculation()

Set up IBM geometry and update namelist.

#### sim *: simulation_module | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### subtract_ahn_mode *: [bool](https://docs.python.org/3/library/functions.html#bool)* *= True*

#### write_files()

Write IBM files.
