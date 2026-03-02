# modular_dales.Surface.surface

### Classes

| [`ConstantFluxesModule`](#modular_dales.Surface.surface.ConstantFluxesModule)([sim, wtsurf, wqsurf, ...])           | Constant heat and moisture fluxes without shear stress (isurf=4).   |
|---------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------|
| [`ConstantFluxesWithShearModule`](#modular_dales.Surface.surface.ConstantFluxesWithShearModule)([sim, wtsurf, ...]) | Constant heat and moisture fluxes with shear stress (isurf=3).      |
| [`ConstantSurfaceTemperatureModule`](#modular_dales.Surface.surface.ConstantSurfaceTemperatureModule)([sim, ...])   | Constant surface temperature boundary condition (isurf=2).          |
| [`SurfaceModule`](#modular_dales.Surface.surface.SurfaceModule)([sim])                                              | Surface simulation module base class.                               |

### *class* modular_dales.Surface.surface.ConstantFluxesModule(sim: simulation_module | [None](https://docs.python.org/3/library/constants.html#None) = None, wtsurf: [float](https://docs.python.org/3/library/functions.html#float) | TimeDependentScalar | [None](https://docs.python.org/3/library/constants.html#None) = None, wqsurf: [float](https://docs.python.org/3/library/functions.html#float) | TimeDependentScalar | [None](https://docs.python.org/3/library/constants.html#None) = None, ps: [float](https://docs.python.org/3/library/functions.html#float) | TimeDependentScalar | [None](https://docs.python.org/3/library/constants.html#None) = None, z0mav: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, z0hav: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, albedoav: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None)

Bases: [`SurfaceModule`](#modular_dales.Surface.surface.SurfaceModule)

Constant heat and moisture fluxes without shear stress (isurf=4).

When added to a simulation, automatically enables isurf=4 in the namelist.
Requires surface config to contain ‘wtsurf’, ‘wqsurf’, and ‘ps’ parameters.

* **Parameters:**
  * **sim** – Parent simulation instance
  * **isurf** – Surface scheme selector (4 for constant fluxes)
  * **wtsurf** – Heat flux (K m/s)
  * **wqsurf** – Moisture flux (kg/kg m/s)
  * **ps** – Surface pressure (Pa)
  * **z0mav** – Momentum roughness length (m)
  * **z0hav** – Heat roughness length (m)
  * **albedoav** – Surface albedo (dimensionless)

## Namelist parameters

*Section* `NAMSURFACE`:
: - `albedoav` (field `albedoav`)
  - `isurf` (field `isurf`)
  - `ps` (field `ps`) (required)
  - `wqsurf` (field `wqsurf`) (required)
  - `wtsurf` (field `wtsurf`) (required)
  - `z0hav` (field `z0hav`) (required)
  - `z0mav` (field `z0mav`) (required)

#### albedoav *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### check_settings()

Validate constant fluxes settings.

#### do_config()

Configure surface settings and namelist.

#### isurf *: [int](https://docs.python.org/3/library/functions.html#int)* *= 4*

#### prepare_calculation()

Prepare data and setup for calculations.

#### prepare_calculations()

No additional preparation needed.

#### ps *: [float](https://docs.python.org/3/library/functions.html#float) | TimeDependentScalar | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### sim *: simulation_module | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### wqsurf *: [float](https://docs.python.org/3/library/functions.html#float) | TimeDependentScalar | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### write_files()

No files to write.

#### wtsurf *: [float](https://docs.python.org/3/library/functions.html#float) | TimeDependentScalar | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### z0hav *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### z0mav *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

### *class* modular_dales.Surface.surface.ConstantFluxesWithShearModule(sim: simulation_module | [None](https://docs.python.org/3/library/constants.html#None) = None, wtsurf: [float](https://docs.python.org/3/library/functions.html#float) | TimeDependentScalar | [None](https://docs.python.org/3/library/constants.html#None) = None, wqsurf: [float](https://docs.python.org/3/library/functions.html#float) | TimeDependentScalar | [None](https://docs.python.org/3/library/constants.html#None) = None, ustin: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, ps: [float](https://docs.python.org/3/library/functions.html#float) | TimeDependentScalar | [None](https://docs.python.org/3/library/constants.html#None) = None, z0mav: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, z0hav: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, albedoav: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None)

Bases: [`SurfaceModule`](#modular_dales.Surface.surface.SurfaceModule)

Constant heat and moisture fluxes with shear stress (isurf=3).

When added to a simulation, automatically enables isurf=3 in the namelist.
Requires surface config to contain ‘wtsurf’, ‘wqsurf’, and ‘ustin’ parameters.
:param sim: Parent simulation instance
:param isurf: Surface scheme selector (3 for constant fluxes with shear)
:param wtsurf: Heat flux (K m/s)
:param wqsurf: Moisture flux (kg/kg m/s)
:param ustin: Friction velocity (m/s)
:param ps: Surface pressure (Pa)
:param z0mav: Momentum roughness length (m)
:param z0hav: Heat roughness length (m)
:param albedoav: Surface albedo (dimensionless)

## Namelist parameters

*Section* `NAMSURFACE`:
: - `albedoav` (field `albedoav`)
  - `isurf` (field `isurf`)
  - `ps` (field `ps`) (required)
  - `ustin` (field `ustin`) (required)
  - `wqsurf` (field `wqsurf`) (required)
  - `wtsurf` (field `wtsurf`) (required)
  - `z0hav` (field `z0hav`) (required)
  - `z0mav` (field `z0mav`) (required)

#### albedoav *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### check_settings()

Validate constant fluxes with shear settings.

#### do_config()

Configure surface settings and namelist.

#### isurf *: [int](https://docs.python.org/3/library/functions.html#int)* *= 3*

#### prepare_calculation()

Prepare data and setup for calculations.

#### prepare_calculations()

No additional preparation needed.

#### ps *: [float](https://docs.python.org/3/library/functions.html#float) | TimeDependentScalar | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### sim *: simulation_module | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### ustin *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### wqsurf *: [float](https://docs.python.org/3/library/functions.html#float) | TimeDependentScalar | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### write_files()

No files to write.

#### wtsurf *: [float](https://docs.python.org/3/library/functions.html#float) | TimeDependentScalar | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### z0hav *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### z0mav *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

### *class* modular_dales.Surface.surface.ConstantSurfaceTemperatureModule(sim: simulation_module | [None](https://docs.python.org/3/library/constants.html#None) = None, thls: [float](https://docs.python.org/3/library/functions.html#float) | TimeDependentScalar | [None](https://docs.python.org/3/library/constants.html#None) = None, z0mav: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, z0hav: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None, ps: [float](https://docs.python.org/3/library/functions.html#float) | TimeDependentScalar | [None](https://docs.python.org/3/library/constants.html#None) = None, albedoav: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None) = None)

Bases: [`SurfaceModule`](#modular_dales.Surface.surface.SurfaceModule)

Constant surface temperature boundary condition (isurf=2).

When added to a simulation, automatically enables isurf=2 in the namelist.
Expects surface config to contain ‘thls’ parameter.

* **Parameters:**
  * **sim** – Parent simulation instance
  * **isurf** – Surface scheme selector (2 for constant surface temperature)
  * **thls** – Surface liquid water potential temperature [K]
  * **z0mav** – Aerodynamic roughness length for momentum [m]
  * **z0hav** – Aerodynamic roughness length for heat and moisture [m]
  * **ps** – Surface pressure (Pa)
  * **albedoav** – Surface albedo (dimensionless)

## Namelist parameters

*Section* `NAMSURFACE`:
: - `albedoav` (field `albedoav`)
  - `isurf` (field `isurf`)
  - `ps` (field `ps`) (required)
  - `thls` (field `thls`) (required)
  - `z0hav` (field `z0hav`) (required)
  - `z0mav` (field `z0mav`) (required)

#### albedoav *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### check_settings()

Validate constant surface temperature settings.

#### do_config()

Configure surface settings and namelist.

#### isurf *: [int](https://docs.python.org/3/library/functions.html#int)* *= 2*

#### prepare_calculation()

Prepare data and setup for calculations.

#### prepare_calculations()

No additional preparation needed.

#### ps *: [float](https://docs.python.org/3/library/functions.html#float) | TimeDependentScalar | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### sim *: simulation_module | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### thls *: [float](https://docs.python.org/3/library/functions.html#float) | TimeDependentScalar | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### write_files()

No files to write.

#### z0hav *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

#### z0mav *: [float](https://docs.python.org/3/library/functions.html#float) | [None](https://docs.python.org/3/library/constants.html#None)* *= None*

### *class* modular_dales.Surface.surface.SurfaceModule(sim: [dales_simulation](modular_dales.dales_simulation.md#modular_dales.dales_simulation) = None)

Bases: `simulation_module`

Surface simulation module base class.
