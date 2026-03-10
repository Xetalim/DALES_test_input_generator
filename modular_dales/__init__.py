"""Public API for the modular_dales package.

This module exposes the main high-level classes and helpers that are
intended for users of the library, while still keeping the internal
package structure available for advanced use.
"""

from .modular.time_dependent_scalars import TimeDependentScalar
from .MODULE_REGISTRY import (
    MODULE_REGISTRY,
    SINGLETON_REGISTRY,
    register_module,
    register_singleton,
)

from .modular.dales_simulation import dales_simulation
from .modular.simulation_module import simulation_module

from .Geometry.GridDales import GridDales, GridDalesOpenBC

from .logging_wrapper import logwrap, setup_logging

from .modular.time_dependent import (
    TimedependentModule,
)
from .vars import VariableDefinition, get_all_vars
from .Atmosphere import AtmosphereModule, AtmosphericProfile, InterpolatedProfile

from .LBC import (
    do_openboundary,
    Nest_in_Dales,
    Nest_in_AtmosphereProfiles,
    NestingTopology,
)

from .Configuration.defaultnamelist import DefaultNamelistModule
from .Configuration.run_and_time import RunModule, TimeModule
from .Configuration.output_modules import EasyOutputModule

from .Geometry.geometry_modification import ModifierClass

from .Surface.surface import (
    SurfaceModule,
    ConstantFluxesModule,
    ConstantFluxesWithShearModule,
    ConstantSurfaceTemperatureModule,
)
from .Surface.LSM.LSM import (
    LSMModule,
    LandUseModification,
    LandUseModifications,
    FromLCZ,
)
from .Surface.LSM.modular_temps_moisture import (
    UniformSkinTemperature,
    UniformSoilTemperature,
    UniformSoilMoisture,
    VaryingSkinTemperature,
    VaryingSoilTemperature,
    VaryingSoilMoisture,
)
from .Surface.LSM.SLuRB.slurb import (
    SLURBModule,
    SLURBModification,
    SLURBModifications,
)

from .Emission.emission import (
    EmissionModule,
    EmissionTracer,
    EmissionPointSource,
)

from .Radiation.radiation import RadiationModule

from .IBM.IBM import (
    IBMModule,
    IBMModification,
    IBMModifications,
    FromAHN,
)

__all__ = [
    # Core simulation framework
    "dales_simulation",
    "simulation_module",
    "TimeDependentScalar",
    "TimedependentModule",
    # Configuration modules
    "DefaultNamelistModule",
    "RunModule",
    "TimeModule",
    "EasyOutputModule",
    # Geometry
    "GridDales",
    "GridDalesOpenBC",
    "ModifierClass",
    # AtmosphereModule, AtmosphericProfile, InterpolatedProfile
    "AtmosphereModule",
    "AtmosphericProfile",
    "InterpolatedProfile",
    # Surface and LSM
    "SurfaceModule",
    "ConstantFluxesModule",
    "ConstantFluxesWithShearModule",
    "ConstantSurfaceTemperatureModule",
    "LSMModule",
    "LandUseModification",
    "LandUseModifications",
    "FromLCZ",
    "UniformSkinTemperature",
    "UniformSoilTemperature",
    "UniformSoilMoisture",
    "VaryingSkinTemperature",
    "VaryingSoilTemperature",
    "VaryingSoilMoisture",
    "SLURBModule",
    "SLURBModification",
    "SLURBModifications",
    # Lateral boundary conditions
    "do_openboundary",
    "Nest_in_Dales",
    "Nest_in_AtmosphereProfiles",
    "NestingTopology",
    # Emissions
    "EmissionModule",
    "EmissionTracer",
    "EmissionPointSource",
    # Radiation
    "RadiationModule",
    # IBM
    "IBMModule",
    "IBMModification",
    "IBMModifications",
    "FromAHN",
    # LBC / nesting
    "NestingTopology",
    "do_openboundary",
    "Nest_in_Dales",
    # variables
    "VariableDefinition",
    "get_all_vars",
    # Registry and logging helpers
    "MODULE_REGISTRY",
    "SINGLETON_REGISTRY",
    "register_module",
    "register_singleton",
    "logwrap",
    "setup_logging",
]
