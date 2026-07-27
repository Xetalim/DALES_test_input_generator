"""Additional non-output settings modules for DALES namelist configuration."""

from dataclasses import dataclass, field
from typing import Optional
import numpy as np

from modular_dales.MODULE_REGISTRY import register_module
from modular_dales.modular.dales_simulation import dales_simulation
from modular_dales.modular.simulation_module import simulation_module


@register_module
@dataclass
class SprayingModule(simulation_module):
    """Settings module for sea-spray and water-spray source terms.

    Supports both index-based and coordinate-based source placement.
    If x_spray/y_spray/z_spray are all provided, coordinates take precedence and
    are mapped to nearest grid centers. If no coordinates are provided, the
    configured i_glob_spray/j_glob_spray/k_glob_spray indices are used.
    """

    sim: Optional["dales_simulation"] = field(default=None, repr=False)
    lwater_spraying: bool = field(
        default=False,
        metadata={
            "nml": "namspraying",
            "key": "lwater_spraying",
            "doc": "Enable prescribed water spraying source.",
        },
    )
    lsalt_spraying: bool = field(
        default=False,
        metadata={
            "nml": "namspraying",
            "key": "lsalt_spraying",
            "doc": "Enable prescribed salt spraying source.",
        },
    )
    i_glob_spray: int = field(
        default=1,
        metadata={
            "nml": "namspraying",
            "key": "i_glob_spray",
            "doc": "Global i-index location for spray source.",
        },
    )
    j_glob_spray: int = field(
        default=1,
        metadata={
            "nml": "namspraying",
            "key": "j_glob_spray",
            "doc": "Global j-index location for spray source.",
        },
    )
    k_glob_spray: int = field(
        default=1,
        metadata={
            "nml": "namspraying",
            "key": "k_glob_spray",
            "doc": "Vertical level index for spray source.",
        },
    )
    x_spray: Optional[float] = field(
        default=None,
        metadata={
            "doc": "Physical x-coordinate [m] of spray source; converted to i_glob_spray using GridModule xt centers."
        },
    )
    y_spray: Optional[float] = field(
        default=None,
        metadata={
            "doc": "Physical y-coordinate [m] of spray source; converted to j_glob_spray using GridModule yt centers."
        },
    )
    z_spray: Optional[float] = field(
        default=None,
        metadata={
            "doc": "Physical z-coordinate [m] of spray source; converted to k_glob_spray using GridModule zt centers."
        },
    )
    water_spray_rate: float = field(
        default=0.0,
        metadata={
            "nml": "namspraying",
            "key": "water_spray_rate",
            "doc": "Injected water spray mass rate.",
        },
    )
    salt_spray_rate: float = field(
        default=0.0,
        metadata={
            "nml": "namspraying",
            "key": "salt_spray_rate",
            "doc": "Injected salt spray mass rate.",
        },
    )
    salinity: float = field(
        default=0.0,
        metadata={
            "nml": "namspraying",
            "key": "salinity",
            "doc": "Salt mass fraction used for water spray.",
        },
    )
    tracer: int = field(
        default=1,
        metadata={
            "nml": "namspraying",
            "key": "tracer",
            "doc": "Tracer index that receives spray source.",
        },
    )
    lsalt_sponge: bool = field(
        default=False,
        metadata={
            "nml": "namspraying",
            "key": "lsalt_sponge",
            "doc": "Enable salt sponge damping around spray source.",
        },
    )
    lcoupled: bool = field(
        default=False,
        metadata={
            "nml": "namspraying",
            "key": "lcoupled",
            "doc": "Enable coupling with external spray forcing.",
        },
    )
    target_mode: int = field(
        default=0,
        metadata={
            "nml": "namspraying",
            "key": "target_mode",
            "doc": "Spray forcing mode selector.",
        },
    )

    def __post_init__(self):
        super().__init__(self.sim)
        self.module_name = "SprayingModule"

    def do_config(self):
        return None

    def prepare_calculation(self):
        provided_coord_count = sum(
            v is not None for v in (self.x_spray, self.y_spray, self.z_spray)
        )
        if provided_coord_count == 0:
            return None

        if provided_coord_count != 3:
            raise ValueError(
                "SprayingModule: provide either all of x_spray/y_spray/z_spray or none."
            )

        if self.grid is None:
            raise ValueError(
                "SprayingModule: x_spray/y_spray/z_spray require GridModule to be configured first."
            )

        if self.x_spray is not None:
            x_centers = np.asarray(self.grid.xt)
            if self.x_spray < float(x_centers[0]) or self.x_spray > float(
                x_centers[-1]
            ):
                raise ValueError(
                    f"SprayingModule: x_spray={self.x_spray} m is outside domain center range [{float(x_centers[0])}, {float(x_centers[-1])}] m."
                )
            self.i_glob_spray = int(np.argmin(np.abs(x_centers - self.x_spray))) + 1

        if self.y_spray is not None:
            y_centers = np.asarray(self.grid.yt)
            if self.y_spray < float(y_centers[0]) or self.y_spray > float(
                y_centers[-1]
            ):
                raise ValueError(
                    f"SprayingModule: y_spray={self.y_spray} m is outside domain center range [{float(y_centers[0])}, {float(y_centers[-1])}] m."
                )
            self.j_glob_spray = int(np.argmin(np.abs(y_centers - self.y_spray))) + 1

        if self.z_spray is not None:
            z_centers = np.asarray(self.grid.zt)
            if self.z_spray < float(z_centers[0]) or self.z_spray > float(
                z_centers[-1]
            ):
                raise ValueError(
                    f"SprayingModule: z_spray={self.z_spray} m is outside domain center range [{float(z_centers[0])}, {float(z_centers[-1])}] m."
                )
            self.k_glob_spray = int(np.argmin(np.abs(z_centers - self.z_spray))) + 1

        return None

    def check_settings(self):
        return None

    def write_files(self):
        return None


@register_module
@dataclass
class LateralSpongeModule(simulation_module):
    """Settings module for lateral sponge damping configuration."""

    sim: Optional["dales_simulation"] = field(default=None, repr=False)
    llateral_sponge: bool = field(
        default=False,
        metadata={
            "nml": "lateral_sponge",
            "key": "llateral_sponge",
            "doc": "Enable lateral sponge damping near side boundaries.",
        },
    )
    nudgedepth: float = field(
        default=0.0,
        metadata={
            "nml": "lateral_sponge",
            "key": "nudgedepth",
            "doc": "Depth/width scale of lateral sponge region in m.",
        },
    )

    def __post_init__(self):
        super().__init__(self.sim)
        self.module_name = "LateralSpongeModule"

    def do_config(self):
        return None

    def prepare_calculation(self):
        return None

    def check_settings(self):
        return None

    def write_files(self):
        return None


@register_module
@dataclass
class GeneralPhysicsModule(simulation_module):
    """Settings module combining thermodynamics and subgrid controls."""

    sim: Optional["dales_simulation"] = field(default=None, repr=False)

    lmoist: Optional[bool] = field(
        default=None,
        metadata={
            "nml": "thermodynamics",
            "key": "lmoist",
            "doc": "Enable moist thermodynamics.",
        },
    )
    chi_half: Optional[bool] = field(
        default=None,
        metadata={
            "nml": "thermodynamics",
            "key": "chi_half",
            "doc": "Use half-level chi/exner treatment where applicable.",
        },
    )
    lconstexner: Optional[bool] = field(
        default=None,
        metadata={
            "nml": "thermodynamics",
            "key": "lconstexner",
            "doc": "Use initial pressure profile in exner function.",
        },
    )
    lbaseexner: Optional[bool] = field(
        default=None,
        metadata={
            "nml": "thermodynamics",
            "key": "lbaseexner",
            "doc": "Use base-state pressure profile in exner function.",
        },
    )
    lnoclouds: Optional[bool] = field(
        default=None,
        metadata={
            "nml": "thermodynamics",
            "key": "lnoclouds",
            "doc": "Disable cloud liquid calculations in thermodynamics.",
        },
    )
    lqlnr: Optional[bool] = field(
        default=None,
        metadata={
            "nml": "thermodynamics",
            "key": "lqlnr",
            "doc": "Toggle non-rain liquid handling option.",
        },
    )

    ldelta: Optional[bool] = field(
        default=None,
        metadata={
            "nml": "NAMSUBGRID",
            "key": "ldelta",
            "doc": "Use delta-based subgrid length scale.",
        },
    )
    lmason: Optional[bool] = field(
        default=None,
        metadata={
            "nml": "NAMSUBGRID",
            "key": "lmason",
            "doc": "Enable Mason wall-damping formulation.",
        },
    )
    cf: Optional[float] = field(
        default=None,
        metadata={
            "nml": "NAMSUBGRID",
            "key": "cf",
            "doc": "Subgrid closure constant cf.",
        },
    )
    cn: Optional[float] = field(
        default=None,
        metadata={
            "nml": "NAMSUBGRID",
            "key": "cn",
            "doc": "Subgrid closure constant cn.",
        },
    )
    Rigc: Optional[float] = field(
        default=None,
        metadata={
            "nml": "NAMSUBGRID",
            "key": "Rigc",
            "doc": "Critical Richardson number for subgrid closure.",
        },
    )
    Prandtl: Optional[float] = field(
        default=None,
        metadata={
            "nml": "NAMSUBGRID",
            "key": "Prandtl",
            "doc": "Subgrid turbulent Prandtl number.",
        },
    )
    lsmagorinsky: Optional[bool] = field(
        default=None,
        metadata={
            "nml": "NAMSUBGRID",
            "key": "lsmagorinsky",
            "doc": "Enable Smagorinsky-type closure.",
        },
    )
    cs: Optional[float] = field(
        default=None,
        metadata={"nml": "NAMSUBGRID", "key": "cs", "doc": "Smagorinsky coefficient."},
    )
    nmason: Optional[int] = field(
        default=None,
        metadata={
            "nml": "NAMSUBGRID",
            "key": "nmason",
            "doc": "Selector for Mason damping variant.",
        },
    )
    sgs_surface_fix: Optional[bool] = field(
        default=None,
        metadata={
            "nml": "NAMSUBGRID",
            "key": "sgs_surface_fix",
            "doc": "Enable near-surface SGS correction.",
        },
    )
    ch1: Optional[float] = field(
        default=None,
        metadata={
            "nml": "NAMSUBGRID",
            "key": "ch1",
            "doc": "Empirical closure coefficient ch1.",
        },
    )
    lanisotrop: Optional[bool] = field(
        default=None,
        metadata={
            "nml": "NAMSUBGRID",
            "key": "lanisotrop",
            "doc": "Enable anisotropic diffusion treatment.",
        },
    )
    lD80R: Optional[bool] = field(
        default=None,
        metadata={
            "nml": "NAMSUBGRID",
            "key": "lD80R",
            "doc": "Enable D80R subgrid option.",
        },
    )

    def __post_init__(self):
        super().__init__(self.sim)
        self.module_name = "GeneralPhysicsModule"

    def do_config(self):
        return None

    def prepare_calculation(self):
        return None

    def check_settings(self):
        return None

    def write_files(self):
        return None
