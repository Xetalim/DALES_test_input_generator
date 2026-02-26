"""Radiation module for solar and thermal radiation handling."""

from dataclasses import dataclass, field
import pathlib
import logging
from typing import Optional

from modular_dales.MODULE_REGISTRY import register_module
from modular_dales.modular.simulation_module import simulation_module
from modular_dales.Surface.surface import SurfaceModule

logger = logging.getLogger(__name__)


@register_module
@dataclass
class RadiationModule(simulation_module):
    """Radiation simulation module.

    Possible values for `iradiation`:
        0: no radiation
        1: full radiation
        2: parameterized radiation
        3: simple surface radiation for land surface model
        4: RRTMG radiation
        5: RTE-RRTMGP radiation
        10: user specified radiation

    Args:
        sim: Parent dales_simulation instance
        iradiation: Radiation type
        ssa: Representative single scattering albedo (0 <= x <= 1)
        ide: Scalar field used as aerosols if laero set to .true.
        laero: .true. for aerosols, .false. for clouds
        lCnstZenith: Switch to apply a fixed solar zenith angle
        ioverlap: Flag for cloud overlap method
        inflglw: Flag for RRTMG longwave input
        iceflglw: Flag for ice particle specification in longwave
        liqflglw: Flag for effect of liquid water in longwave
        inflgsw: Flag for RRTMG shortwave input
        iceflgsw: Flag for ice particle specification in shortwave
        liqflgsw: Flag for effect of liquid water in shortwave
        iyear: Year of the simulation
        ocean: Switch to calculate radiation over ocean
        nbatch: Number of batch of vertical columns sent to RTE-RRTMGP routines
        usepade: Use Pade coefficients for cloud optical properties instead of lookup tables
        doclearsky: Use clear sky radiation in the calculation
    """

    sim: Optional["simulation_module"] = field(default=None, repr=False)
    iradiation: Optional[int] = field(
        default=None, metadata={"nml": "PHYSICS", "key": "IRADIATION"}
    )
    # NAMDE parameters
    ssa: Optional[float] = field(default=None, metadata={"nml": "NAMDE", "key": "ssa"})
    laero: Optional[float] = field(
        default=None, metadata={"nml": "NAMDE", "key": "laero"}
    )
    ide: Optional[int] = field(default=None, metadata={"nml": "NAMDE", "key": "ide"})

    # NAMRADIATION parameters
    lCnstZenith: Optional[bool] = field(
        default=None, metadata={"nml": "NAMRADIATION", "key": "lCnstZenith"}
    )
    ioverlap: Optional[int] = field(
        default=None, metadata={"nml": "NAMRADIATION", "key": "ioverlap"}
    )
    inflglw: Optional[int] = field(
        default=None, metadata={"nml": "NAMRADIATION", "key": "inflglw"}
    )
    iceflglw: Optional[int] = field(
        default=None, metadata={"nml": "NAMRADIATION", "key": "iceflglw"}
    )
    liqflglw: Optional[int] = field(
        default=None, metadata={"nml": "NAMRADIATION", "key": "liqflglw"}
    )
    inflgsw: Optional[int] = field(
        default=None, metadata={"nml": "NAMRADIATION", "key": "inflgsw"}
    )
    iceflgsw: Optional[int] = field(
        default=None, metadata={"nml": "NAMRADIATION", "key": "iceflgsw"}
    )
    liqflgsw: Optional[int] = field(
        default=None, metadata={"nml": "NAMRADIATION", "key": "liqflgsw"}
    )
    iyear: Optional[int] = field(
        default=None, metadata={"nml": "NAMRADIATION", "key": "iyear"}
    )
    ocean: Optional[bool] = field(
        default=None, metadata={"nml": "NAMRADIATION", "key": "ocean"}
    )

    # NAMRTERRMTGP parameters
    nbatch: Optional[int] = field(
        default=None, metadata={"nml": "NAMRTERRTMGP", "key": "nbatch"}
    )
    usepade: Optional[bool] = field(
        default=None, metadata={"nml": "NAMRTERRTMGP", "key": "usepade"}
    )
    doclearsky: Optional[bool] = field(
        default=None, metadata={"nml": "NAMRTERRTMGP", "key": "doclearsky"}
    )

    timerad: int = field(
        default=60,
        metadata={
            "nml": "PHYSICS",
            "key": "timerad",
            "required": True,
            "serialize": True,
        },
        init=True,
    )

    surface_module: Optional["SurfaceModule"] = field(
        default=None, init=False, repr=False, metadata={"serialize": False}
    )

    def __post_init__(self):
        super().__init__(self.sim)
        self.module_name = "RadiationModule"

    def do_config(self):
        """Ensure radiation configuration is set."""
        return None

    def prepare_calculation(self):
        """No additional preparation needed."""
        if self.module_exists(SurfaceModule):
            self.surface_module = self.retrieve_module(SurfaceModule)
        else:
            raise ValueError("RadiationModule requires a SurfaceModule.")

        return None

    def check_settings(self):
        """Validate constant fluxes settings."""
        pass

    def write_files(self):
        iradiation = self.iradiation or 0
        if iradiation != 0 and self.surface_module is not None:
            if getattr(self.surface_module, "albedoav", None) is None:
                raise ValueError(
                    "RadiationModule: albedoav must be set in surface config for radiation to work properly"
                )
        if iradiation != 0 and self.surface_module is None:
            raise ValueError(
                "RadiationModule: albedoav must be set in surface config for radiation to work properly"
            )
        exp_id = self.exp_id
        machine_conf = self.sim.machine_conf
        if iradiation == 1:
            self.sim.required_files[f"backrad.inp.{exp_id:03d}"] = (
                pathlib.Path(machine_conf["case_conf"]["SOURCE_PATH"])
                / "cases"
                / "example"
                / "backrad.inp.001"
            ).as_posix()
        if iradiation == 4:
            # this is an RRTMG case, we need RRTMG_LW and RRTMG_SW and backrad.inp.001.nc
            self.sim.required_files[f"backrad.inp.{exp_id:03d}.nc"] = (
                pathlib.Path.cwd() / "extra_data/backrad.inp.001.nc"
            ).as_posix()
            self.sim.required_files["rrtmg_lw.nc"] = (
                pathlib.Path(machine_conf["case_conf"]["SOURCE_PATH"])
                / "external"
                / "RRTMG"
                / "RRTMG_LW"
                / "data"
                / "rrtmg_lw.nc"
            ).as_posix()
            self.sim.required_files["rrtmg_sw.nc"] = (
                pathlib.Path(machine_conf["case_conf"]["SOURCE_PATH"])
                / "external"
                / "RRTMG"
                / "RRTMG_SW"
                / "data"
                / "rrtmg_sw.nc"
            ).as_posix()
        elif iradiation == 5:
            # this is RTE_RRTMG, we need all data from RTE_RRTMG
            self.sim.required_files[f"backrad.inp.{exp_id:03d}.nc"] = (
                pathlib.Path.cwd() / "extra_data/backrad.inp.001.nc"
            ).as_posix()
            self.sim.required_files["rrtmg_lw.nc"] = (
                pathlib.Path(machine_conf["case_conf"]["SOURCE_PATH"])
                / "external"
                / "RRTMG"
                / "RRTMG_LW"
                / "data"
                / "rrtmg_lw.nc"
            ).as_posix()
            self.sim.required_files["rrtmg_sw.nc"] = (
                pathlib.Path(machine_conf["case_conf"]["SOURCE_PATH"])
                / "external"
                / "RRTMG"
                / "RRTMG_SW"
                / "data"
                / "rrtmg_sw.nc"
            ).as_posix()
            for file in (
                pathlib.Path(machine_conf["case_conf"]["SOURCE_PATH"])
                / "external"
                / "rrtmgp-data"
            ).glob("*.nc"):
                self.sim.required_files[file.name] = file.as_posix()
        return None
