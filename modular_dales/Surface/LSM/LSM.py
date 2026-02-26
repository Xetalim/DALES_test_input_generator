import logging
import pathlib
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Optional, Union

from modular_dales.modular.simulation_module import simulation_module
from modular_dales.MODULE_REGISTRY import register_module, register_singleton
from modular_dales.Surface.LSM.LSM_output_dales import LSM_output_dales, LsmModifier
from modular_dales.Surface.LSM.modular_temps_moisture import (
    SoilTemperatureMoistureFromHarmonie,
    UniformSkinTemperature,
    UniformSoilMoisture,
    UniformSoilTemperature,
    VaryingSkinTemperature,
    VaryingSoilMoisture,
    VaryingSoilTemperature,
)
from modular_dales.Surface.LSM import plot_lsm
from modular_dales.Surface.LSM.translation_tables import landuse_types
from modular_dales.Surface.surface import SurfaceModule

from .SLuRB.slurb import SLURBModule

logger = logging.getLogger(__name__)


@register_module
@dataclass
class LandUseModification:
    """Single land use modification."""

    geometry: str
    """Geometry type (all, circle_idx, rectangle_idx, etc.)"""
    type: str
    """Land use type (grs, urb, fbd, etc.)"""
    params: Dict[str, Any] = field(default_factory=dict)
    """Geometry-specific parameters"""


@register_module
@dataclass
class LandUseModifications:
    """Collection of land use modifications."""

    modifications: List[LandUseModification] = field(
        default_factory=list, metadata={"serialize": True}
    )

    def __add__(self, modification: LandUseModification) -> "LandUseModifications":
        """Add a modification."""
        self.modifications.append(modification)
        return self

    def __iadd__(self, modification: LandUseModification) -> "LandUseModifications":
        """In-place addition."""
        return self.__add__(modification)

    def apply_to_config(self, config: Dict[str, Any]) -> None:
        """Apply modifications to config dictionary."""
        if "land_use_modifications" not in config:
            config["land_use_modifications"] = []
        for mod in self.modifications:
            config["land_use_modifications"].append(mod.to_dict())


@register_singleton
@register_module
@dataclass
class FromLCZ:
    """LCZ (Local Climate Zone) classification approach for LSM."""


@register_module
@dataclass
class LSMModule(SurfaceModule):
    """Land Surface Model simulation module.

    When added to a simulation, automatically enables LSM by setting isurf=11
    in the namelist (unless overridden by a surface module like ConstantFluxesModule).

    Args:
        sim: Parent simulation instance
        isurf: Surface scheme selector (11 for LSM)
        ps: Surface pressure (Pa)
        z0mav: Momentum roughness length (m)
        z0hav: Heat roughness length (m)
        albedoav: Albedo (dimensionless)
        iinterp_t: Interpolation method for temperature (1-4)
        iinterp_theta: Interpolation method for potential temperature (1-4)
        dz_soil: List of soil layer thicknesses (m)
    """

    sim: Optional["simulation_module"] = field(default=None, repr=False)
    isurf: int = field(
        default=11, init=False, metadata={"nml": "NAMSURFACE", "key": "isurf"}
    )
    ps: Optional[float] = field(
        default=None,
        metadata={
            "nml": "NAMSURFACE",
            "key": "ps",
            "required": True,
            "serialize": True,
        },
    )
    z0mav: Optional[float] = field(
        default=None,
        metadata={
            "nml": "NAMSURFACE",
            "key": "z0mav",
            "required": True,
            "serialize": True,
        },
    )
    z0hav: Optional[float] = field(
        default=None,
        metadata={
            "nml": "NAMSURFACE",
            "key": "z0hav",
            "required": True,
            "serialize": True,
        },
    )
    albedoav: Optional[float] = field(
        default=None,
        metadata={"nml": "NAMSURFACE", "key": "albedoav", "serialize": True},
    )
    iinterp_t: Optional[int] = field(
        default=None,
        metadata={
            "nml": "NAMLSM",
            "key": "iinterp_t",
            "required": True,
            "serialize": True,
        },
    )
    iinterp_theta: Optional[int] = field(
        default=None,
        metadata={
            "nml": "NAMLSM",
            "key": "iinterp_theta",
            "required": True,
            "serialize": True,
        },
    )
    kmax_soil: Optional[int] = field(
        default=4,
        metadata={
            "nml": "DOMAIN",
            "key": "kmax_soil",
            "required": True,
            "serialize": True,
        },
    )
    dz_soil: Optional[List[float]] = field(
        default=None,
        metadata={
            "nml": "NAMLSM",
            "key": "dz_soil",
            "required": True,
            "serialize": True,
        },
    )
    lheterogeneous: bool = field(
        default=True,
        metadata={"nml": "NAMLSM", "key": "lheterogeneous", "serialize": False},
    )
    nlu: int = field(
        default=0,
        metadata={"nml": "NAMLSM", "key": "nlu", "serialize": False},
        init=False,
        repr=False,
    )
    land_use_modifications: LandUseModifications = field(
        default_factory=LandUseModifications,
        init=True,
        repr=True,
        metadata={"serialize": True},
    )
    from_lcz: Optional[FromLCZ] = field(
        default=None,
        init=True,
        repr=False,
        metadata={"serialize": True},
    )
    skin_temperature: Optional[
        Union[UniformSkinTemperature, VaryingSkinTemperature]
    ] = field(
        default=None,
        init=True,
        repr=True,
        metadata={"serialize": True},
    )
    soil_temperature: Optional[
        Union[
            UniformSoilTemperature,
            VaryingSoilTemperature,
            SoilTemperatureMoistureFromHarmonie,
        ]
    ] = field(
        default=None,
        init=True,
        repr=True,
        metadata={"serialize": True},
    )
    soil_moisture: Optional[
        Union[
            UniformSoilMoisture,
            VaryingSoilMoisture,
            SoilTemperatureMoistureFromHarmonie,
        ]
    ] = field(
        default=None,
        init=True,
        repr=True,
        metadata={"serialize": True},
    )
    lsm_writer: Optional[LSM_output_dales] = field(
        default=None, init=False, repr=False, metadata={"serialize": False}
    )
    slurb_module: Optional[SLURBModule] = field(
        default=None, init=False, repr=False, metadata={"serialize": False}
    )

    def __post_init__(self):
        super().__init__(self.sim)
        self.module_name = "LSMModule"

    def do_config(self):
        """Configure namelist and surface configuration for LSM."""

        # Validate iinterp settings
        for param_name in ["iinterp_t", "iinterp_theta"]:
            value = getattr(self, param_name)
            if value is None:
                raise ValueError(f"LSMModule requires '{param_name}' parameter")
            if not isinstance(value, int):
                raise ValueError(f"{param_name} must be an integer")
            if not 1 <= value <= 4:
                raise ValueError(
                    f"{param_name} must be an integer between 1 and 4 (1=arithmetic mean, 2=geometric mean, 3=harmonic mean, 4=max)"
                )

        logger.info("LSMModule: NAMSURFACE/NAMLSM configured from dataclass fields")

    def __add__(self, obj) -> "LSMModule":
        """Add configurations to LSM module.

        Args:
            obj: LandUseModification, LandUseModifications or FromLCZ

        Returns:
            self for chaining
        """
        if isinstance(obj, LandUseModification):
            self.land_use_modifications += obj
        elif isinstance(obj, LandUseModifications):
            for mod in obj.modifications:
                self.land_use_modifications += mod
        elif isinstance(obj, FromLCZ):
            self.from_lcz = obj
        elif isinstance(obj, (UniformSkinTemperature, VaryingSkinTemperature)):
            self.skin_temperature = obj
        elif isinstance(obj, (UniformSoilTemperature, VaryingSoilTemperature)):
            self.soil_temperature = obj
        elif isinstance(
            obj,
            (
                UniformSoilMoisture,
                VaryingSoilMoisture,
            ),
        ):
            self.soil_moisture = obj
        elif isinstance(obj, SoilTemperatureMoistureFromHarmonie):
            self.soil_moisture = obj
            self.soil_temperature = obj
        else:
            raise TypeError(
                f"Expected LandUseModification/FromLCZ/SkinTemperatures/SoilTemperatures/SoilMoistures, got {type(obj)}"
            )
        return self

    def __iadd__(self, modification) -> "LSMModule":
        """In-place addition of modification."""
        return self.__add__(modification)

    def prepare_calculation(self):
        return self.prepare_calculations()

    def prepare_calculations(self):
        """Set up LSM with land use types and soil configuration."""

        self.exists_soil_temp_moisture_skin_temp()

        if self.module_exists(SLURBModule):
            self.slurb_module = self.retrieve_module(SLURBModule)

        # Determine land use types
        if self.from_lcz is None:
            lu_types = landuse_types.lu_types_depac.copy()
            if self.slurb_module:
                lu_types["slb"] = landuse_types.slb
        else:
            lu_types = landuse_types.lu_types_ifs

        # Initialize LSM writer
        self.lsm_writer = LSM_output_dales(
            self.grid,
            lu_types,
            soil_levels=self.kmax_soil,
        )

        # Fill standard geometry
        self.lsm_writer.standard_fill_geometry_modification()

        # Configure based on LCZ or standard approach
        if self.from_lcz is None:
            # Apply land use modifications
            modifier = LsmModifier(self.lsm_writer, self.grid)
            for modification in self.land_use_modifications.modifications:
                modifier.parse_yaml_name(asdict(modification))
            self.lsm_writer = modifier.lsm_input
            self.lsm_writer.init_lutypes_ifs()
            self.lsm_writer.recalculate_remaining_cover()
        else:
            self.lsm_writer.from_lcz()

            # Apply SLURB parameters if SLURB is enabled and module is present
            if self.slurb_module is not None:
                self.slurb_module.init_generator()  # Ensure SLURB generator is initialized
                self.lsm_writer.apply_slurb_parameters_lcz(
                    self.slurb_module.slb_generator
                )
        # Apply soil temperature profiles, skin temperature, soil moisture profiles
        self.apply_soil_temp_moisture_skin_temp()
        # Finalize
        self.lsm_writer.trim_landuse()

        self.nlu = self.lsm_writer.nlu + 1

    def exists_soil_temp_moisture_skin_temp(self):
        # check if soil temp exists
        if self.soil_temperature is None:
            raise ValueError(
                "Soil temperature profiles must be specified for LSMModule. "
                "Add a SoilTemperatures object to LSMModule: lsm += SoilTemperatures()"
            )
        # check if soil moisture exists
        if self.soil_moisture is None:
            raise ValueError(
                "Soil moisture profiles must be specified for LSMModule. "
                "Add a SoilMoistures object to LSMModule: lsm += SoilMoistures()"
            )
        # check if skin temp exists
        if self.skin_temperature is None:
            raise ValueError(
                "Skin temperature profiles must be specified for LSMModule. "
                "Add a SkinTemperatures object to LSMModule: lsm += SkinTemperatures()"
            )

    def apply_soil_temp_moisture_skin_temp(self):
        if self.lsm_writer is None:
            raise ValueError(
                "Can't apply soil temperature/moisture/skin temp without initializing LSM writer first"
            )
        if isinstance(self.soil_moisture, (UniformSoilMoisture)):
            self.lsm_writer.set_uniform_soil_moisture(self.soil_moisture.soil_moisture)
        elif isinstance(self.soil_moisture, (VaryingSoilMoisture)):
            self.lsm_writer.set_soil_moisture_array(self.soil_moisture.soil_moisture)
        elif isinstance(self.soil_moisture, (SoilTemperatureMoistureFromHarmonie)):
            self.lsm_writer.set_soil_moisture_array(
                self.soil_moisture.get_soil_moisture_array(self.grid, self.dz_soil)
            )
        else:
            raise ValueError(
                "Invalid soil moisture configuration. Must be UniformSoilMoisture, VaryingSoilMoisture, or SoilTemperatureMoistureFromHarmonie."
            )

        if isinstance(
            self.soil_temperature,
            (UniformSoilTemperature),
        ):
            self.lsm_writer.set_uniform_soil_temperature(
                self.soil_temperature.soil_temperature
            )
        elif isinstance(
            self.soil_temperature,
            (VaryingSoilTemperature),
        ):
            self.lsm_writer.set_soil_temperature_array(
                self.soil_temperature.soil_temperature
            )
        elif isinstance(self.soil_temperature, (SoilTemperatureMoistureFromHarmonie)):
            self.lsm_writer.set_soil_temperature_array(
                self.soil_temperature.get_soil_temperature_array(
                    self.grid, self.dz_soil
                )
            )
        else:
            raise ValueError(
                "Invalid soil temperature configuration. Must be UniformSoilTemperature, VaryingSoilTemperature, or SoilTemperatureMoistureFromHarmonie."
            )

        if isinstance(self.skin_temperature, (UniformSkinTemperature)):
            self.lsm_writer.set_skin_temperature(
                self.skin_temperature.skin_temperature, lu_type="all"
            )
        elif isinstance(self.skin_temperature, (VaryingSkinTemperature)):
            self.lsm_writer.set_skin_temperature_array(
                self.skin_temperature.skin_temperature, lu_type="all"
            )
        else:
            raise ValueError(
                "Invalid skin temperature configuration. Must be UniformSkinTemperature or VaryingSkinTemperature."
            )

    def check_settings(self):
        """Check LSM settings validity."""
        self.sim.required_files["van_genuchten_parameters.nc"] = (
            pathlib.Path(self.sim.machine_conf["case_conf"]["SOURCE_PATH"])
            / "data"
            / "van_genuchten_parameters.nc"
        )

    def write_files(self):
        """Write LSM input files and generate plots."""

        self.lsm_writer.save_netcdf(self.output_path / "input", self.exp_id)
        plot_lsm.plot_lsm_cover(
            pathlib.Path(self.output_path) / "input" / f"lsm.inp_{self.exp_id:03d}.nc",
            pathlib.Path(self.output_path) / "profiles",
        )
