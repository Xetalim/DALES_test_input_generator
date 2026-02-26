"""Immersed Boundary Method (IBM) module for building and geometry representation."""

from dataclasses import dataclass, field, asdict
import logging
from typing import Any, Dict, List, Optional
import netCDF4
import numpy as np

from modular_dales.Geometry import (
    ModifierClass,
    GridDales,
)
from modular_dales.modular import simulation_module
from modular_dales.MODULE_REGISTRY import register_module, register_singleton
from .download_AHN import get_process_ahn

logger = logging.getLogger(__name__)


@register_singleton
@register_module
@dataclass
class FromAHN:
    """AHN (Actueel Hoogtebestand Nederland) approach for IBM.
    Add this to IBMModule to generate IBM geometry from AHN data."""


@register_module
@dataclass
class IBMModification:
    """Single IBM modification. Set Geometry type and parameters to define the modification."""

    geometry: Optional[str] = None
    """Geometry type"""
    height: Optional[float] = None
    """Building height"""
    params: Dict[str, Any] = field(default_factory=dict)
    """Geometry-specific parameters"""


class IBMCreatorClass(ModifierClass):
    """
    A class for creating and managing IBM (Immersed Boundary Method) modifications to DALES grid geometry.

    This class handles the creation and output of boundary condition heights for immersed boundary
    method simulations. It maintains a 2D array of boundary condition heights that can be modified
    based on geometry specifications and exported to NetCDF format.

    Attributes:
        bc_height (np.ndarray): A 2D array storing boundary condition heights, initialized with zeros
                                matching the grid mesh dimensions.
    """

    def __init__(self, grid: GridDales):
        super().__init__(grid)
        self.bc_height = np.zeros_like(self.meshx)

    def do_modification(self, geometry, modification):
        self.bc_height[geometry] = modification["height"]

    def output_nc(self, filename):
        with netCDF4.Dataset(filename, "w") as nc:

            nc.createDimension("x", len(self.x))
            nc.createDimension("y", len(self.y))

            var_x = nc.createVariable("x", float, "x")
            var_y = nc.createVariable("y", float, "y")

            var_x[:] = self.x[:]
            var_y[:] = self.y[:]

            dims = ["y", "x"]
            bc_height_var = nc.createVariable("bc_height", float, dims)
            bc_height_var[:, :] = self.bc_height[:, :]


@register_module
@dataclass
class IBMModule(simulation_module):
    """IBM Module for DALES simulation.

    This module handles Immersed Boundary Method (IBM) configuration and setup for DALES
    (Dutch Atmospheric Large-Eddy Simulation) simulations. It manages IBM modifications,
    AHN (Actueel Hoogtebestand Nederland) terrain data integration, and generates IBM
    geometry files.

    Attributes:
        sim (Optional[simulation_module]): Reference to parent simulation module.
        ibm_modifications (List[IBMModification]): List of IBM geometry modifications to apply.
        ibm_generator (Any): IBM geometry creator instance, initialized during preparation.
        from_ahn (Optional[FromAHN]): AHN terrain data configuration object.
        ahn_zeroes_buffer (int): Buffer size for AHN zero elevation handling. Default: 5.
        subtract_ahn_mode (bool): Whether to subtract AHN elevation data. Default: True.
        apply_ibm (bool): Enable/disable IBM application. Default: True.
        ibas_prf (int): PRF advection scheme Default: 2.
        iadv_mom (int): Advection scheme for momentum. Default: 2.

    Methods:
        __post_init__: Initialize module and set module name.
        __add__: Add IBM modification or AHN configuration to module.
        __iadd__: In-place addition operator for fluent API.
        prepare_calculation: Set up IBM geometry, integrate AHN data if configured, and apply modifications.
        check_settings: --
        write_files: Export IBM geometry to NetCDF file.
    """

    sim: Optional["simulation_module"] = field(default=None, repr=False)
    ibm_modifications: List[IBMModification] = field(
        default_factory=list, init=True, repr=False, metadata={"serialize": True}
    )
    ibm_generator: Any = field(
        default=None, init=False, repr=False, metadata={"serialize": False}
    )
    from_ahn: Optional[FromAHN] = field(
        default=None,
        init=True,
        repr=False,
        metadata={"serialize": True},
    )
    ahn_zeroes_buffer: int = field(
        default=5,
        metadata={
            "serialize": True,
        },
        init=True,
        repr=False,
    )
    subtract_ahn_mode: bool = field(
        default=True,
        metadata={
            "serialize": True,
        },
        init=True,
        repr=False,
    )
    apply_ibm: bool = field(
        default=True,
        metadata={
            "nml": "NAMIBM",
            "key": "lapply_ibm",
            "serialize": True,
            "required": True,
        },
        init=True,
        repr=False,
    )
    ibas_prf: int = field(
        default=2,
        metadata={
            "nml": "DYNAMICS",
            "key": "IBAS_PRF",
            "serialize": True,
            "required": True,
        },
        init=True,
        repr=False,
    )
    iadv_mom: int = field(
        default=2,
        metadata={
            "nml": "DYNAMICS",
            "key": "IADV_MOM",
            "serialize": True,
            "required": True,
        },
        init=True,
        repr=False,
    )

    def __post_init__(self):
        super().__init__(self.sim)
        self.module_name = "IBMModule"

    def __add__(self, obj) -> "IBMModule":
        """Add IBM modification to module.

        Args:
            modification: IBMModification to add

        Returns:
            self for chaining
        """
        if isinstance(obj, FromAHN):
            self.from_ahn = obj
        elif isinstance(obj, IBMModification):
            self.ibm_modifications.append(obj)
        else:
            raise TypeError(
                f"Cannot add object of type {type(obj)} to IBMModule. Expected FromAHN or IBMModification."
            )
        return self

    def __iadd__(self, obj) -> "IBMModule":
        """In-place addition."""
        return self.__add__(obj)

    def prepare_calculation(self):
        """Set up IBM geometry and update namelist."""

        self.ibm_generator = IBMCreatorClass(self.grid)

        if self.from_ahn:
            ds = get_process_ahn(
                self.grid,
                zeroes_buffer=self.ahn_zeroes_buffer,
                subtract_ahn_mode=self.subtract_ahn_mode,
            )
            self.ibm_generator.bc_height[:, :] = ds[:, :]
        # Apply IBM modifications
        for modification in self.ibm_modifications:
            self.ibm_generator.parse_yaml_name(asdict(modification))

        logger.info("IBMModule: IBM configured and namelist updated")

    def check_settings(self):
        """Check IBM settings validity."""

        pass

    def write_files(self):
        """Write IBM files."""
        if self.ibm_generator is not None:
            self.ibm_generator.output_nc(
                self.output_path / "input" / f"ibm.inp_{self.exp_id:03d}.nc"
            )


@register_module
@dataclass
class IBMModifications:
    """Collection of IBM modifications."""

    modifications: List[IBMModification] = field(
        default_factory=list, metadata={"serialize": True}
    )

    def __add__(self, modification: IBMModification) -> "IBMModifications":
        """Add a modification."""
        self.modifications.append(modification)
        return self

    def __iadd__(self, modification: IBMModification) -> "IBMModifications":
        """In-place addition."""
        return self.__add__(modification)

    def apply_to_config(self, config: Dict[str, Any]) -> None:
        """Apply modifications to config dictionary."""
        if "ibm_modifications" not in config:
            config["ibm_modifications"] = []
        for mod in self.modifications:
            mod_dict = mod.to_dict()
            if mod_dict is not None:
                config["ibm_modifications"].append(mod_dict)
