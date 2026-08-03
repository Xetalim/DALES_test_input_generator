"""Module to enforce cross-section output for periodic DALES precursor runs."""

from dataclasses import dataclass, field
import logging
from typing import Optional

from modular_dales.modular.simulation_module import simulation_module
from modular_dales.MODULE_REGISTRY import register_module

logger = logging.getLogger(__name__)


@register_module
@dataclass
class PeriodicPrecursorCrossSections(simulation_module):
    """Force edge/middle cross-sections and a selected crossxy layer.

    This module configures ``namcrosssection`` so a periodic precursor run writes
    all cross-sections needed for turbulent-perturbation extraction:

    - ``crossyz`` at west edge, center, east edge
    - ``crossxz`` at south edge, center, north edge
    - ``crossxy`` at one configurable vertical level
    """

    sim: Optional["simulation_module"] = field(default=None, repr=False)
    top_layer_index: int = field(
        default=-1,
        repr=True,
        metadata={"serialize": True},
        init=True,
    )
    lcross: bool = field(
        default=True,
        metadata={"nml": "namcrosssection", "key": "lcross", "required": True},
        init=False,
    )

    def __post_init__(self):
        super().__init__(self.sim)
        self.module_name = "periodic_precursor_crosssections"

    def prepare_calculation(self):
        if self.grid is None:
            raise ValueError(
                "PeriodicPrecursorCrossSections requires an initialized simulation grid."
            )

        ix_west = 0
        ix_east = self.grid.itot - 1
        iy_south = 0
        iy_north = self.grid.jtot - 1
        ix_middle = self.grid.itot // 2
        iy_middle = self.grid.jtot // 2

        k_idx = int(self.top_layer_index)
        if k_idx < 0:
            k_idx += self.grid.kmax
        if k_idx < 0 or k_idx >= self.grid.kmax:
            raise IndexError(
                f"top_layer_index={self.top_layer_index} out of bounds for kmax={self.grid.kmax}."
            )

        existing_height = list(
            self.nml.get("namcrosssection", {}).get("crossheight", [])
        )
        existing_plane = list(self.nml.get("namcrosssection", {}).get("crossplane", []))
        existing_ortho = list(self.nml.get("namcrosssection", {}).get("crossortho", []))

        crossheight = list(dict.fromkeys(existing_height + [k_idx + 1]))
        crossplane = list(
            dict.fromkeys(existing_plane + [iy_south + 1, iy_middle + 1, iy_north + 1])
        )
        crossortho = list(
            dict.fromkeys(existing_ortho + [ix_west + 1, ix_middle + 1, ix_east + 1])
        )

        self.set_nml_section("namcrosssection", "crossheight", crossheight)
        self.set_nml_section("namcrosssection", "crossplane", crossplane)
        self.set_nml_section("namcrosssection", "crossortho", crossortho)

        logger.info(
            "Configured precursor cross-sections: crossortho=%s crossplane=%s crossheight=%s",
            crossortho,
            crossplane,
            crossheight,
        )

    def write_files(self):
        pass
