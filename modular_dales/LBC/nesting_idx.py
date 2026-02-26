from dataclasses import dataclass
from typing import Optional

from modular_dales.Geometry import GridDales


@dataclass
class NestingIndices:
    ix_west: int
    ix_east: int
    iy_south: int
    iy_north: int
    supergrid_x0: float
    supergrid_y0: float
    subgrid_x0: float
    subgrid_y0: float
    supergrid: Optional["GridDales"] = None
