import logging
from dataclasses import dataclass, field
from typing import Any, Protocol, Union

import numpy as np

from modular_dales.Geometry import GridDales
from modular_dales.MODULE_REGISTRY import register_module, register_special_serializing

logger = logging.getLogger(__name__)
logger.debug("Entered module: %s", __name__)


@register_special_serializing
@register_module
@dataclass
class AllGeometry:
    """Select the full 2D horizontal domain."""

    def to_mask(self, modifier: "ModifierClass") -> np.ndarray:
        return np.ones_like(modifier.meshx, dtype=bool)


@register_special_serializing
@register_module
@dataclass
class CircleRealGeometry:
    """Circular mask in real-space coordinates (meters)."""

    x0: float
    y0: float
    size: float

    def to_mask(self, modifier: "ModifierClass") -> np.ndarray:
        return (modifier.meshx - self.x0) ** 2 + (
            modifier.meshy - self.y0
        ) ** 2 <= self.size**2


@register_special_serializing
@register_module
@dataclass
class RectangleRealGeometry:
    """Rectangular mask in real-space coordinates (meters)."""

    minx: float
    maxx: float
    miny: float
    maxy: float

    def to_mask(self, modifier: "ModifierClass") -> np.ndarray:
        return (
            (modifier.meshx >= self.minx)
            & (modifier.meshx <= self.maxx)
            & (modifier.meshy >= self.miny)
            & (modifier.meshy <= self.maxy)
        )


@register_special_serializing
@register_module
@dataclass
class RectangleIdxGeometry:
    """Rectangular mask in index-space coordinates."""

    minx: int
    maxx: int
    miny: int
    maxy: int

    def to_mask(self, modifier: "ModifierClass") -> np.ndarray:
        return (
            (modifier.idxmesh >= self.minx)
            & (modifier.idxmesh <= self.maxx)
            & (modifier.idymesh >= self.miny)
            & (modifier.idymesh <= self.maxy)
        )


@register_special_serializing
@register_module
@dataclass
class CircleIdxGeometry:
    """Circular mask in index-space coordinates."""

    idx0: int
    idy0: int
    size: float

    def to_mask(self, modifier: "ModifierClass") -> np.ndarray:
        return (modifier.idxmesh - self.idx0) ** 2 + (
            modifier.idymesh - self.idy0
        ) ** 2 <= self.size**2


@register_special_serializing
@register_module
@dataclass
class MaskGeometry:
    """Explicit boolean mask geometry.

    The provided mask must have shape ``(jtot, itot)`` after conversion to a
    numpy array and will be interpreted as a boolean selection.
    """

    mask: Any

    def to_mask(self, modifier: "ModifierClass") -> np.ndarray:
        mask_arr = np.asarray(self.mask)
        if mask_arr.shape != modifier.meshx.shape:
            raise ValueError(
                "MaskGeometry shape mismatch: "
                f"expected {modifier.meshx.shape}, got {mask_arr.shape}"
            )
        return mask_arr.astype(bool)


class GeometryLike(Protocol):
    """Protocol for typed geometry objects that produce a boolean mask."""

    def to_mask(self, modifier: "ModifierClass") -> np.ndarray: ...


GeometrySpec = Union[
    AllGeometry,
    CircleRealGeometry,
    RectangleRealGeometry,
    RectangleIdxGeometry,
    CircleIdxGeometry,
    MaskGeometry,
]


@dataclass
class GeometricModification:
    """Base class for modifications that target a typed horizontal geometry."""

    geometry: GeometrySpec = field(default_factory=AllGeometry)


class ModifierClass:
    """
    A class for generating and applying geometric modifications to DALES grid data.

    This class provides utilities to create various geometric shapes (circles, rectangles)
    in both real-space and index-space coordinates on a DALES grid, and to apply
    modifications based on these geometries.

    Attributes:
        x (ndarray): X-coordinates of the grid.
        y (ndarray): Y-coordinates of the grid.
        meshx (ndarray): 2D mesh grid of X-coordinates.
        meshy (ndarray): 2D mesh grid of Y-coordinates.
        idxmesh (ndarray): 2D mesh grid of X-indices.
        idymesh (ndarray): 2D mesh grid of Y-indices.

    Args:
        grid (GridDales): A GridDales object containing grid coordinate information.
    """

    def __init__(self, grid: "GridDales"):
        self.x = grid.xt
        self.y = grid.yt
        self.meshx, self.meshy = np.meshgrid(self.x, self.y)
        self.idxmesh, self.idymesh = np.meshgrid(
            np.arange(len(self.x)), np.arange(len(self.y))
        )
        self.grid = grid

    def apply_modification(self, modification: GeometricModification):
        """Apply a typed modification object to this modifier."""

        if not isinstance(modification, GeometricModification):
            raise TypeError(
                "Modification must inherit GeometricModification, got "
                f"{type(modification)}"
            )
        geometry_obj = modification.geometry
        if not hasattr(geometry_obj, "to_mask") or not callable(geometry_obj.to_mask):
            raise TypeError(
                "geometry must implement to_mask(modifier), got "
                f"{type(geometry_obj)}"
            )
        geometry = geometry_obj.to_mask(self)
        self.do_modification(geometry, modification)

    def do_modification(self, geometry, modification):
        raise NotImplementedError(
            "Can't call do_modification of superclass, call the subclass instead"
        )
