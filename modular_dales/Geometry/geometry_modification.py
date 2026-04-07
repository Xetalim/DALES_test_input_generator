import logging

import numpy as np

from modular_dales.Geometry import GridDales

logger = logging.getLogger(__name__)
logger.debug("Entered module: %s", __name__)


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

    def allGeometry(self):
        return np.ones_like(self.meshx, dtype=bool)

    def circleGeometry_realspace(self, x0, y0, size):
        return (self.meshx - x0) ** 2 + (self.meshy - y0) ** 2 <= size**2

    def rectangleGeometry_realspace(self, minx, maxx, miny, maxy):
        return (
            (self.meshx >= minx)
            & (self.meshx <= maxx)
            & (self.meshy >= miny)
            & (self.meshy <= maxy)
        )

    def rectangleGeometry_idxspace(self, minx, maxx, miny, maxy):
        return (
            (self.idxmesh >= minx)
            & (self.idxmesh <= maxx)
            & (self.idymesh >= miny)
            & (self.idymesh <= maxy)
        )

    def circleGeometry_idxspace(self, idx0, idy0, size):
        return (self.idxmesh - idx0) ** 2 + (self.idymesh - idy0) ** 2 <= size**2

    def parse_yaml_name(self, modification):
        """
        Parse YAML modification configuration and apply geometric transformation.

        Executes the specified geometry function based on the modification type,
        then applies the modification operation to the generated geometry.

        Args:
            modification (dict): Configuration dictionary containing:
                - "geometry" (str): Type of geometry to generate. Must be one of:
                  "circle_real", "all", "rectangle_real", "rectangle_idx", "circle_idx"
                - "params" (dict): Parameter dictionary passed to the geometry function

        Returns:
            None
        """
        param_dic = modification["params"]
        dic = {
            "circle_real": lambda: self.circleGeometry_idxspace(**param_dic),
            "all": lambda: self.allGeometry(),
            "rectangle_real": lambda: self.rectangleGeometry_realspace(**param_dic),
            "rectangle_idx": lambda: self.rectangleGeometry_idxspace(**param_dic),
            "circle_idx": lambda: self.circleGeometry_idxspace(**param_dic),
        }

        name = modification["geometry"]
        geometry_function = dic[name]
        geometry = geometry_function()

        self.do_modification(geometry, modification)

    def do_modification(self, geometry, modification):
        raise NotImplementedError(
            "Can't call do_modification of superclass, call the subclass instead"
        )
