from helper_scripts.LBC.dales_nest.timestep0 import boundaries_timestep0
from helper_scripts.grids import GridDalesOpenBC, nesting_idx
from helper_scripts.LBC.dales_nest.get_timesteps0 import load_any_boundary_var

import xarray as xr


from pathlib import Path
import glob
import re
import numpy as np
import xarray as xr
from typing import Tuple, Dict, Generic
import logging
from helper_scripts.logging_wrapper import logwrap

logger = logging.getLogger(__name__)
logger.debug("Entered module: %s", __name__)


@logwrap
def get_all_dales_boundaries(input_json, grid: GridDalesOpenBC, indices: nesting_idx):
    """
    This function gets DALES boundaries from a host DALES simulation.

    :param input_json: Configuration dict for LBC
    :param grid: Grid object describing the output grid for open boundaries, which is different from the normal DALES grid.
    :type grid: GridDalesOpenBC
    :param indices: Indices describing where in the supergrid the edges of the current grid lie.
    :type indices: nesting_idx
    """
    ds_boundaries_timestep0 = boundaries_timestep0(input_json, grid, indices)
    # Get later time steps from corresponding coarse simulation output
    boundary_dict = {
        "west": (
            Path(input_json["outpath_coarse"]) / ".." / "crossyz.nc",
            {"xm": indices.ix_west + 1, "xt": indices.ix_west + 1},
        ),
        "east": (
            Path(input_json["outpath_coarse"]) / ".." / "crossyz.nc",
            {"xm": indices.ix_east + 1, "xt": indices.ix_east + 1},
        ),
        "south": (
            Path(input_json["outpath_coarse"]) / ".." / "crossxz.nc",
            {"ym": indices.iy_south + 1, "yt": indices.iy_south + 1},
        ),
        "north": (
            Path(input_json["outpath_coarse"]) / ".." / "crossxz.nc",
            {"ym": indices.iy_north + 1, "yt": indices.iy_north + 1},
        ),
        "top": (
            Path(input_json["outpath_coarse"]) / ".." / "crossxy.nc",
            {"zt": grid.kmax},
        ),
    }
    all_ls = []
    for boundary, (boundaryfile, sel_index) in boundary_dict.items():
        with xr.open_dataset(boundaryfile) as ds:
            for var in [
                "u",
                "v",
                "w",
                "thl",
                "qt",
                "e12",
                *input_json["tracernames"],
            ]:
                if var == "e12":
                    var_postfix = "0"
                else:
                    var_postfix = ""
                all_ls.append(
                    load_any_boundary_var(
                        ds.sel(sel_index),
                        var,
                        boundary=boundary,
                        grid=grid,
                        indices=indices,
                        var_postfix=var_postfix,
                    )
                )
    return xr.concat([ds_boundaries_timestep0, xr.merge(all_ls)], dim="time")
