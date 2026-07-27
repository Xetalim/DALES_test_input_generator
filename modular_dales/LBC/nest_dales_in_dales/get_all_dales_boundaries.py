import glob
import logging
from typing import TYPE_CHECKING

import xarray as xr

from modular_dales.Geometry import GridDalesOpenBC
from modular_dales.logging_wrapper import logwrap


from modular_dales.LBC.nest_dales_in_dales.load_any_boundary_var import (
    get_boundary_dict,
    load_any_boundary_var,
)
from modular_dales.LBC.nest_dales_in_dales.timestep0 import boundaries_timestep0

if TYPE_CHECKING:
    from modular_dales.LBC.nesting_idx import NestingIndices

logger = logging.getLogger(__name__)
logger.debug("Entered module: %s", __name__)


@logwrap
def get_all_dales_boundaries(
    input_json, grid: GridDalesOpenBC, indices: "NestingIndices"
):
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
    boundary_dict = get_boundary_dict(input_json["outpath_coarse"], grid, indices)
    all_ls = []
    for boundary, (boundaryfile, sel_index) in boundary_dict.items():
        with xr.open_mfdataset(glob.glob(boundaryfile.as_posix()), join="outer") as ds:
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
