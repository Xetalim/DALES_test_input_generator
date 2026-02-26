from typing import TYPE_CHECKING
import logging

import xarray as xr


from modular_dales.Geometry.GridDales import GridDalesOpenBC
from modular_dales.LBC.nest_dales_in_dales.load_any_boundary_var import (
    load_any_boundary_var,
)
from modular_dales.logging_wrapper import logwrap
from modular_dales.LBC.nest_dales_in_dales.load_any_boundary_var import (
    get_boundary_dict,
)

if TYPE_CHECKING:
    from modular_dales.LBC.nesting_idx import NestingIndices

logger = logging.getLogger(__name__)
logger.debug("Entered module: %s", __name__)


@logwrap
def boundaries_timestep0(input_json, grid: GridDalesOpenBC, indices: "NestingIndices"):
    """
    Docstring for timestep0

    :param input_json: Configuration dict for LBC
    :param grid: Grid object describing the output grid for open boundaries, which is different from the normal DALES grid.
    :type grid: GridDalesOpenBC
    :param indices: Indices describing where in the supergrid the edges of the current grid lie.
    :type indices: nesting_idx

    This function returns the boundary fields for the initial timestep.
    As DALES currently doesn't output crossections at t=0, we must get the boundary information from initfields.nc if we start our
    simulation at the same time as the host simulation.
    """
    if input_json["time0"] == input_json["start"]:
        all_ls = []
        with xr.open_mfdataset(
            f"{input_json['inpath_coarse']}initfields.inp.*.nc"
        ) as ds:
            for boundary in ["west", "east", "north", "south", "top"]:
                for var in ["u", "v", "w", "thl", "qt", "e12"]:
                    all_ls.append(
                        load_any_boundary_var(
                            ds,
                            var,
                            boundary,
                            grid,
                            indices,
                            isel=True,
                            expand_dims=True,
                            expand_dims_time0=input_json["time0"],
                            var_postfix="0",
                            move_x0y0=False,
                        )
                    )
                if len(input_json["tracernames"]) > 0:
                    sv_boundary = []
                    for tracername in input_json["tracernames"]:
                        sv_boundary.append(
                            xr.zeros_like(all_ls[-1]).rename(f"{tracername}{boundary}")
                        )
                    all_ls.append(xr.merge(sv_boundary))

    # Get initial boundary fields from previous simulation, specifically,
    # the last time step in the output of the previous simulation
    else:

        boundary_dict = get_boundary_dict(
            input_json["outpath_coarse_old"], grid, indices
        )

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
                    sel_index = sel_index
                    all_ls.append(
                        load_any_boundary_var(
                            ds.sel(sel_index),
                            var,
                            boundary=boundary,
                            grid=grid,
                            indices=indices,
                            isel={"time": -1},
                            expand_dims=True,
                            expand_dims_time0=input_json["start"],
                            var_postfix=var_postfix,
                        )
                    )

    return xr.merge(all_ls)
