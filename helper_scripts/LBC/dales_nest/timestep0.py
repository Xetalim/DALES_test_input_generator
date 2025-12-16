from helper_scripts.LBC.dales_nest.get_timesteps0 import (
    get_e120,
    get_qt0,
    get_thl0,
    get_u0,
    get_v0,
    get_w0,
    load_var,
)
from helper_scripts.grids import GridDalesOpenBC, nesting_idx
from helper_scripts.logging_wrapper import logwrap

import pandas as pd
import xarray as xr

import logging

logger = logging.getLogger(__name__)
logger.debug("Entered module: %s", __name__)


from pathlib import Path


@logwrap
def timestep0(input_json, grid: GridDalesOpenBC, indices: nesting_idx):
    if input_json["time0"] == input_json["start"]:
        all_ls = []
        with xr.open_mfdataset(
            f"{input_json['inpath_coarse']}initfields.inp.*.nc"
        ) as ds:
            for boundary in ["west", "east", "north", "south", "top"]:
                for var in ["u", "v", "w", "thl", "qt", "e12"]:
                    all_ls.append(
                        load_var(
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

    # Get initial boundary fields from previous simulation
    else:
        boundary_dict = {
            "west": (
                Path(input_json["outpath_coarse_old"]) / ".." / "crossyz.nc",
                {"xm": indices.indices + 1, "xt": indices.indices + 1},
            ),
            "east": (
                Path(input_json["outpath_coarse_old"]) / ".." / "crossyz.nc",
                {"xm": indices.ix_east + 1, "xt": indices.ix_east + 1},
            ),
            "south": (
                Path(input_json["outpath_coarse_old"]) / ".." / "crossxz.nc",
                {"ym": indices.iy_south + 1, "yt": indices.iy_south + 1},
            ),
            "north": (
                Path(input_json["outpath_coarse_old"]) / ".." / "crossxz.nc",
                {"ym": indices.iy_north + 1, "yt": indices.iy_north + 1},
            ),
            "top": (
                Path(input_json["outpath_coarse_old"]) / ".." / "crossxy.nc",
                {"zt": grid.kmax, "zm": grid.kmax},
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
                        load_var(
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
