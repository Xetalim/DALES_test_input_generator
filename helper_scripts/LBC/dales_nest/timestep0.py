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


import pandas as pd
import xarray as xr


from pathlib import Path


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
            "west": sorted(
                list(
                    Path(input_json["outpath_coarse_old"]).glob(
                        f"crossyz.{indices.ix_west+1:04d}*.nc"
                    )
                )
            ),
            "east": sorted(
                list(
                    Path(input_json["outpath_coarse_old"]).glob(
                        f"crossyz.{indices.ix_east+1:04d}*.nc"
                    )
                )
            ),
            "south": sorted(
                list(
                    Path(input_json["outpath_coarse_old"]).glob(
                        f"crossxz.{indices.iy_south+1:04d}*.nc"
                    )
                )
            ),
            "north": sorted(
                list(
                    Path(input_json["outpath_coarse_old"]).glob(
                        f"crossxz.{indices.iy_north+1:04d}*.nc"
                    )
                )
            ),
            "top": sorted(
                list(
                    Path(input_json["outpath_coarse_old"]).glob(
                        f"crossxy.{grid.kmax:04d}*.nc"
                    )
                )
            ),
        }
        all_ls = []
        for boundary, boundaryfiles in boundary_dict.items():
            with xr.open_mfdataset(
                boundaryfiles,
                join="outer",
                # chunks={"time": input_json["tchunk"]},
            ) as ds:
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
                            ds,
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
