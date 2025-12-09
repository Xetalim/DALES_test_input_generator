from helper_scripts.LBC.dales_nest.timestep0 import timestep0
from helper_scripts.grids import GridDalesOpenBC, nesting_idx
from helper_scripts.LBC.dales_nest.get_timesteps0 import load_var

import xarray as xr


from pathlib import Path
import glob
import re
import numpy as np
import xarray as xr
from typing import Tuple, Dict, Generic


def get_all_dales_boundaries(input_json, grid: GridDalesOpenBC, indices: nesting_idx):
    ds_timestep0 = timestep0(input_json, grid, indices)
    # Get later time steps from corresponding coarse simulation output
    # West boundary
    boundary_dict = {
        "west": sorted(
            list(
                Path(input_json["outpath_coarse"]).glob(
                    f"crossyz.{indices.ix_west+1:04d}*.nc"
                )
            )
        ),
        "east": sorted(
            list(
                Path(input_json["outpath_coarse"]).glob(
                    f"crossyz.{indices.ix_east+1:04d}*.nc"
                )
            )
        ),
        "south": sorted(
            list(
                Path(input_json["outpath_coarse"]).glob(
                    f"crossxz.{indices.iy_south+1:04d}*.nc"
                )
            )
        ),
        "north": sorted(
            list(
                Path(input_json["outpath_coarse"]).glob(
                    f"crossxz.{indices.iy_north+1:04d}*.nc"
                )
            )
        ),
        "top": sorted(
            list(
                Path(input_json["outpath_coarse"]).glob(f"crossxy.{grid.kmax:04d}*.nc")
            )
        ),
    }
    all_ls = []
    for boundary, (boundaryfiles, sel_index) in boundary_dict.items():
        with xr.open_mfdataset(
            boundaryfiles,
            join="outer",
            # chunks={"time": input_json["tchunk"]},
            combine="by_coords",
            coords="all",
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
                        var_postfix=var_postfix,
                    )
                )
    return xr.concat([ds_timestep0, xr.merge(all_ls)], dim="time")
