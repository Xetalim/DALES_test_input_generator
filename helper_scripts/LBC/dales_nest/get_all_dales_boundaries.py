from helper_scripts.LBC.dales_nest.timestep0 import timestep0
from helper_scripts.grids import GridDalesOpenBC


import xarray as xr


from pathlib import Path
import glob
import re
import numpy as np
import xarray as xr
from typing import Tuple, Dict, Generic


def get_all_dales_boundaries(
    input_json, grid: GridDalesOpenBC, ix_west, ix_east, iy_south, iy_north
):
    (
        ueast0,
        uwest0,
        unorth0,
        usouth0,
        utop0,
        veast0,
        vwest0,
        vnorth0,
        vsouth0,
        vtop0,
        weast0,
        wwest0,
        wnorth0,
        wsouth0,
        wtop0,
        thleast0,
        thlwest0,
        thlnorth0,
        thlsouth0,
        thltop0,
        qteast0,
        qtwest0,
        qtnorth0,
        qtsouth0,
        qttop0,
        e12east0,
        e12west0,
        e12north0,
        e12south0,
        e12top0,
        sveast0,
        svwest0,
        svnorth0,
        svsouth0,
        svtop0,
    ) = timestep0(input_json, grid, ix_west, ix_east, iy_south, iy_north)
    # Get later time steps from corresponding coarse simulation output
    # West boundary

    boundary_dict_dic = {
        "west": ["y", "z"],
        "east": ["y", "z"],
        "north": ["x", "z"],
        "south": ["x", "z"],
        "top": ["x", "y"],
    }
    var_dims_dic = {
        "u": {"x": "xm", "y": "yt", "z": "zt"},
        "v": {"x": "xt", "y": "ym", "z": "zt"},
        "w": {"x": "xt", "y": "yt", "z": "zm"},
        "default": {"x": "xt", "y": "yt", "z": "zt"},
    }

    interp_grid_dic = {
        "xt": grid.xt - grid.x0,
        "xm": grid.xm - grid.x0,
        "yt": grid.yt - grid.y0,
        "ym": grid.ym - grid.y0,
        "zt": grid.zt,
        "zm": grid.zm,
    }
    assign_grid_dic = {
        "xt": grid.xt,
        "xm": grid.xm,
        "yt": grid.yt,
        "ym": grid.ym,
        "zt": grid.zt,
        "zm": grid.zm,
    }

    interp_kwargs_dic = {"w": {"fill_value": "extrapolate"}}

    def load_var(
        ds,
        var: str,
        boundary: str,
    ):
        if var in var_dims_dic:
            var_dims = var_dims_dic[var]
        else:
            var_dims = var_dims_dic["default"]

        var_boundary_dims = [boundary_dict_dic[boundary][dim] for dim in var_dims]

        interp_coords = {dim: interp_grid_dic[dim] for dim in var_boundary_dims}
        assign_coords = {dim: assign_grid_dic[dim] for dim in var_boundary_dims}

        if var in interp_kwargs_dic:
            interp_kwargs = interp_kwargs_dic[var]
        else:
            interp_kwargs = None

        return (
            ds[var]
            .interp(interp_coords, kwargs=interp_kwargs)
            .rename(f"{var}{boundary}")
            .assign_coords(assign_coords)
        )

    with xr.open_mfdataset(
        sorted(
            list(
                Path(input_json["outpath_coarse"]).glob(f"crossyz.{ix_west+1:04d}*.nc")
            )
        ),
        join="outer",
        # chunks={"time": input_json["tchunk"]},
        combine="by_coords",
        coords="all",
        # concat_dim=["yt", "ym"],
    ) as ds:
        ds.load()
        print(ds.yt)
        uwest = (
            ds["u"]
            .interp(yt=grid.yt - grid.y0)
            .rename("uwest")
            .assign_coords(yt=grid.yt, zt=grid.zt)
        )
        # uwest = load_var(ds, "u", "west")
        uwest = xr.concat([uwest0, uwest], dim="time")
        # vwest = load_var(ds, "v", "west")
        vwest = (
            ds["v"]
            .interp(ym=grid.ym - grid.y0)
            .rename("vwest")
            .assign_coords(ym=grid.ym, zt=grid.zt)
        )
        vwest = xr.concat([vwest0, vwest], dim="time")
        wwest = (
            ds["w"]
            .interp(
                yt=grid.yt - grid.y0,
                zm=grid.zm,
                kwargs={"fill_value": "extrapolate"},
            )
            .rename("wwest")
            .assign_coords(yt=grid.yt, zm=grid.zm)
        )
        wwest = xr.concat([wwest0, wwest], dim="time")
        thlwest = (
            ds["thl"]
            .interp(yt=grid.yt - grid.y0)
            .rename("thlwest")
            .assign_coords(yt=grid.yt, zt=grid.zt)
        )
        thlwest = xr.concat([thlwest0, thlwest], dim="time")
        qtwest = (
            ds["qt"]
            .interp(yt=grid.yt - grid.y0)
            .rename("qtwest")
            .assign_coords(yt=grid.yt, zt=grid.zt)
        )
        qtwest = xr.concat([qtwest0, qtwest], dim="time")
        e12west = (
            ds["e120"]
            .interp(yt=grid.yt - grid.y0)
            .rename("e12west")
            .assign_coords(yt=grid.yt, zt=grid.zt)
        )
        e12west = xr.concat([e12west0, e12west], dim="time")
        if len(input_json["tracernames"]) > 0:
            svwest = []
            for tracername in input_json["tracernames"]:
                svwest.append(
                    ds[f"{tracername}"]
                    .interp(yt=grid.yt - grid.y0)
                    .rename(f"{tracername}west")
                    .assign_coords(yt=grid.yt, zt=grid.zt)
                )
            # for isv in range(input_json['nsv']):
            #   with xr.open_mfdataset(f"{path}sv{isv+1:03}*",chunks={"time": input_json['tchunk']}) as ds:
            #     svwest.append(ds[f"sv{isv+1:03}"].interp(yt=grid.yt-grid.y0+input_json['y_offset']).rename('svwest').assign_coords(yt=grid.yt,zt=grid.zt))
            svwest = xr.merge(svwest)
            svwest = xr.concat([svwest0, svwest], dim="time")
    # east boundary
    with xr.open_mfdataset(
        sorted(
            list(
                Path(input_json["outpath_coarse"]).glob(f"crossyz.{ix_east+1:04d}*.nc")
            )
        ),
        join="outer",
        # chunks={"time": input_json["tchunk"]},
    ) as ds:
        ueast = (
            ds["u"]
            .interp(yt=grid.yt - grid.y0)
            .rename("ueast")
            .assign_coords(yt=grid.yt, zt=grid.zt)
        )
        ueast = xr.concat([ueast0, ueast], dim="time")
        veast = (
            ds["v"]
            .interp(ym=grid.ym - grid.y0)
            .rename("veast")
            .assign_coords(ym=grid.ym, zt=grid.zt)
        )
        veast = xr.concat([veast0, veast], dim="time")
        weast = (
            ds["w"]
            .interp(
                yt=grid.yt - grid.y0,
                zm=grid.zm,
                kwargs={"fill_value": "extrapolate"},
            )
            .rename("weast")
            .assign_coords(yt=grid.yt, zm=grid.zm)
        )
        weast = xr.concat([weast0, weast], dim="time")
        thleast = (
            ds["thl"]
            .interp(yt=grid.yt - grid.y0)
            .rename("thleast")
            .assign_coords(yt=grid.yt, zt=grid.zt)
        )
        thleast = xr.concat([thleast0, thleast], dim="time")
        qteast = (
            ds["qt"]
            .interp(yt=grid.yt - grid.y0)
            .rename("qteast")
            .assign_coords(yt=grid.yt, zt=grid.zt)
        )
        qteast = xr.concat([qteast0, qteast], dim="time")
        e12east = (
            ds["e120"]
            .interp(yt=grid.yt - grid.y0)
            .rename("e12east")
            .assign_coords(yt=grid.yt, zt=grid.zt)
        )
        e12east = xr.concat([e12east0, e12east], dim="time")
        if len(input_json["tracernames"]) > 0:
            sveast = []
            for tracername in input_json["tracernames"]:
                sveast.append(
                    ds[f"{tracername}"]
                    .interp(yt=grid.yt - grid.y0)
                    .rename(f"{tracername}east")
                    .assign_coords(yt=grid.yt, zt=grid.zt)
                )
            # for isv in range(input_json['nsv']):
            #   with xr.open_mfdataset(f"{path}sv{isv+1:03}*",chunks={"time": input_json['tchunk']}) as ds:
            #     sveast.append(ds[f"sv{isv+1:03}"].interp(yt=grid.yt-grid.y0+input_json['y_offset']).rename('sveast').assign_coords(yt=grid.yt,zt=grid.zt))
            sveast = xr.merge(sveast)
            sveast = xr.concat([sveast0, sveast], dim="time")
    # south boundary
    with xr.open_mfdataset(
        sorted(
            list(
                Path(input_json["outpath_coarse"]).glob(f"crossxz.{iy_south+1:04d}*.nc")
            )
        ),
        join="outer",
        # chunks={"time": input_json["tchunk"]},
    ) as ds:
        usouth = (
            ds["u"]
            .interp(xm=grid.xm - grid.x0)
            .rename("usouth")
            .assign_coords(xm=grid.xm, zt=grid.zt)
        )
        usouth = xr.concat([usouth0, usouth], dim="time")
        vsouth = (
            ds["v"]
            .interp(xt=grid.xt - grid.x0)
            .rename("vsouth")
            .assign_coords(xt=grid.xt, zt=grid.zt)
        )
        vsouth = xr.concat([vsouth0, vsouth], dim="time")
        wsouth = (
            ds["w"]
            .interp(
                xt=grid.xt - grid.x0,
                zm=grid.zm,
                kwargs={"fill_value": "extrapolate"},
            )
            .rename("wsouth")
            .assign_coords(xt=grid.xt, zm=grid.zm)
        )
        wsouth = xr.concat([wsouth0, wsouth], dim="time")
        thlsouth = (
            ds["thl"]
            .interp(xt=grid.xt - grid.x0)
            .rename("thlsouth")
            .assign_coords(xt=grid.xt, zt=grid.zt)
        )
        thlsouth = xr.concat([thlsouth0, thlsouth], dim="time")
        qtsouth = (
            ds["qt"]
            .interp(xt=grid.xt - grid.x0)
            .rename("qtsouth")
            .assign_coords(xt=grid.xt, zt=grid.zt)
        )
        qtsouth = xr.concat([qtsouth0, qtsouth], dim="time")
        e12south = (
            ds["e120"]
            .interp(xt=grid.xt - grid.x0)
            .rename("e12south")
            .assign_coords(xt=grid.xt, zt=grid.zt)
        )
        e12south = xr.concat([e12south0, e12south], dim="time")
        if len(input_json["tracernames"]) > 0:
            svsouth = []
            for tracername in input_json["tracernames"]:
                svsouth.append(
                    ds[f"{tracername}"]
                    .interp(xt=grid.xt - grid.x0)
                    .rename(f"{tracername}south")
                    .assign_coords(xt=grid.xt, zt=grid.zt)
                )
            # for isv in range(input_json['nsv']):
            #   with xr.open_mfdataset(f"{path}sv{isv+1:03}*",chunks={"time": input_json['tchunk']}) as ds:
            #     svsouth.append(ds[f"sv{isv+1:03}"].interp(xt=grid.xt-grid.x0+input_json['x_offset']).rename('svsouth').assign_coords(xt=grid.xt,zt=grid.zt))
            svsouth = xr.merge(svsouth)
            svsouth = xr.concat([svsouth0, svsouth], dim="time")
    # north boundary
    with xr.open_mfdataset(
        sorted(
            list(
                Path(input_json["outpath_coarse"]).glob(f"crossxz.{iy_north+1:04d}*.nc")
            )
        ),
        join="outer",
        # chunks={"time": input_json["tchunk"]},
    ) as ds:
        unorth = (
            ds["u"]
            .interp(xm=grid.xm - grid.x0)
            .rename("unorth")
            .assign_coords(xm=grid.xm, zt=grid.zt)
        )
        unorth = xr.concat([unorth0, unorth], dim="time")
        vnorth = (
            ds["v"]
            .interp(xt=grid.xt - grid.x0)
            .rename("vnorth")
            .assign_coords(xt=grid.xt, zt=grid.zt)
        )
        vnorth = xr.concat([vnorth0, vnorth], dim="time")
        wnorth = (
            ds["w"]
            .interp(
                xt=grid.xt - grid.x0,
                zm=grid.zm,
                kwargs={"fill_value": "extrapolate"},
            )
            .rename("wnorth")
            .assign_coords(xt=grid.xt, zm=grid.zm)
        )
        wnorth = xr.concat([wnorth0, wnorth], dim="time")
        thlnorth = (
            ds["thl"]
            .interp(xt=grid.xt - grid.x0)
            .rename("thlnorth")
            .assign_coords(xt=grid.xt, zt=grid.zt)
        )
        thlnorth = xr.concat([thlnorth0, thlnorth], dim="time")
        qtnorth = (
            ds["qt"]
            .interp(xt=grid.xt - grid.x0)
            .rename("qtnorth")
            .assign_coords(xt=grid.xt, zt=grid.zt)
        )
        qtnorth = xr.concat([qtnorth0, qtnorth], dim="time")
        e12north = (
            ds["e120"]
            .interp(xt=grid.xt - grid.x0)
            .rename("e12north")
            .assign_coords(xt=grid.xt, zt=grid.zt)
        )
        e12north = xr.concat([e12north0, e12north], dim="time")
        if len(input_json["tracernames"]) > 0:
            svnorth = []
            for tracername in input_json["tracernames"]:
                svnorth.append(
                    ds[f"{tracername}"]
                    .interp(
                        xt=grid.xt - grid.x0,
                    )
                    .rename(f"{tracername}north")
                    .assign_coords(xt=grid.xt, zt=grid.zt)
                )
            # for isv in range(input_json['nsv']):
            #   with xr.open_mfdataset(f"{path}sv{isv+1:03}*",chunks={"time": input_json['tchunk']}) as ds:
            #     svnorth.append(ds[f"sv{isv+1:03}"].interp(xt=grid.xt-grid.x0+input_json['x_offset']).rename('svnorth').assign_coords(xt=grid.xt,zt=grid.zt))
            svnorth = xr.merge(svnorth)
            svnorth = xr.concat([svnorth0, svnorth], dim="time")
    # top boundary
    with xr.open_mfdataset(
        sorted(
            list(
                Path(input_json["outpath_coarse"]).glob(f"crossxy.{grid.kmax:04d}*.nc")
            )
        ),
        join="outer",
        # chunks={"time": input_json["tchunk"]},
    ) as ds:
        utop = (
            ds["u"]
            .interp(xm=grid.xm - grid.x0, yt=grid.yt - grid.y0)
            .rename("utop")
            .assign_coords(xm=grid.xm, yt=grid.yt)
        )
        utop = xr.concat([utop0, utop], dim="time")
        vtop = (
            ds["v"]
            .interp(xt=grid.xt - grid.x0, ym=grid.ym - grid.y0)
            .rename("vtop")
            .assign_coords(xt=grid.xt, ym=grid.ym)
        )
        vtop = xr.concat([vtop0, vtop], dim="time")
        wtop = (
            ds["w"]
            .interp(xt=grid.xt - grid.x0, yt=grid.yt - grid.y0)
            .rename("wtop")
            .assign_coords(xt=grid.xt, yt=grid.yt)
        )
        wtop = xr.concat([wtop0, wtop], dim="time")
        thltop = (
            ds["thl"]
            .interp(xt=grid.xt - grid.x0, yt=grid.yt - grid.y0)
            .rename("thltop")
            .assign_coords(xt=grid.xt, yt=grid.yt)
        )
        thltop = xr.concat([thltop0, thltop], dim="time")
        qttop = (
            ds["qt"]
            .interp(xt=grid.xt - grid.x0, yt=grid.yt - grid.y0)
            .rename("qttop")
            .assign_coords(xt=grid.xt, yt=grid.yt)
        )
        qttop = xr.concat([qttop0, qttop], dim="time")
        e12top = (
            ds["e120"]
            .interp(xt=grid.xt - grid.x0, yt=grid.yt - grid.y0)
            .rename("e12top")
            .assign_coords(xt=grid.xt, yt=grid.yt)
        )
        e12top = xr.concat([e12top0, e12top], dim="time")
        if len(input_json["tracernames"]) > 0:
            svtop = []
            for tracername in input_json["tracernames"]:
                svtop.append(
                    ds[f"{tracername}"]
                    .interp(
                        xt=grid.xt - grid.x0,
                        yt=grid.yt - grid.y0,
                    )
                    .rename(f"{tracername}top")
                    .assign_coords(xt=grid.xt, yt=grid.yt)
                )
            # for isv in range(input_json['nsv']):
            #   with xr.open_mfdataset(f"{path}sv{isv+1:03}*",chunks={"time": input_json['tchunk']}) as ds:
            #     svtop.append(ds[f"sv{isv+1:03}"].interp(xt=grid.xt-grid.x0+input_json['x_offset'],yt=grid.yt+input_json['y_offset']).rename('svtop').assign_coords(xt=grid.xt,yt=grid.yt))
            svtop = xr.merge(svtop)
            svtop = xr.concat([svtop0, svtop], dim="time")
    if len(input_json["tracernames"]) > 0:
        return (
            uwest,
            vwest,
            wwest,
            thlwest,
            qtwest,
            e12west,
            svwest,
            ueast,
            veast,
            weast,
            thleast,
            qteast,
            e12east,
            sveast,
            usouth,
            vsouth,
            wsouth,
            thlsouth,
            qtsouth,
            e12south,
            svsouth,
            unorth,
            vnorth,
            wnorth,
            thlnorth,
            qtnorth,
            e12north,
            svnorth,
            utop,
            vtop,
            wtop,
            thltop,
            qttop,
            e12top,
            svtop,
        )
    else:
        return (
            uwest,
            vwest,
            wwest,
            thlwest,
            qtwest,
            e12west,
            [],
            ueast,
            veast,
            weast,
            thleast,
            qteast,
            e12east,
            [],
            usouth,
            vsouth,
            wsouth,
            thlsouth,
            qtsouth,
            e12south,
            [],
            unorth,
            vnorth,
            wnorth,
            thlnorth,
            qtnorth,
            e12north,
            [],
            utop,
            vtop,
            wtop,
            thltop,
            qttop,
            e12top,
            [],
        )
