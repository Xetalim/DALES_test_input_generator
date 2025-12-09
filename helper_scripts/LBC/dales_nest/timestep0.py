from helper_scripts.LBC.dales_nest.get_timesteps0 import (
    get_e120,
    get_qt0,
    get_thl0,
    get_u0,
    get_v0,
    get_w0,
)
from helper_scripts.grids import GridDalesOpenBC


import pandas as pd
import xarray as xr


from pathlib import Path


def timestep0(input_json, grid: GridDalesOpenBC, ix_west, ix_east, iy_south, iy_north):
    if input_json["time0"] == input_json["start"]:
        with xr.open_mfdataset(
            f"{input_json['inpath_coarse']}initfields.inp.*.nc"
        ) as ds:
            # West boundary
            ueast0, uwest0, unorth0, usouth0, utop0 = get_u0(
                input_json, grid, ix_west, ix_east, iy_south, iy_north, ds
            )
            veast0, vwest0, vnorth0, vsouth0, vtop0 = get_v0(
                input_json, grid, ix_west, ix_east, iy_south, iy_north, ds
            )
            weast0, wwest0, wnorth0, wsouth0, wtop0 = get_w0(
                input_json, grid, ix_west, ix_east, iy_south, iy_north, ds
            )
            thleast0, thlwest0, thlnorth0, thlsouth0, thltop0 = get_thl0(
                input_json, grid, ix_west, ix_east, iy_south, iy_north, ds
            )
            qteast0, qtwest0, qtnorth0, qtsouth0, qttop0 = get_qt0(
                input_json, grid, ix_west, ix_east, iy_south, iy_north, ds
            )
            e12east0, e12west0, e12north0, e12south0, e12top0 = get_e120(
                input_json, grid, ix_west, ix_east, iy_south, iy_north, ds
            )
            if len(input_json["tracernames"]) > 0:
                sveast0, svwest0, svnorth0, svsouth0, svtop0 = [], [], [], [], []
                for tracername in input_json["tracernames"]:
                    svwest0.append(xr.zeros_like(e12west0).rename(f"{tracername}west"))
                    sveast0.append(xr.zeros_like(e12east0).rename(f"{tracername}east"))
                    svsouth0.append(
                        xr.zeros_like(e12south0).rename(f"{tracername}south")
                    )
                    svnorth0.append(
                        xr.zeros_like(e12north0).rename(f"{tracername}north")
                    )
                    svtop0.append(xr.zeros_like(e12top0).rename(f"{tracername}top"))
                svwest0 = xr.merge(svwest0)
                sveast0 = xr.merge(sveast0)
                svsouth0 = xr.merge(svsouth0)
                svnorth0 = xr.merge(svnorth0)
                svtop0 = xr.merge(svtop0)
    # Get initial boundary fields from previous simulation
    else:
        # West boundary
        with xr.open_mfdataset(
            sorted(
                list(
                    Path(input_json["outpath_coarse_old"]).glob(
                        f"crossyz.{ix_west+1:04d}*.nc"
                    )
                )
            ),
            join="outer",
            # chunks={"time": input_json["tchunk"]},
        ) as ds:
            uwest0 = (
                ds["u"]
                .isel(time=-1, drop=True)
                .interp(yt=grid.yt - grid.y0)
                .rename("uwest")
                .assign_coords(yt=grid.yt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            vwest0 = (
                ds["v"]
                .isel(time=-1, drop=True)
                .interp(ym=grid.ym - grid.y0)
                .rename("vwest")
                .assign_coords(ym=grid.ym, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            wwest0 = (
                ds["w"]
                .isel(time=-1, drop=True)
                .interp(
                    yt=grid.yt,
                    zm=grid.zm,
                    kwargs={"fill_value": "extrapolate"},
                )
                .rename("wwest")
                .assign_coords(yt=grid.yt, zm=grid.zm)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            thlwest0 = (
                ds["thl"]
                .isel(time=-1, drop=True)
                .interp(yt=grid.yt - grid.y0)
                .rename("thlwest")
                .assign_coords(yt=grid.yt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            qtwest0 = (
                ds["qt"]
                .isel(time=-1, drop=True)
                .interp(yt=grid.yt - grid.y0)
                .rename("qtwest")
                .assign_coords(yt=grid.yt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            e12west0 = (
                ds["e120"]
                .isel(time=-1, drop=True)
                .interp(yt=grid.yt - grid.y0)
                .rename("e12west")
                .assign_coords(yt=grid.yt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            if len(input_json["tracernames"]) > 0:
                svwest0 = []
                for tracername in input_json["tracernames"]:
                    svwest0.append(
                        ds[f"{tracername}"]
                        .isel(time=-1, drop=True)
                        .interp(
                            yt=grid.yt,
                        )
                        .rename(f"{tracername}west")
                        .assign_coords(yt=grid.yt, zt=grid.zt)
                        .expand_dims(
                            {"time": [pd.Timestamp(input_json["start"])]}, axis=0
                        )
                    )
                svwest0 = xr.merge(svwest0)
        # east boundary
        with xr.open_mfdataset(
            sorted(
                list(
                    Path(input_json["outpath_coarse_old"]).glob(
                        f"crossyz.{ix_east+1:04d}*.nc"
                    )
                )
            ),
            join="outer",
            # chunks={"time": input_json["tchunk"]},
        ) as ds:
            ueast0 = (
                ds["u"]
                .isel(time=-1, drop=True)
                .interp(yt=grid.yt - grid.y0)
                .rename("ueast")
                .assign_coords(yt=grid.yt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            veast0 = (
                ds["v"]
                .isel(time=-1, drop=True)
                .interp(ym=grid.ym - grid.y0)
                .rename("veast")
                .assign_coords(ym=grid.ym, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            weast0 = (
                ds["w"]
                .isel(time=-1, drop=True)
                .interp(
                    yt=grid.yt,
                    zm=grid.zm,
                    kwargs={"fill_value": "extrapolate"},
                )
                .rename("weast")
                .assign_coords(yt=grid.yt, zm=grid.zm)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            thleast0 = (
                ds["thl"]
                .isel(time=-1, drop=True)
                .interp(yt=grid.yt - grid.y0)
                .rename("thleast")
                .assign_coords(yt=grid.yt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            qteast0 = (
                ds["qt"]
                .isel(time=-1, drop=True)
                .interp(yt=grid.yt - grid.y0)
                .rename("qteast")
                .assign_coords(yt=grid.yt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            e12east0 = (
                ds["e120"]
                .isel(time=-1, drop=True)
                .interp(yt=grid.yt - grid.y0)
                .rename("e12east")
                .assign_coords(yt=grid.yt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            if len(input_json["tracernames"]) > 0:
                sveast0 = []
                for tracername in input_json["tracernames"]:
                    sveast0.append(
                        ds[f"{tracername}"]
                        .isel(time=-1, drop=True)
                        .interp(
                            yt=grid.yt,
                        )
                        .rename(f"{tracername}east")
                        .assign_coords(yt=grid.yt, zt=grid.zt)
                        .expand_dims(
                            {"time": [pd.Timestamp(input_json["start"])]}, axis=0
                        )
                    )
            sveast0 = xr.merge(sveast0)
        # south boundary
        with xr.open_mfdataset(
            sorted(
                list(
                    Path(input_json["outpath_coarse_old"]).glob(
                        f"crossxz.{iy_south+1:04d}*.nc"
                    )
                )
            ),
            join="outer",
            # chunks={"time": input_json["tchunk"]},
        ) as ds:
            usouth0 = (
                ds["u"]
                .isel(time=-1, drop=True)
                .interp(xm=grid.xm - grid.x0)
                .rename("usouth")
                .assign_coords(xm=grid.xm, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            vsouth0 = (
                ds["v"]
                .isel(time=-1, drop=True)
                .interp(xt=grid.xt - grid.x0)
                .rename("vsouth")
                .assign_coords(xt=grid.xt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            wsouth0 = (
                ds["w"]
                .isel(time=-1, drop=True)
                .interp(
                    xt=grid.xt,
                    zm=grid.zm,
                    kwargs={"fill_value": "extrapolate"},
                )
                .rename("wsouth")
                .assign_coords(xt=grid.xt, zm=grid.zm)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            thlsouth0 = (
                ds["thl"]
                .isel(time=-1, drop=True)
                .interp(xt=grid.xt - grid.x0)
                .rename("thlsouth")
                .assign_coords(xt=grid.xt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            qtsouth0 = (
                ds["qt"]
                .isel(time=-1, drop=True)
                .interp(xt=grid.xt - grid.x0)
                .rename("qtsouth")
                .assign_coords(xt=grid.xt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            e12south0 = (
                ds["e120"]
                .isel(time=-1, drop=True)
                .interp(xt=grid.xt - grid.x0)
                .rename("e12south")
                .assign_coords(xt=grid.xt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            if len(input_json["tracernames"]) > 0:
                svsouth0 = []
                for tracername in input_json["tracernames"]:
                    svsouth0.append(
                        ds[f"{tracername}"]
                        .isel(time=-1, drop=True)
                        .interp(
                            xt=grid.xt,
                        )
                        .rename(f"{tracername}south")
                        .assign_coords(xt=grid.xt, zt=grid.zt)
                        .expand_dims(
                            {"time": [pd.Timestamp(input_json["start"])]}, axis=0
                        )
                    )
                svsouth0 = xr.merge(svsouth0)
        # north boundary
        with xr.open_mfdataset(
            sorted(
                list(
                    Path(input_json["outpath_coarse_old"]).glob(
                        f"crossxz.{iy_north+1:04d}*.nc"
                    )
                )
            ),
            join="outer",
            # chunks={"time": input_json["tchunk"]},
        ) as ds:
            unorth0 = (
                ds["u"]
                .isel(time=-1, drop=True)
                .interp(xm=grid.xm - grid.x0)
                .rename("unorth")
                .assign_coords(xm=grid.xm, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            vnorth0 = (
                ds["v"]
                .isel(time=-1, drop=True)
                .interp(xt=grid.xt - grid.x0)
                .rename("vnorth")
                .assign_coords(xt=grid.xt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            wnorth0 = (
                ds["w"]
                .isel(time=-1, drop=True)
                .interp(
                    xt=grid.xt,
                    zm=grid.zm,
                    kwargs={"fill_value": "extrapolate"},
                )
                .rename("wnorth")
                .assign_coords(xt=grid.xt, zm=grid.zm)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            thlnorth0 = (
                ds["thl"]
                .isel(time=-1, drop=True)
                .interp(xt=grid.xt - grid.x0)
                .rename("thlnorth")
                .assign_coords(xt=grid.xt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            qtnorth0 = (
                ds["qt"]
                .isel(time=-1, drop=True)
                .interp(xt=grid.xt - grid.x0)
                .rename("qtnorth")
                .assign_coords(xt=grid.xt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            e12north0 = (
                ds["e120"]
                .isel(time=-1, drop=True)
                .interp(xt=grid.xt - grid.x0)
                .rename("e12north")
                .assign_coords(xt=grid.xt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            if len(input_json["tracernames"]) > 0:
                svnorth0 = []
                for tracername in input_json["tracernames"]:
                    svnorth0.append(
                        ds[f"{tracername}"]
                        .isel(time=-1, drop=True)
                        .interp(
                            xt=grid.xt,
                        )
                        .rename(f"{tracername}north")
                        .assign_coords(xt=grid.xt, zt=grid.zt)
                        .expand_dims(
                            {"time": [pd.Timestamp(input_json["start"])]}, axis=0
                        )
                    )
                svnorth0 = xr.merge(svnorth0)
        # top boundary
        with xr.open_mfdataset(
            sorted(
                list(
                    Path(input_json["outpath_coarse_old"]).glob(
                        f"crossxy.{grid.kmax:04d}*.nc"
                    )
                )
            ),
            join="outer",
            # chunks={"time": input_json["tchunk"]},
        ) as ds:
            utop0 = (
                ds["u"]
                .isel(time=-1, drop=True)
                .interp(
                    xm=grid.xm,
                    yt=grid.yt,
                )
                .rename("utop")
                .assign_coords(xm=grid.xm, yt=grid.yt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            vtop0 = (
                ds["v"]
                .isel(time=-1, drop=True)
                .interp(
                    xt=grid.xt,
                    ym=grid.ym,
                )
                .rename("vtop")
                .assign_coords(xt=grid.xt, ym=grid.ym)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            wtop0 = (
                ds["w"]
                .isel(time=-1, drop=True)
                .interp(
                    xt=grid.xt,
                    yt=grid.yt,
                )
                .rename("wtop")
                .assign_coords(xt=grid.xt, yt=grid.yt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            thltop0 = (
                ds["thl"]
                .isel(time=-1, drop=True)
                .interp(
                    xt=grid.xt,
                    yt=grid.yt,
                )
                .rename("thltop")
                .assign_coords(xt=grid.xt, yt=grid.yt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            qttop0 = (
                ds["qt"]
                .isel(time=-1, drop=True)
                .interp(
                    xt=grid.xt,
                    yt=grid.yt,
                )
                .rename("qttop")
                .assign_coords(xt=grid.xt, yt=grid.yt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            e12top0 = (
                ds["e120"]
                .isel(time=-1, drop=True)
                .interp(
                    xt=grid.xt,
                    yt=grid.yt,
                )
                .rename("e12top")
                .assign_coords(xt=grid.xt, yt=grid.yt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            if len(input_json["tracernames"]) > 0:
                svtop0 = []
                for tracername in input_json["tracernames"]:
                    svtop0.append(
                        ds[f"{tracername}"]
                        .isel(time=-1, drop=True)
                        .interp(
                            xt=grid.xt,
                            yt=grid.yt,
                        )
                        .rename(f"{tracername}top")
                        .assign_coords(xt=grid.xt, yt=grid.yt)
                        .expand_dims(
                            {"time": [pd.Timestamp(input_json["start"])]}, axis=0
                        )
                    )
                svtop0 = xr.merge(svtop0)
    # return (
    #     ueast0,
    #     uwest0,
    #     unorth0,
    #     usouth0,
    #     utop0,
    #     veast0,
    #     vwest0,
    #     vnorth0,
    #     vsouth0,
    #     vtop0,
    #     weast0,
    #     wwest0,
    #     wnorth0,
    #     wsouth0,
    #     wtop0,
    #     thleast0,
    #     thlwest0,
    #     thlnorth0,
    #     thlsouth0,
    #     thltop0,
    #     qteast0,
    #     qtwest0,
    #     qtnorth0,
    #     qtsouth0,
    #     qttop0,
    #     e12east0,
    #     e12west0,
    #     e12north0,
    #     e12south0,
    #     e12top0,
    #     sveast0,
    #     svwest0,
    #     svnorth0,
    #     svsouth0,
    #     svtop0,
    # )

    if len(input_json["tracernames"]) > 0:
        return (
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
        )
    else:
        return (
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
            [],
            [],
            [],
            [],
            [],
        )
