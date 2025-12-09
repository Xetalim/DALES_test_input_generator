import pandas as pd
from helper_scripts.grids import GridDalesOpenBC


def get_u0(input_json, grid: GridDalesOpenBC, ix_west, ix_east, iy_south, iy_north, ds):
    uwest0 = (
        ds["u0"]
        .isel(xm=ix_west, drop=True)
        .interp(yt=grid.yt)
        .rename("uwest")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(yt=grid.yt, zt=grid.zt)
    )
    ueast0 = (
        ds["u0"]
        .isel(xm=ix_east, drop=True)
        .interp(yt=grid.yt)
        .rename("ueast")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(yt=grid.yt, zt=grid.zt)
    )
    usouth0 = (
        ds["u0"]
        .isel(yt=iy_south, drop=True)
        .interp(xm=grid.xm)
        .rename("usouth")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(xm=grid.xm, zt=grid.zt)
    )
    unorth0 = (
        ds["u0"]
        .isel(yt=iy_north, drop=True)
        .interp(xm=grid.xm)
        .rename("unorth")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(xm=grid.xm, zt=grid.zt)
    )
    utop0 = (
        ds["u0"]
        .isel(zt=grid.kmax - 1, drop=True)
        .interp(
            xm=grid.xm,
            yt=grid.yt,
        )
        .rename("utop")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(xm=grid.xm, yt=grid.yt)
    )
    return ueast0, uwest0, unorth0, usouth0, utop0


def get_v0(input_json, grid, ix_west, ix_east, iy_south, iy_north, ds):
    veast0 = (
        ds["v0"]
        .isel(xt=ix_east, drop=True)
        .interp(ym=grid.ym)
        .rename("veast")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(ym=grid.ym, zt=grid.zt)
    )
    vwest0 = (
        ds["v0"]
        .isel(xt=ix_west, drop=True)
        .interp(ym=grid.ym)
        .rename("vwest")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(ym=grid.ym, zt=grid.zt)
    )
    vnorth0 = (
        ds["v0"]
        .isel(ym=iy_north, drop=True)
        .interp(xt=grid.xt)
        .rename("vnorth")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(xt=grid.xt, zt=grid.zt)
    )
    vsouth0 = (
        ds["v0"]
        .isel(ym=iy_south, drop=True)
        .interp(xt=grid.xt)
        .rename("vsouth")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(xt=grid.xt, zt=grid.zt)
    )
    vtop0 = (
        ds["v0"]
        .isel(zt=grid.kmax - 1, drop=True)
        .interp(
            xt=grid.xt,
            ym=grid.ym,
        )
        .rename("vtop")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(xt=grid.xt, ym=grid.ym)
    )
    return veast0, vwest0, vnorth0, vsouth0, vtop0


def get_w0(input_json, grid: GridDalesOpenBC, ix_west, ix_east, iy_south, iy_north, ds):
    weast0 = (
        ds["w0"]
        .isel(xt=ix_east, drop=True)
        .interp(
            yt=grid.yt,
            zm=grid.zm,
            kwargs={"fill_value": "extrapolate"},
        )
        .rename("weast")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(yt=grid.yt, zm=grid.zm)
    )
    wwest0 = (
        ds["w0"]
        .isel(xt=ix_west, drop=True)
        .interp(
            yt=grid.yt,
            zm=grid.zm,
            kwargs={"fill_value": "extrapolate"},
        )
        .rename("wwest")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(yt=grid.yt, zm=grid.zm)
    )
    wnorth0 = (
        ds["w0"]
        .isel(yt=iy_north, drop=True)
        .interp(
            xt=grid.xt,
            zm=grid.zm,
            kwargs={"fill_value": "extrapolate"},
        )
        .rename("wnorth")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(xt=grid.xt, zm=grid.zm)
    )
    wsouth0 = (
        ds["w0"]
        .isel(yt=iy_south, drop=True)
        .interp(
            xt=grid.xt,
            zm=grid.zm,
            kwargs={"fill_value": "extrapolate"},
        )
        .rename("wsouth")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(xt=grid.xt, zm=grid.zm)
    )
    wtop0 = (
        ds["w0"]
        .isel(zm=grid.kmax, drop=True)
        .interp(
            xt=grid.xt,
            yt=grid.yt,
        )
        .rename("wtop")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(xt=grid.xt, yt=grid.yt)
    )
    return weast0, wwest0, wnorth0, wsouth0, wtop0


def get_thl0(input_json, grid, ix_west, ix_east, iy_south, iy_north, ds):
    thleast0 = (
        ds["thl0"]
        .isel(xt=ix_east, drop=True)
        .interp(yt=grid.yt)
        .rename("thleast")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(yt=grid.yt, zt=grid.zt)
    )
    thlwest0 = (
        ds["thl0"]
        .isel(xt=ix_west, drop=True)
        .interp(yt=grid.yt)
        .rename("thlwest")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(yt=grid.yt, zt=grid.zt)
    )
    thlnorth0 = (
        ds["thl0"]
        .isel(yt=iy_north, drop=True)
        .interp(xt=grid.xt)
        .rename("thlnorth")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(xt=grid.xt, zt=grid.zt)
    )
    thlsouth0 = (
        ds["thl0"]
        .isel(yt=iy_south, drop=True)
        .interp(xt=grid.xt)
        .rename("thlsouth")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(xt=grid.xt, zt=grid.zt)
    )
    thltop0 = (
        ds["thl0"]
        .isel(zt=grid.kmax - 1, drop=True)
        .interp(
            xt=grid.xt,
            yt=grid.yt,
        )
        .rename("thltop")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(xt=grid.xt, yt=grid.yt)
    )
    return thleast0, thlwest0, thlnorth0, thlsouth0, thltop0


def get_qt0(
    input_json, grid: GridDalesOpenBC, ix_west, ix_east, iy_south, iy_north, ds
):
    qteast0 = (
        ds["qt0"]
        .isel(xt=ix_east, drop=True)
        .interp(yt=grid.yt)
        .rename("qteast")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(yt=grid.yt, zt=grid.zt)
    )
    qtwest0 = (
        ds["qt0"]
        .isel(xt=ix_west, drop=True)
        .interp(yt=grid.yt)
        .rename("qtwest")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(yt=grid.yt, zt=grid.zt)
    )
    qtnorth0 = (
        ds["qt0"]
        .isel(yt=iy_north, drop=True)
        .interp(xt=grid.xt)
        .rename("qtnorth")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(xt=grid.xt, zt=grid.zt)
    )
    qtsouth0 = (
        ds["qt0"]
        .isel(yt=iy_south, drop=True)
        .interp(xt=grid.xt)
        .rename("qtsouth")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(xt=grid.xt, zt=grid.zt)
    )
    qttop0 = (
        ds["qt0"]
        .isel(zt=grid.kmax - 1, drop=True)
        .interp(
            xt=grid.xt,
            yt=grid.yt,
        )
        .rename("qttop")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(xt=grid.xt, yt=grid.yt)
    )
    return qteast0, qtwest0, qtnorth0, qtsouth0, qttop0


def get_e120(
    input_json, grid: GridDalesOpenBC, ix_west, ix_east, iy_south, iy_north, ds
):
    e12east0 = (
        ds["e120"]
        .isel(xt=ix_east, drop=True)
        .interp(yt=grid.yt)
        .rename("e12east")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(yt=grid.yt, zt=grid.zt)
    )
    e12west0 = (
        ds["e120"]
        .isel(xt=ix_west, drop=True)
        .interp(yt=grid.yt)
        .rename("e12west")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(yt=grid.yt, zt=grid.zt)
    )
    e12north0 = (
        ds["e120"]
        .isel(yt=iy_north, drop=True)
        .interp(xt=grid.xt)
        .rename("e12north")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(xt=grid.xt, zt=grid.zt)
    )
    e12south0 = (
        ds["e120"]
        .isel(yt=iy_south, drop=True)
        .interp(xt=grid.xt)
        .rename("e12south")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(xt=grid.xt, zt=grid.zt)
    )
    e12top0 = (
        ds["e120"]
        .isel(zt=grid.kmax - 1, drop=True)
        .interp(
            xt=grid.xt,
            yt=grid.yt,
        )
        .rename("e12top")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(xt=grid.xt, yt=grid.yt)
    )

    return e12east0, e12west0, e12north0, e12south0, e12top0
