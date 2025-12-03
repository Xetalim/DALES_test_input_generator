import pandas as pd
import xarray as xr
from pathlib import Path


def get_u0(input_json, grid, ix_west, ix_east, iy_south, iy_north, ds):
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


def get_w0(input_json, grid, ix_west, ix_east, iy_south, iy_north, ds):
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


def get_qt0(input_json, grid, ix_west, ix_east, iy_south, iy_north, ds):
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


def get_e120(input_json, grid, ix_west, ix_east, iy_south, iy_north, ds):
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


def get_all_dales_boundaries(input_json, grid, ix_west, ix_east, iy_south, iy_north):
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
    path = Path(input_json["outpath_coarse"]) / f"crossyz.{ix_west+2:04d}"
    with xr.open_mfdataset(path.glob("*"), chunks={"time": input_json["tchunk"]}) as ds:
        uwest = (
            ds["uyz"]
            .interp(yt=grid.yt)
            .rename("uwest")
            .assign_coords(yt=grid.yt, zt=grid.zt)
        )
        uwest = xr.concat([uwest0, uwest], dim="time")
        vwest = (
            ds["vyz"]
            .interp(ym=grid.ym)
            .rename("vwest")
            .assign_coords(ym=grid.ym, zt=grid.zt)
        )
        vwest = xr.concat([vwest0, vwest], dim="time")
        wwest = (
            ds["wyz"]
            .interp(
                yt=grid.yt,
                zm=grid.zm,
                kwargs={"fill_value": "extrapolate"},
            )
            .rename("wwest")
            .assign_coords(yt=grid.yt, zm=grid.zm)
        )
        wwest = xr.concat([wwest0, wwest], dim="time")
        thlwest = (
            ds["thlyz"]
            .interp(yt=grid.yt)
            .rename("thlwest")
            .assign_coords(yt=grid.yt, zt=grid.zt)
        )
        thlwest = xr.concat([thlwest0, thlwest], dim="time")
        qtwest = (
            ds["qtyz"]
            .interp(yt=grid.yt)
            .rename("qtwest")
            .assign_coords(yt=grid.yt, zt=grid.zt)
        )
        qtwest = xr.concat([qtwest0, qtwest], dim="time")
        e12west = (
            ds["e120yz"]
            .interp(yt=grid.yt)
            .rename("e12west")
            .assign_coords(yt=grid.yt, zt=grid.zt)
        )
        e12west = xr.concat([e12west0, e12west], dim="time")
        if len(input_json["tracernames"]) > 0:
            svwest = []
            for tracername in input_json["tracernames"]:
                svwest.append(
                    ds[f"{tracername}yz"]
                    .interp(yt=grid.yt)
                    .rename(f"{tracername}west")
                    .assign_coords(yt=grid.yt, zt=grid.zt)
                )
            # for isv in range(input_json['nsv']):
            #   with xr.open_mfdataset(f"{path}sv{isv+1:03}*",chunks={"time": input_json['tchunk']}) as ds:
            #     svwest.append(ds[f"sv{isv+1:03}"].interp(yt=grid.yt+input_json['y_offset']).rename('svwest').assign_coords(yt=grid.yt,zt=grid.zt))
            svwest = xr.merge(svwest)
            svwest = xr.concat([svwest0, svwest], dim="time")
    # east boundary
    path = Path(input_json["outpath_coarse"]) / f"crossyz.{ix_east+2:04d}"
    with xr.open_mfdataset(path.glob("*"), chunks={"time": input_json["tchunk"]}) as ds:
        ueast = (
            ds["uyz"]
            .interp(yt=grid.yt)
            .rename("ueast")
            .assign_coords(yt=grid.yt, zt=grid.zt)
        )
        ueast = xr.concat([ueast0, ueast], dim="time")
        veast = (
            ds["vyz"]
            .interp(ym=grid.ym)
            .rename("veast")
            .assign_coords(ym=grid.ym, zt=grid.zt)
        )
        veast = xr.concat([veast0, veast], dim="time")
        weast = (
            ds["wyz"]
            .interp(
                yt=grid.yt,
                zm=grid.zm,
                kwargs={"fill_value": "extrapolate"},
            )
            .rename("weast")
            .assign_coords(yt=grid.yt, zm=grid.zm)
        )
        weast = xr.concat([weast0, weast], dim="time")
        thleast = (
            ds["thlyz"]
            .interp(yt=grid.yt)
            .rename("thleast")
            .assign_coords(yt=grid.yt, zt=grid.zt)
        )
        thleast = xr.concat([thleast0, thleast], dim="time")
        qteast = (
            ds["qtyz"]
            .interp(yt=grid.yt)
            .rename("qteast")
            .assign_coords(yt=grid.yt, zt=grid.zt)
        )
        qteast = xr.concat([qteast0, qteast], dim="time")
        e12east = (
            ds["e120yz"]
            .interp(yt=grid.yt)
            .rename("e12east")
            .assign_coords(yt=grid.yt, zt=grid.zt)
        )
        e12east = xr.concat([e12east0, e12east], dim="time")
        if len(input_json["tracernames"]) > 0:
            sveast = []
            for tracername in input_json["tracernames"]:
                sveast.append(
                    ds[f"{tracername}yz"]
                    .interp(yt=grid.yt)
                    .rename(f"{tracername}east")
                    .assign_coords(yt=grid.yt, zt=grid.zt)
                )
            # for isv in range(input_json['nsv']):
            #   with xr.open_mfdataset(f"{path}sv{isv+1:03}*",chunks={"time": input_json['tchunk']}) as ds:
            #     sveast.append(ds[f"sv{isv+1:03}"].interp(yt=grid.yt+input_json['y_offset']).rename('sveast').assign_coords(yt=grid.yt,zt=grid.zt))
            sveast = xr.merge(sveast)
            sveast = xr.concat([sveast0, sveast], dim="time")
    # south boundary
    path = Path(input_json["outpath_coarse"]) / f"crossxz.{iy_south+2:04d}"
    with xr.open_mfdataset(path.glob("*"), chunks={"time": input_json["tchunk"]}) as ds:
        usouth = (
            ds["uxz"]
            .interp(xm=grid.xm)
            .rename("usouth")
            .assign_coords(xm=grid.xm, zt=grid.zt)
        )
        usouth = xr.concat([usouth0, usouth], dim="time")
        vsouth = (
            ds["vxz"]
            .interp(xt=grid.xt)
            .rename("vsouth")
            .assign_coords(xt=grid.xt, zt=grid.zt)
        )
        vsouth = xr.concat([vsouth0, vsouth], dim="time")
        wsouth = (
            ds["wxz"]
            .interp(
                xt=grid.xt,
                zm=grid.zm,
                kwargs={"fill_value": "extrapolate"},
            )
            .rename("wsouth")
            .assign_coords(xt=grid.xt, zm=grid.zm)
        )
        wsouth = xr.concat([wsouth0, wsouth], dim="time")
        thlsouth = (
            ds["thlxz"]
            .interp(xt=grid.xt)
            .rename("thlsouth")
            .assign_coords(xt=grid.xt, zt=grid.zt)
        )
        thlsouth = xr.concat([thlsouth0, thlsouth], dim="time")
        qtsouth = (
            ds["qtxz"]
            .interp(xt=grid.xt)
            .rename("qtsouth")
            .assign_coords(xt=grid.xt, zt=grid.zt)
        )
        qtsouth = xr.concat([qtsouth0, qtsouth], dim="time")
        e12south = (
            ds["e120xz"]
            .interp(xt=grid.xt)
            .rename("e12south")
            .assign_coords(xt=grid.xt, zt=grid.zt)
        )
        e12south = xr.concat([e12south0, e12south], dim="time")
        if len(input_json["tracernames"]) > 0:
            svsouth = []
            for tracername in input_json["tracernames"]:
                svsouth.append(
                    ds[f"{tracername}xz"]
                    .interp(xt=grid.xt)
                    .rename(f"{tracername}south")
                    .assign_coords(xt=grid.xt, zt=grid.zt)
                )
            # for isv in range(input_json['nsv']):
            #   with xr.open_mfdataset(f"{path}sv{isv+1:03}*",chunks={"time": input_json['tchunk']}) as ds:
            #     svsouth.append(ds[f"sv{isv+1:03}"].interp(xt=grid.xt+input_json['x_offset']).rename('svsouth').assign_coords(xt=grid.xt,zt=grid.zt))
            svsouth = xr.merge(svsouth)
            svsouth = xr.concat([svsouth0, svsouth], dim="time")
    # north boundary
    path = Path(input_json["outpath_coarse"]) / f"crossxz.{iy_north+2:04d}"
    with xr.open_mfdataset(path.glob("*"), chunks={"time": input_json["tchunk"]}) as ds:
        unorth = (
            ds["uxz"]
            .interp(xm=grid.xm)
            .rename("unorth")
            .assign_coords(xm=grid.xm, zt=grid.zt)
        )
        unorth = xr.concat([unorth0, unorth], dim="time")
        vnorth = (
            ds["vxz"]
            .interp(xt=grid.xt)
            .rename("vnorth")
            .assign_coords(xt=grid.xt, zt=grid.zt)
        )
        vnorth = xr.concat([vnorth0, vnorth], dim="time")
        wnorth = (
            ds["wxz"]
            .interp(
                xt=grid.xt,
                zm=grid.zm,
                kwargs={"fill_value": "extrapolate"},
            )
            .rename("wnorth")
            .assign_coords(xt=grid.xt, zm=grid.zm)
        )
        wnorth = xr.concat([wnorth0, wnorth], dim="time")
        thlnorth = (
            ds["thlxz"]
            .interp(xt=grid.xt)
            .rename("thlnorth")
            .assign_coords(xt=grid.xt, zt=grid.zt)
        )
        thlnorth = xr.concat([thlnorth0, thlnorth], dim="time")
        qtnorth = (
            ds["qtxz"]
            .interp(xt=grid.xt)
            .rename("qtnorth")
            .assign_coords(xt=grid.xt, zt=grid.zt)
        )
        qtnorth = xr.concat([qtnorth0, qtnorth], dim="time")
        e12north = (
            ds["e120xz"]
            .interp(xt=grid.xt)
            .rename("e12north")
            .assign_coords(xt=grid.xt, zt=grid.zt)
        )
        e12north = xr.concat([e12north0, e12north], dim="time")
        if len(input_json["tracernames"]) > 0:
            svnorth = []
            for tracername in input_json["tracernames"]:
                svnorth.append(
                    ds[f"{tracername}xz"]
                    .interp(
                        xt=grid.xt,
                    )
                    .rename(f"{tracername}north")
                    .assign_coords(xt=grid.xt, zt=grid.zt)
                )
            # for isv in range(input_json['nsv']):
            #   with xr.open_mfdataset(f"{path}sv{isv+1:03}*",chunks={"time": input_json['tchunk']}) as ds:
            #     svnorth.append(ds[f"sv{isv+1:03}"].interp(xt=grid.xt+input_json['x_offset']).rename('svnorth').assign_coords(xt=grid.xt,zt=grid.zt))
            svnorth = xr.merge(svnorth)
            svnorth = xr.concat([svnorth0, svnorth], dim="time")
    # top boundary
    path = Path(input_json["outpath_coarse"]) / f"crossxy.{grid.kmax:04d}"
    with xr.open_mfdataset(path.glob("*"), chunks={"time": input_json["tchunk"]}) as ds:
        utop = (
            ds["uxy"]
            .interp(xm=grid.xm, yt=grid.yt)
            .rename("utop")
            .assign_coords(xm=grid.xm, yt=grid.yt)
        )
        utop = xr.concat([utop0, utop], dim="time")
        vtop = (
            ds["vxy"]
            .interp(xt=grid.xt, ym=grid.ym)
            .rename("vtop")
            .assign_coords(xt=grid.xt, ym=grid.ym)
        )
        vtop = xr.concat([vtop0, vtop], dim="time")
        wtop = (
            ds["wxy"]
            .interp(xt=grid.xt, yt=grid.yt)
            .rename("wtop")
            .assign_coords(xt=grid.xt, yt=grid.yt)
        )
        wtop = xr.concat([wtop0, wtop], dim="time")
        thltop = (
            ds["thlxy"]
            .interp(xt=grid.xt, yt=grid.yt)
            .rename("thltop")
            .assign_coords(xt=grid.xt, yt=grid.yt)
        )
        thltop = xr.concat([thltop0, thltop], dim="time")
        qttop = (
            ds["qtxy"]
            .interp(xt=grid.xt, yt=grid.yt)
            .rename("qttop")
            .assign_coords(xt=grid.xt, yt=grid.yt)
        )
        qttop = xr.concat([qttop0, qttop], dim="time")
        e12top = (
            ds["e120xy"]
            .interp(xt=grid.xt, yt=grid.yt)
            .rename("e12top")
            .assign_coords(xt=grid.xt, yt=grid.yt)
        )
        e12top = xr.concat([e12top0, e12top], dim="time")
        if len(input_json["tracernames"]) > 0:
            svtop = []
            for tracername in input_json["tracernames"]:
                svtop.append(
                    ds[f"{tracername}xy"]
                    .interp(
                        xt=grid.xt,
                        yt=grid.yt,
                    )
                    .rename(f"{tracername}top")
                    .assign_coords(xt=grid.xt, yt=grid.yt)
                )
            # for isv in range(input_json['nsv']):
            #   with xr.open_mfdataset(f"{path}sv{isv+1:03}*",chunks={"time": input_json['tchunk']}) as ds:
            #     svtop.append(ds[f"sv{isv+1:03}"].interp(xt=grid.xt+input_json['x_offset'],yt=grid.yt+input_json['y_offset']).rename('svtop').assign_coords(xt=grid.xt,yt=grid.yt))
            svtop = xr.merge(svtop)
            svtop = xr.concat([svtop0, svtop], dim="time")
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


def timestep0(input_json, grid, ix_west, ix_east, iy_south, iy_north):
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
        path = Path(input_json["outpath_coarse_old"]) / f"crossyz.{ix_west+2:04d}"
        with xr.open_mfdataset(
            path.glob("*"), chunks={"time": input_json["tchunk"]}
        ) as ds:
            uwest0 = (
                ds["uyz"]
                .isel(time=-1, drop=True)
                .interp(yt=grid.yt)
                .rename("uwest")
                .assign_coords(yt=grid.yt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            vwest0 = (
                ds["vyz"]
                .isel(time=-1, drop=True)
                .interp(ym=grid.ym)
                .rename("vwest")
                .assign_coords(ym=grid.ym, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            wwest0 = (
                ds["wyz"]
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
                ds["thlyz"]
                .isel(time=-1, drop=True)
                .interp(yt=grid.yt)
                .rename("thlwest")
                .assign_coords(yt=grid.yt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            qtwest0 = (
                ds["qtyz"]
                .isel(time=-1, drop=True)
                .interp(yt=grid.yt)
                .rename("qtwest")
                .assign_coords(yt=grid.yt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            e12west0 = (
                ds["e120yz"]
                .isel(time=-1, drop=True)
                .interp(yt=grid.yt)
                .rename("e12west")
                .assign_coords(yt=grid.yt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            if len(input_json["tracernames"]) > 0:
                svwest0 = []
                for tracername in input_json["tracernames"]:
                    svwest0.append(
                        ds[f"{tracername}yz"]
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
        path = Path(input_json["outpath_coarse_old"]) / f"crossyz.{ix_east+2:04d}"
        with xr.open_mfdataset(
            path.glob("*"), chunks={"time": input_json["tchunk"]}
        ) as ds:
            ueast0 = (
                ds["uyz"]
                .isel(time=-1, drop=True)
                .interp(yt=grid.yt)
                .rename("ueast")
                .assign_coords(yt=grid.yt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            veast0 = (
                ds["vyz"]
                .isel(time=-1, drop=True)
                .interp(ym=grid.ym)
                .rename("veast")
                .assign_coords(ym=grid.ym, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            weast0 = (
                ds["wyz"]
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
                ds["thlyz"]
                .isel(time=-1, drop=True)
                .interp(yt=grid.yt)
                .rename("thleast")
                .assign_coords(yt=grid.yt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            qteast0 = (
                ds["qtyz"]
                .isel(time=-1, drop=True)
                .interp(yt=grid.yt)
                .rename("qteast")
                .assign_coords(yt=grid.yt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            e12east0 = (
                ds["e120yz"]
                .isel(time=-1, drop=True)
                .interp(yt=grid.yt)
                .rename("e12east")
                .assign_coords(yt=grid.yt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            if len(input_json["tracernames"]) > 0:
                sveast0 = []
                for tracername in input_json["tracernames"]:
                    sveast0.append(
                        ds[f"{tracername}yz"]
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
        path = Path(input_json["outpath_coarse_old"]) / f"crossxz.{iy_south+2:04d}"
        with xr.open_mfdataset(
            path.glob("*"), chunks={"time": input_json["tchunk"]}
        ) as ds:
            usouth0 = (
                ds["uxz"]
                .isel(time=-1, drop=True)
                .interp(xm=grid.xm)
                .rename("usouth")
                .assign_coords(xm=grid.xm, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            vsouth0 = (
                ds["vxz"]
                .isel(time=-1, drop=True)
                .interp(xt=grid.xt)
                .rename("vsouth")
                .assign_coords(xt=grid.xt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            wsouth0 = (
                ds["wxz"]
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
                ds["thlxz"]
                .isel(time=-1, drop=True)
                .interp(xt=grid.xt)
                .rename("thlsouth")
                .assign_coords(xt=grid.xt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            qtsouth0 = (
                ds["qtxz"]
                .isel(time=-1, drop=True)
                .interp(xt=grid.xt)
                .rename("qtsouth")
                .assign_coords(xt=grid.xt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            e12south0 = (
                ds["e120xz"]
                .isel(time=-1, drop=True)
                .interp(xt=grid.xt)
                .rename("e12south")
                .assign_coords(xt=grid.xt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            if len(input_json["tracernames"]) > 0:
                svsouth0 = []
                for tracername in input_json["tracernames"]:
                    svsouth0.append(
                        ds[f"{tracername}xz"]
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
        path = Path(input_json["outpath_coarse_old"]) / f"crossxz.{iy_north+2:04d}"
        with xr.open_mfdataset(
            path.glob("*"), chunks={"time": input_json["tchunk"]}
        ) as ds:
            unorth0 = (
                ds["uxz"]
                .isel(time=-1, drop=True)
                .interp(xm=grid.xm)
                .rename("unorth")
                .assign_coords(xm=grid.xm, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            vnorth0 = (
                ds["vxz"]
                .isel(time=-1, drop=True)
                .interp(xt=grid.xt)
                .rename("vnorth")
                .assign_coords(xt=grid.xt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            wnorth0 = (
                ds["wxz"]
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
                ds["thlxz"]
                .isel(time=-1, drop=True)
                .interp(xt=grid.xt)
                .rename("thlnorth")
                .assign_coords(xt=grid.xt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            qtnorth0 = (
                ds["qtxz"]
                .isel(time=-1, drop=True)
                .interp(xt=grid.xt)
                .rename("qtnorth")
                .assign_coords(xt=grid.xt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            e12north0 = (
                ds["e120xz"]
                .isel(time=-1, drop=True)
                .interp(xt=grid.xt)
                .rename("e12north")
                .assign_coords(xt=grid.xt, zt=grid.zt)
                .expand_dims({"time": [pd.Timestamp(input_json["start"])]}, axis=0)
            )
            if len(input_json["tracernames"]) > 0:
                svnorth0 = []
                for tracername in input_json["tracernames"]:
                    svnorth0.append(
                        ds[f"{tracername}xz"]
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
        path = Path(input_json["outpath_coarse_old"]) / f"crossxy.{grid.kmax:04d}"
        with xr.open_mfdataset(
            path.glob("*"), chunks={"time": input_json["tchunk"]}
        ) as ds:
            utop0 = (
                ds["uxy"]
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
                ds["vxy"]
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
                ds["wxy"]
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
                ds["thlxy"]
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
                ds["qtxy"]
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
                ds["e120xy"]
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
                        ds[f"{tracername}xy"]
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
