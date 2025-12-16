import pandas as pd
from helper_scripts.grids import GridDalesOpenBC, nesting_idx
from typing import Union, List, Dict
from helper_scripts.logging_wrapper import logwrap
import logging

logger = logging.getLogger(__name__)
logger.debug("Entered module: %s", __name__)


@logwrap
def load_var(
    ds,
    var: str,
    boundary: str,
    grid: GridDalesOpenBC,
    indices: nesting_idx,
    isel: Union[bool, None] = False,
    expand_dims=False,
    expand_dims_time0=None,
    var_postfix="",
    move_x0y0=True,
):
    var_dims_dic = {
        "u": {"x": "xm", "y": "yt", "z": "zt"},
        "v": {"x": "xt", "y": "ym", "z": "zt"},
        "w": {"x": "xt", "y": "yt", "z": "zm"},
        "default": {"x": "xt", "y": "yt", "z": "zt"},
    }
    assign_grid_dic = {
        "xt": grid.xt,
        "xm": grid.xm,
        "yt": grid.yt,
        "ym": grid.ym,
        "zt": grid.zt,
        "zm": grid.zm,
    }
    if move_x0y0:
        interp_grid_dic = {
            "xt": grid.xt - grid.x0 + (indices.subgrid_x0 - indices.supergrid_x0),
            "xm": grid.xm - grid.x0 + (indices.subgrid_x0 - indices.supergrid_x0),
            "yt": grid.yt - grid.y0 + (indices.subgrid_y0 - indices.supergrid_y0),
            "ym": grid.ym - grid.y0 + (indices.subgrid_y0 - indices.supergrid_y0),
            "zt": grid.zt,
            "zm": grid.zm,
        }
    else:
        interp_grid_dic = assign_grid_dic

    boundary_var_interp_dic = {
        "west": {"w": ["y", "z"], "default": ["y"]},
        "east": {"w": ["y", "z"], "default": ["y"]},
        "south": {"w": ["x", "z"], "default": ["x"]},
        "north": {"w": ["x", "z"], "default": ["x"]},
        "top": {"default": ["x", "y"]},
    }
    boundary_var_assign_dic = {
        "west": {"default": ["y", "z"]},
        "east": {"default": ["y", "z"]},
        "south": {"default": ["x", "z"]},
        "north": {"default": ["x", "z"]},
        "top": {"default": ["x", "y"]},
    }

    boundary_var_isel_dim_dic = {
        "west": {"default": "x"},
        "east": {"default": "x"},
        "south": {"default": "y"},
        "north": {"default": "y"},
        "top": {"default": "z"},
    }

    boundary_var_isel_index_dic = {
        "east": {"default": indices.ix_east},
        "west": {"default": indices.ix_west},
        "south": {"default": indices.iy_south},
        "north": {"default": indices.iy_north},
        "top": {"default": grid.kmax - 1, "w": grid.kmax},
    }

    interp_kwargs_dic = {"w": {"fill_value": "extrapolate"}}

    # mapping from variable to (mapping from dimension name to dimension (x->xt))
    # we get var_dims : (x->xt,y->ym...)
    if var in var_dims_dic:
        var_dims = var_dims_dic[var]
    else:
        var_dims = var_dims_dic["default"]

    # which dimension to use to interpolate. for all variables this is y, x, or xy, but w has yz, xz and xy
    # at the end, var_interp_dims is a list of xt or xm or zt and so forth
    var_interp_dims = []
    if var in boundary_var_interp_dic[boundary]:
        for el in boundary_var_interp_dic[boundary][var]:
            var_interp_dims.append(var_dims[el])
    else:
        for el in boundary_var_interp_dic[boundary]["default"]:
            var_interp_dims.append(var_dims[el])

    # which dimension to use to assign coords to. for all variables this is yz, xz, or xy
    # at the end, var_assign_dims is a list of xt or xm or zt and so forth
    var_assign_dims = []
    if var in boundary_var_assign_dic[boundary]:
        for el in boundary_var_assign_dic[boundary][var]:
            var_assign_dims.append(var_dims[el])

    else:
        for el in boundary_var_assign_dic[boundary]["default"]:
            var_assign_dims.append(var_dims[el])

    # which dimension to use to select boundary indices from. for all variables this is yz, xz, or xy
    # at the end, var_assign_dims is a list of xt or xm or zt and so forth
    if var in boundary_var_isel_dim_dic[boundary]:
        var_isel_dim = var_dims[boundary_var_isel_dim_dic[boundary][var]]
    else:
        var_isel_dim = var_dims[boundary_var_isel_dim_dic[boundary]["default"]]

    if var in boundary_var_isel_index_dic[boundary]:
        var_isel_index = boundary_var_isel_index_dic[boundary][var]
    else:
        var_isel_index = boundary_var_isel_index_dic[boundary]["default"]

    # the precise coordinates to which to interpolate. for interpolating we subtract x0 and y0 as dales doesn't use those for output
    # for assigning coords, we do use x0 and y0 so we have a different dict for this.
    interp_coords = {dim: interp_grid_dic[dim] for dim in var_interp_dims}
    assign_coords = {dim: assign_grid_dic[dim] for dim in var_assign_dims}

    if var in interp_kwargs_dic:
        interp_kwargs = interp_kwargs_dic[var]
    else:
        interp_kwargs = None

    if isinstance(isel, bool):
        if isel:
            isel_dict = {var_isel_dim: var_isel_index}
        else:
            isel_dict = {}
    elif isinstance(isel, dict):
        isel_dict = isel
    else:
        raise TypeError(f"Unknown type for isel {isel}")

    if expand_dims:
        if expand_dims_time0 is None:
            raise ValueError("Missing time0!")
        if not isel:
            logger.Warning("Are you sure you want to do this?")
        interpolated_ds = (
            ds[f"{var}{var_postfix}"]
            .isel(isel_dict, drop=True)
            .interp(interp_coords, kwargs=interp_kwargs)
            .rename(f"{var}{boundary}")
            .assign_coords(assign_coords)
        ).expand_dims({"time": [pd.Timestamp(expand_dims_time0)]}, axis=0)
    else:
        interpolated_ds = (
            ds[f"{var}{var_postfix}"]
            .isel(isel_dict, drop=True)
            .interp(interp_coords, kwargs=interp_kwargs)
            .rename(f"{var}{boundary}")
            .assign_coords(assign_coords)
        )

    if interpolated_ds.isnull().sum() > 0:
        raise ValueError(f"Found NaN in dataset {boundary} {var}")
    return interpolated_ds


def get_u0(input_json, grid: GridDalesOpenBC, indices: nesting_idx, ds):
    uwest0 = (
        ds["u0"]
        .isel(xm=indices.ix_west, drop=True)
        .interp(yt=grid.yt)
        .rename("uwest")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(yt=grid.yt, zt=grid.zt)
    )
    uwest0 = load_var(
        ds,
        "u",
        "west",
        grid,
        indices,
        isel=True,
        expand_dims=True,
        expand_dims_time0=input_json["time0"],
        var_postfix="0",
    )
    ueast0 = (
        ds["u0"]
        .isel(xm=indices.ix_east, drop=True)
        .interp(yt=grid.yt)
        .rename("ueast")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(yt=grid.yt, zt=grid.zt)
    )
    usouth0 = (
        ds["u0"]
        .isel(yt=indices.iy_south, drop=True)
        .interp(xm=grid.xm)
        .rename("usouth")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(xm=grid.xm, zt=grid.zt)
    )
    unorth0 = (
        ds["u0"]
        .isel(yt=indices.iy_north, drop=True)
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


def get_v0(input_json, grid, indices: nesting_idx, ds):
    veast0 = (
        ds["v0"]
        .isel(xt=indices.ix_east, drop=True)
        .interp(ym=grid.ym)
        .rename("veast")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(ym=grid.ym, zt=grid.zt)
    )
    vwest0 = (
        ds["v0"]
        .isel(xt=indices.ix_west, drop=True)
        .interp(ym=grid.ym)
        .rename("vwest")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(ym=grid.ym, zt=grid.zt)
    )
    vnorth0 = (
        ds["v0"]
        .isel(ym=indices.iy_north, drop=True)
        .interp(xt=grid.xt)
        .rename("vnorth")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(xt=grid.xt, zt=grid.zt)
    )
    vsouth0 = (
        ds["v0"]
        .isel(ym=indices.iy_south, drop=True)
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


def get_w0(input_json, grid: GridDalesOpenBC, indices: nesting_idx, ds):
    weast0 = (
        ds["w0"]
        .isel(xt=indices.ix_east, drop=True)
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
        .isel(xt=indices.ix_west, drop=True)
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
        .isel(yt=indices.iy_north, drop=True)
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
        .isel(yt=indices.iy_south, drop=True)
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


def get_thl0(input_json, grid, indices: nesting_idx, ds):
    thleast0 = (
        ds["thl0"]
        .isel(xt=indices.ix_east, drop=True)
        .interp(yt=grid.yt)
        .rename("thleast")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(yt=grid.yt, zt=grid.zt)
    )
    thlwest0 = (
        ds["thl0"]
        .isel(xt=indices.ix_west, drop=True)
        .interp(yt=grid.yt)
        .rename("thlwest")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(yt=grid.yt, zt=grid.zt)
    )
    thlnorth0 = (
        ds["thl0"]
        .isel(yt=indices.iy_north, drop=True)
        .interp(xt=grid.xt)
        .rename("thlnorth")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(xt=grid.xt, zt=grid.zt)
    )
    thlsouth0 = (
        ds["thl0"]
        .isel(yt=indices.iy_south, drop=True)
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


def get_qt0(input_json, grid: GridDalesOpenBC, indices: nesting_idx, ds):
    qteast0 = (
        ds["qt0"]
        .isel(xt=indices.ix_east, drop=True)
        .interp(yt=grid.yt)
        .rename("qteast")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(yt=grid.yt, zt=grid.zt)
    )
    qtwest0 = (
        ds["qt0"]
        .isel(xt=indices.ix_west, drop=True)
        .interp(yt=grid.yt)
        .rename("qtwest")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(yt=grid.yt, zt=grid.zt)
    )
    qtnorth0 = (
        ds["qt0"]
        .isel(yt=indices.iy_north, drop=True)
        .interp(xt=grid.xt)
        .rename("qtnorth")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(xt=grid.xt, zt=grid.zt)
    )
    qtsouth0 = (
        ds["qt0"]
        .isel(yt=indices.iy_south, drop=True)
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


def get_e120(input_json, grid: GridDalesOpenBC, indices: nesting_idx, ds):
    e12east0 = (
        ds["e120"]
        .isel(xt=indices.ix_east, drop=True)
        .interp(yt=grid.yt)
        .rename("e12east")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(yt=grid.yt, zt=grid.zt)
    )
    e12west0 = (
        ds["e120"]
        .isel(xt=indices.ix_west, drop=True)
        .interp(yt=grid.yt)
        .rename("e12west")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(yt=grid.yt, zt=grid.zt)
    )
    e12north0 = (
        ds["e120"]
        .isel(yt=indices.iy_north, drop=True)
        .interp(xt=grid.xt)
        .rename("e12north")
        .expand_dims({"time": [pd.Timestamp(input_json["time0"])]}, axis=0)
        .assign_coords(xt=grid.xt, zt=grid.zt)
    )
    e12south0 = (
        ds["e120"]
        .isel(yt=indices.iy_south, drop=True)
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
