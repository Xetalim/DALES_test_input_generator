# Crops the HARMONIE data to the required time and spatial extends.
# Transforms pressure coordinates into height levels.
# Transforms HARMONIE prognostic variables to DALES prognostic variables
import numpy as np
import xarray as xr
import dask
import warnings
from helper_scripts.LBC.Transform import Transform
import helper_scripts.LBC.hybrid_levels as hybrid_levels
from helper_scripts.grids import GridDalesOpenBC
from helper_scripts.logging_wrapper import logwrap
from helper_scripts.LBC.helper import (
    calcBaseprof,
    differentiate,
    interp_z,
    load_data,
)

import logging

logger = logging.getLogger(__name__)
logger.debug("Entered module: %s", __name__)


xr.set_options(use_new_combine_kwarg_defaults=True)
DATA_LAT_VAR = "latitude"
DATA_LON_VAR = "longitude"
# Constants
# Simulation constants (modglobal)
p0 = 1e5  # Reference pressure
Rd = 287.04  # Gas constant for dry air
Rv = 461.5  # Gas constant for water vapor
cp = 1004.0  # Specific heat at constant pressure (dry air)
Lv = 2.53e6  # Latent heat for vaporisation
grav = 9.81  # Gravitational constant
kappa = 0.4  # Von Karman constant


@logwrap
def prep_harmonie(input_json, grid: GridDalesOpenBC):
    variables = ["ua", "va", "wa", "ta", "hus", "clw", "ps", "tas", "huss"]
    if "synturb" in input_json:
        variables = variables + ["tke", "tauu", "tauv", "cb", "hfss"]
    data, transform, x_sw, y_sw = create_xarray_dataset(
        input_json, grid, variables
    )

    data, = dask.optimize(data)#.chunk({"x":"auto","y":"auto","time":5,"lev":"auto"}))
    # Change transform parameters to new DALES origin and update transform
    transform = update_transform(transform, x_sw, y_sw)

    # Calculate pressure levels
    p = calculate_pressure(input_json, data)
    # print(p.mean(dim=["x", "y"]).values)
    # print(data.lev.values)
    # drop NAN if levels from harmonie is less than 90
    data = data.assign({"p": p})  # .dropna(dim="lev")
    # print(data.assign({"p": p}).p.values)
    # print(data.assign({"p": p}).p.values)
    # Add missing surface fields to 3d fields
    variables = ["uas", "vas", "was", "clws"]
    if "synturb" in input_json:
        variables.append("tkes")
    data = data.assign({var: xr.zeros_like(data["msl"]) for var in variables})

    if "synturb" in input_json:
        turbulence_data_dic = {
            "tauu": data["tauu"],
            "tauv": data["tauv"],
            "hfss": data["hfss"],
            "cb": data["cb"],
        }

    # Concatenate surface and 3D fields
    variables = ["ua", "va", "wa", "ta", "hus", "clw", "p"]
    if "synturb" in input_json:
        variables.append("tke")
    data = merge_steps(data, variables)

    # Calculate 3D height levels

    z = calculate_3d_height_levels(data)
    data = data.assign(
        {
            "z3d": z
        }
    )
    # Get reference height levels (mean of height field first time step) and crop to grid.zsize
    z_int = get_ref_height_crop(input_json, grid, data)

    data, = dask.optimize(data)

    z_int = z_int.compute()

    print(z_int)

    
    # Interpolate data to reference height levels
    data = interpolate_ref_height(input_json, data, z_int)

    # make sure z_int is also dimension z now..
    z_int = z_int.rename({"lev": "z"})

    # Calculate qt
    data = data.assign({"qt": data["clwc"] + data["q"]})
    
    # Calculate base profiles and exnr function
    ps_exnr, exnrs, thls_exnr, exnr = calc_base_exner(input_json, grid, data, z_int)


    input_json["ps"] = float(ps_exnr)  # save the ps value used for the profile
    input_json["thls"] = float(thls_exnr)  # and thls
    exnr = xr.DataArray(
        np.concatenate([[exnrs], exnr]),
        dims=["z"],
        coords={"z": z_int},
        name="exnr",
        attrs={"thls": thls_exnr, "ps": ps_exnr},
    )

    data = data.assign({"exnr": exnr})

    # Calculate liquid potential temperature and total specific humidity
    thl = data["t"] / exnr - Lv * data["clwc"] / (cp * exnr)
    data = data.assign({"thl": thl})
    # Calculate turbulence parameters
    if "synturb" in input_json:
        calculate_turbulence_vars(
            grid,
            data,
            turbulence_data_dic,
            z_int,
            ps_exnr,
            exnrs,
            thls_exnr,
            thl,
        )
    # Organize data, rename variables and drop non DALES prognostic variables
    data = (
        data.rename({"wz": "w"})
        .drop(["t", "p", "clwc", "q"])
        .assign(
            {
                "transform": xr.DataArray(
                    [], name="Lambert_Conformal", attrs=transform.parameters
                )
            }
        )
    )
    data, = dask.optimize(data)
    return data, transform


@logwrap
def merge_steps(data, variables):
    dic = {
        "ua": "u",
        "va": "v",
        "wa": "wz",
        "ta": "t",
        "hus": "q",
        "clw": "clwc",
        "ps": "msl",
        "tas": "2t",
        "huss": "2sh",
        "p": "p",
    }
    data = data[[dic[var] for var in variables]]

    return data


@logwrap
def get_ref_height_crop(input_json, grid: GridDalesOpenBC, data):
    if input_json["start"] == input_json["time0"]:  # Define reference height levels
        z_int = (
            data["z3d"].isel({"time": 0}, drop=True)
            # .sel(x=slice(0, grid.xsize), y=slice(0, grid.ysize))
            .mean(dim=["x", "y"])  # [::-1]
        )
    else:  # Take reference height levels from exnr.inp.xxx
        exnr = np.loadtxt(input_json["exnr_file"], skiprows=1)
        z_int = exnr[:, 0]
    return z_int


@logwrap
def calculate_turbulence_vars(
    grid: GridDalesOpenBC,
    data,
    turbulence_data_dic,
    z_int,
    ps_exnr,
    exnrs,
    thls_exnr,
    thl,
):

    # Calculate inversion height from maximum curvature and with cloud base as a backup
    zi_min = 200
    zi_max = 4000
    thlmean = (
        thl
        #    .sel(x=slice(0, grid.xsize), y=slice(0, grid.ysize))
        .mean(dim=["x", "y"])
    )
    cbmean = (
        xr.where(turbulence_data_dic["cb"] > 0.0, turbulence_data_dic["cb"], np.NaN)
        # .sel(x=slice(0, grid.xsize), y=slice(0, grid.ysize))
        .mean(dim=["x", "y"])
    )
    its = 0
    d2thl = []
    d1thl = []
    rhobf = calcBaseprof(z_int, thls_exnr, ps_exnr, pref0=p0)
    for tchunk in thlmean.chunks[0]:
        ite = its + tchunk
        d1thl.append(
            dask.delayed(differentiate)(
                thlmean.isel({"time": np.arange(its, ite)}), "z", 1, acc=6
            )
        )
        d2thl.append(
            dask.delayed(differentiate)(
                thlmean.isel({"time": np.arange(its, ite)}), "z", 2, acc=6
            )
        )
        its = ite
    d1thl = dask.delayed(xr.concat)(d2thl, dim="time").compute()
    d2thl = dask.delayed(xr.concat)(d1thl, dim="time").compute()
    zi = d2thl.where(d1thl > 0).sel(z=slice(zi_min, zi_max)).idxmax("z").fillna(cbmean)
    rhobs = rhobf[0] - z_int[0] * (rhobf[1] - rhobf[0]) / (z_int[1] - z_int[0])
    ustar = np.sqrt(np.maximum(turbulence_data_dic["tauu"], 0) / rhobs).rename("ustar")
    vstar = np.sqrt(np.maximum(turbulence_data_dic["tauv"], 0) / rhobs).rename("vstar")
    wthls = (turbulence_data_dic["hfss"] / (exnrs * rhobs * cp)).rename("wthls")
    data = data.assign({"ustar": ustar, "vstar": vstar, "wthls": wthls, "zi": zi})


@logwrap
def calc_base_exner(input_json, grid: GridDalesOpenBC, data, z_int):
    if input_json["start"] == input_json["time0"]:  # Calculate exnr function
        z_min = data.z.argmin(dim="z")
        tas_exnr = (
            data["t"].isel({"time": 0, "z": z_min}, drop=True)
            # .sel(x=slice(0, grid.xsize), y=slice(0, grid.ysize))
            .mean(dim=["x", "y"])
        )

        ps_exnr = (
            data["p"].isel({"time": 0, "z": z_min}, drop=True)
            # .sel(x=slice(0, grid.xsize), y=slice(0, grid.ysize))
            .mean(dim=["x", "y"])
        )
        exnrs = (ps_exnr / p0) ** (Rd / cp)
        thls_exnr = tas_exnr / exnrs
        # somehow this function doesn't like being dasked. max size of array going in is (lev,) or (1,) so shouldn't be too big of a problem?
        # if ever you have performance issues, look at among others this function..
        rhobf = calcBaseprof(z_int.compute(), thls_exnr.compute(), ps_exnr.compute(), pref0=p0)
        p_exnr = (
            rhobf[1:]
            * Rd
            * data["t"].isel(time=0, z=slice(1, None, None))
            # .sel(x=slice(0, grid.xsize), y=slice(0, grid.ysize))
            .mean(dim=["x", "y"])
            * (
                1
                + (Rv / Rd - 1) * data["qt"][0, 1:]
                # .sel(x=slice(0, grid.xsize), y=slice(0, grid.ysize))
                .mean(dim=["x", "y"]).values
                - Rv / Rd * data["clwc"][0, 1:]
                # .sel(x=slice(0, grid.xsize), y=slice(0, grid.ysize))
                .mean(dim=["x", "y"])
            )
        )  # Ideal gas law
        exnr = (p_exnr / p0) ** (Rd / cp)
    else:  # Read exnr.inp.xxx
        with open(input_json["exnr_file"], "r") as file:
            line0 = file.readline()
        thls_exnr = float(line0.split(",")[1].split("thls = ")[-1])
        ps_exnr = float(line0.split(",")[2].split("ps = ")[-1])
        exnr = np.loadtxt(input_json["exnr_file"], skiprows=1)
        exnrs = exnr[0, 1]
        exnr = exnr[1:, 1]
    return ps_exnr, exnrs, thls_exnr, exnr

@logwrap
def vertical_interp_all(ds, z3d, z_new):
    """
    Interpolate all variables in ds along a 4D vertical coordinate z3d
    onto target heights z_new.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with dims (time, x, y, lev)
    z3d : xarray.DataArray
        4D array with actual heights (time, x, y, lev)
    z_new : 1D array
        Target vertical levels (physical heights)

    Returns
    -------
    ds_interp : xarray.Dataset
        Dataset with all variables interpolated to z_new
    """
    z_target = z_new

    def _interp_func(var_col, z_col, z_target=z_target):
        # var_col, z_col: 1D arrays of shape (lev,)
        return np.interp(z_target, z_col, var_col, left=np.nan, right=np.nan)

    interp_vars = {}
    for var in ds.data_vars:
        da = ds[var]
        # Apply vertical interpolation along lev dimension
        da_interp = xr.apply_ufunc(
            _interp_func,
            da,
            z3d,
            input_core_dims=[["lev"], ["lev"]],
            output_core_dims=[["lev"]],
            vectorize=True,
            dask="parallelized",
            output_dtypes=[da.dtype],
        )
        da_interp = da_interp.assign_coords(lev=z_target)
        interp_vars[var] = da_interp

    # Build new dataset
    ds_interp = xr.Dataset(
        interp_vars, coords={k: v for k, v in ds.coords.items() if k != "lev"}
    )
    ds_interp = ds_interp.assign_coords(lev=z_target)

    return ds_interp.rename({"lev": "z"})

@logwrap
def interpolate_ref_height(input_json, data, z_int):
    data_intz = []
    its = 0
    new_shape = [data.sizes["time"], len(z_int), data.sizes["y"], data.sizes["x"]]
    new_coords = {
        "time": data.coords["time"],
        "z": z_int,
        "y": data.coords["y"],
        "x": data.coords["x"],
    }
    new_dims = ["time", "z", "y", "x"]
    variables = ["ua", "va", "wa", "ta", "p", "hus", "clw"]
    dic = {
        "ua": "u",
        "va": "v",
        "wa": "wz",
        "ta": "t",
        "hus": "q",
        "clw": "clwc",
        "ps": "msl",
        "tas": "2t",
        "huss": "2sh",
        "p": "p",
    }
    if "synturb" in input_json:
        variables.append("tke")
    logger.debug("Checking if data is ascending..")
    z_col = data["z3d"].isel(time=0, x=0, y=0, lev=slice(0, 2))
    z0 = z_col.isel(lev=0)
    z1 = z_col.isel(lev=1)
    logger.debug("Checking if data is ascending.. SUCCEEDED")
    # descending = float(z1) < float(z0)
    if float(z1) < float(z0):
        data = data.isel(lev=slice(None, None, -1))
        logger.debug("Successfully inverted data!")
    else:
        pass

    data = vertical_interp_all(data, data["z3d"], z_int)
    return data


@logwrap
def calculate_3d_height_levels(data):
    rho = (
        data["p"]
        / (
            Rd
            * data["t"]
            * (1 + (Rv / Rd - 1) * (data["q"] + data["clwc"]) - Rv / Rd * data["clwc"])
        )
    )
    rhoh = 0.5 * (rho.assign_coords({"lev": rho["lev"].values - 1}) + rho)

    pdiff = data["p"].diff(dim="lev", label="lower")
    z = (
        (-pdiff / (rhoh * grav))
        .cumsum(dim="lev")
    ).reindex(lev=data.lev, fill_value=0) # make sure 0 is included
    z = z - z.isel(lev=0)

    # print(z.compute())
    return z


@logwrap
def calculate_pressure(input_json, data):
    # right now this function uses 90 hybrid model levels.
    # hybrid_coeff = np.loadtxt(f"{input_json['inpath']}H43_65lev.txt")
    hybrid_A = xr.DataArray(
        hybrid_levels.ahalf, dims=["lev"], coords={"lev": np.arange(1, 92)}
    )
    hybrid_B = xr.DataArray(
        hybrid_levels.bhalfs, dims=["lev"], coords={"lev": np.arange(1, 92)}
    )
    ph = (hybrid_A + hybrid_B * data["msl"]).transpose("time", "lev", "y", "x")
    # calculate on pressure levels
    p = 0.5 * (ph.assign_coords({"lev": ph["lev"] - 1}) + ph)
    return p.sel(lev=data["lev"])


@logwrap
def update_transform(transform, x_sw, y_sw):
    return transform
    transform.parameters["false_easting"] = transform.parameters["false_easting"] - x_sw
    transform.parameters["false_northing"] = (
        transform.parameters["false_northing"] - y_sw
    )
    proj4 = ""
    for param in transform.parameters["proj4"][1:].split("+"):
        line = "+" + param
        if "x_0" in param:
            line = f"+x_0={transform.parameters['false_easting']} "
        if "y_0" in param:
            line = f"+y_0={transform.parameters['false_northing']} "
        proj4 = proj4 + line
    transform.parameters["proj4"] = proj4.rstrip()
    transform = Transform(transform.parameters)
    return transform

@logwrap
def create_xarray_dataset(input_json, grid: GridDalesOpenBC, variables):
    data = []
    # Get time epochs

    # Open data and crop data

    # first get the time and transform data....
    # ["ua", "va", "wa", "ta", "hus", "clw", "ps", "tas", "huss"]
    x_sw, y_sw = grid.x0, grid.y0
    var = variables[0]
    ds_ml = xr.open_mfdataset(
        input_json["HARMONIE_ml_glob"],
        decode_coords="all",
        parallel=True,
        chunks={"time":input_json["tchunk"],"lev":-1},
        # chunks={"x": "auto", "y": "auto", "time": "auto", "lev": -1},
    ).drop_duplicates(dim="time")
    transform, _, _, time = get_transform_time(input_json, var, ds_ml)

    # print(ds_ml.lev)
    # print(time.values)

    ds_sfc = xr.open_mfdataset(
        input_json["HARMONIE_sfc_glob"],
        decode_coords="all",
        parallel=True,
        chunks={"time":input_json["tchunk"]},
        # chunks={"x": "auto", "y": "auto", "time": "auto", "lev": -1},
    )

    logger.debug("Succesfully read in HARMONIE_ml_glob")

    ds_sfc = ds_sfc.drop_duplicates(dim="time")
    
    # interpolate surface fluxes to higher time resolution, without assuming correct sorting as I've had problems with this before.
    ds_sfc = ds_sfc.interp(
        time=time,
        assume_sorted=False,
        kwargs={"fill_value": "extrapolate"},
    )
    logger.debug("Succesfully read in and interpolated HARMONIE_sfc_glob")
    
    for var_raw in variables:
        var = {
            "ua": "u",
            "va": "v",
            "wa": "wz",
            "ta": "t",
            "hus": "q",
            "clw": "clwc",
            "ps": "msl",
            "tas": "2t",
            "huss": "2sh",
        }[var_raw]
        logger.debug(f"Reading in variable {var}")

        # Crop data to time and spatial range, using harmonie spatial resolution or filter as buffer
        dx = ds_ml["x"][1] - ds_ml["x"][0]
        dy = ds_ml["y"][1] - ds_ml["y"][0]
        
        if "filter" in input_json:  # add some extra width for gaussian filtering
            buffer = 4 * input_json["filter"]["sigma"]
        else:
            buffer = dx

        # Interpolate fluxes and surface levels to same time
        if var in ["tauu", "tauv", "hfss", "msl", "2t", "2sh"]:
            data.append(ds_sfc[var])
        else:
            data.append(
                ds_ml[var]
            )
    # Merge into xarray dataset
    logger.debug("Succesfully read in vars, merging now...")
    data = xr.merge(data, compat="override", join="outer").sel(  # also step TODO
        time=time.sortby("time").sel(time=slice(input_json["start"],input_json["end"])),
        x=slice(int(x_sw - buffer), int(x_sw + grid.xsize + buffer)),
        y=slice(int(y_sw - buffer), int(y_sw + grid.ysize + buffer)),  # TODO INT
    )

    return data, transform, x_sw, y_sw


@logwrap
def get_transform_time(input_json, var, ds):
    # sometimes dask doesn't recognise that rioxarray has been imported. We import rioxarray to make sure we can acess rio.crs..
    import rioxarray

    # Read transform information and transform lat/lon of southwest corner to harmonie x/y
    proj = ds.rio.crs.to_proj4()
    transform = Transform({"proj4": proj})
    x_sw, y_sw = transform.latlon_to_xy(input_json["lat_sw"], input_json["lon_sw"])
    # Round to 5 meters to avoid numerical error in coordinates
    x_sw = np.round(x_sw, 0)
    y_sw = np.round(y_sw, 0)
    time = ds["time"]
    return transform, x_sw, y_sw, time
