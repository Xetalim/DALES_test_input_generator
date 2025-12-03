# Crops the HARMONIE data to the required time and spatial extends.
# Transforms pressure coordinates into height levels.
# Transforms HARMONIE prognostic variables to DALES prognostic variables
import numpy as np
import xarray as xr
import dask
import warnings
from helper_scripts.LBC.Transform import Transform
import helper_scripts.LBC.hybrid_levels as hybrid_levels
from helper_scripts.LBC.helper import (
    calcBaseprof,
    differentiate,
    interp_z,
    load_data,
)

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


def prep_harmonie(input_json, grid):
    variables = ["ua", "va", "wa", "ta", "hus", "clw", "ps", "tas", "huss"]
    if "synturb" in input_json:
        variables = variables + ["tke", "tauu", "tauv", "cb", "hfss"]
    data, transform, x_sw, y_sw = create_xarray_dataset_POLYTOPE(
        input_json, grid, variables
    )
    # Change transform parameters to new DALES origin and update transform
    transform = update_transform(transform, x_sw, y_sw)

    # Calculate pressure levels
    p = calculate_pressure(input_json, data)
    print(p.mean(dim=["x", "y"]).values)
    print(data.lev.values)
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
            "z3d": xr.concat(z, dim="lev")
            .chunk({"lev": data.sizes["lev"] + 1})
            .transpose("time", "lev", "y", "x")
        }
    )
    # Get reference height levels (mean of height field first time step) and crop to grid.zsize
    z_int = get_ref_height_crop(input_json, grid, data)
    # Interpolate data to reference height levels
    data = interpolate_ref_height(input_json, data, z_int)

    # Calculate qt
    data = data.assign({"qt": data["clwc"] + data["q"]})
    # Calculate base profiles and exnr function
    ps_exnr, exnrs, thls_exnr, exnr = calc_base_exner(input_json, grid, data, z_int)

    print(f"ps = {ps_exnr}")
    print(f"thls = {thls_exnr}")
    input_json["ps"] = ps_exnr  # save the ps value used for the profile
    input_json["thls"] = thls_exnr  # and thls
    exnr = xr.DataArray(
        np.concatenate([exnrs[None], exnr]),
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
    return data, transform


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
    data = xr.merge(
        [
            xr.concat(
                [
                    data[dic[var]].assign_coords({"lev": data["lev"]})
                    # .expand_dims({"lev": [data.sizes["lev"] + 1]}, axis=1),
                    # data[dic[var] + "s"].expand_dims(
                    #     {"lev": [data.sizes["lev"] + 1]}, axis=1
                    # ),
                ],
                dim="lev",
            ).chunk({"lev": data.sizes["lev"] + 1})
            for var in variables
        ]
    )

    return data


def get_ref_height_crop(input_json, grid, data):
    if input_json["start"] == input_json["time0"]:  # Define reference height levels
        z_int = (
            data["z3d"]
            .isel({"time": 0}, drop=True)
            # .sel(x=slice(0, grid.xsize), y=slice(0, grid.ysize))
            .mean(dim=["x", "y"])  # [::-1]
            .values
        )
        try:
            arg = np.argwhere(z_int < grid.zsize).max()
        except ValueError:
            arg = -2
        z_int = z_int[: arg + 1]
    else:  # Take reference height levels from exnr.inp.xxx
        exnr = np.loadtxt(input_json["exnr_file"], skiprows=1)
        z_int = exnr[:, 0]
    return z_int


def calculate_turbulence_vars(
    grid,
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


def calc_base_exner(input_json, grid, data, z_int):
    if input_json["start"] == input_json["time0"]:  # Calculate exnr function
        z_min = data.z.argmin(dim="z")
        tas_exnr = (
            data["t"]
            .isel({"time": 0, "z": z_min}, drop=True)
            # .sel(x=slice(0, grid.xsize), y=slice(0, grid.ysize))
            .mean(dim=["x", "y"])
            .values
        )
        print(data.z.isel(z=z_min))
        print(data.z.isel(z=0))
        ps_exnr = (
            data["p"]
            .isel({"time": 0, "z": z_min}, drop=True)
            # .sel(x=slice(0, grid.xsize), y=slice(0, grid.ysize))
            .mean(dim=["x", "y"])
            .values
        )
        exnrs = (ps_exnr / p0) ** (Rd / cp)
        thls_exnr = tas_exnr / exnrs
        rhobf = calcBaseprof(z_int, thls_exnr, ps_exnr, pref0=p0)
        p_exnr = (
            rhobf[1:]
            * Rd
            * data["t"][0, 1:]
            # .sel(x=slice(0, grid.xsize), y=slice(0, grid.ysize))
            .mean(dim=["x", "y"]).values
            * (
                1
                + (Rv / Rd - 1) * data["qt"][0, 1:]
                # .sel(x=slice(0, grid.xsize), y=slice(0, grid.ysize))
                .mean(dim=["x", "y"]).values
                - Rv / Rd * data["clwc"][0, 1:]
                # .sel(x=slice(0, grid.xsize), y=slice(0, grid.ysize))
                .mean(dim=["x", "y"]).values
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
    for untranslated_var in variables:
        var = dic[untranslated_var]
        its = 0
        var_intz = []
        # Loop over time chunks (allows for parallel calculation in time index)
        for tchunk in data.chunks["time"]:
            ite = its + tchunk
            data_slice = dask.delayed(load_data)(
                data[var], {"time": np.arange(its, ite)}, drop=False
            )
            z_slice = dask.delayed(load_data)(
                data["z3d"], {"time": np.arange(its, ite)}, drop=False
            )
            var_intz.append(dask.delayed(interp_z)(z_slice, data_slice, z_int))
            its = ite
        # Concatenate data along time chunks and convert back to xarray
        var_intz = dask.delayed(np.concatenate)(var_intz, axis=0)
        var_intz = xr.DataArray(
            dask.array.from_delayed(var_intz, new_shape, dtype=float),
            dims=new_dims,
            coords=new_coords,
            name=var,
            attrs=data[var].attrs,
        ).chunk({"time": input_json["tchunk"]})
        data_intz.append(var_intz)
    # Store interpolated height data in a DataSet
    data = xr.merge(data_intz).assign({"x": data.x, "y": data.y})
    return data


def calculate_3d_height_levels(data):
    rho = data["p"] / (
        Rd
        * data["t"]
        * (1 + (Rv / Rd - 1) * (data["q"] + data["clwc"]) - Rv / Rd * data["clwc"])
    )
    rhoh = 0.5 * (rho.assign_coords({"lev": rho["lev"].values - 1}) + rho)
    z = [xr.zeros_like(data["p"].isel(lev=1, drop=True).rename("z3d"))]
    for k in np.arange(data.sizes["lev"] - 2, -1, -1):
        z = [
            z[0]
            - (data["p"].isel(lev=k, drop=True) - data["p"].isel(lev=k + 1, drop=True))
            / (rhoh.isel(lev=k, drop=True) * grav)
        ] + z

    return z


def calculate_pressure(input_json, data):

    # hybrid_coeff = np.loadtxt(f"{input_json['inpath']}H43_65lev.txt")
    hybrid_A = xr.DataArray(
        hybrid_levels.ahalf, dims=["lev"], coords={"lev": np.arange(1, 92)}
    )
    hybrid_B = xr.DataArray(
        hybrid_levels.bhalfs, dims=["lev"], coords={"lev": np.arange(1, 92)}
    )
    ph = (hybrid_A + hybrid_B * data["msl"]).transpose("time", "lev", "y", "x")
    # calculate on pressure levels
    p = 0.5 * (ph.assign_coords({"lev": ph["lev"].values - 1}) + ph)
    return p.sel(lev=data["lev"].values)


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


def create_xarray_dataset(input_json, grid, variables):
    data = []
    # Get time epochs

    # Open data and crop data
    x_sw, y_sw = grid.x0, grid.y0
    # first get the time and transform data....
    var = variables[0]
    with xr.open_mfdataset(
        f"{input_json['inpath']}{var}_*.nc", chunks={"time": input_json["tchunk"]}
    ) as ds:
        transform, _, _, time = get_transform_time(input_json, var, ds)

    for var in variables:
        with xr.open_mfdataset(
            f"{input_json['inpath']}{var}_*.nc", chunks={"time": input_json["tchunk"]}
        ) as ds:
            # Interpolate fluxes to same time
            if var in ["tauu", "tauv", "hfss"]:
                ds = ds.interp(
                    time=time, assume_sorted=True, kwargs={"fill_value": "extrapolate"}
                ).chunk({"time": input_json["tchunk"]})
            # Crop data to time and spatial range, using harmonie spatial resolution or filter as buffer
            dx = ds["x"][1].values - ds["x"][0].values
            dy = ds["y"][1].values - ds["y"][0].values
            if "filter" in input_json:  # add some extra width for gaussian filtering
                buffer = 4 * input_json["filter"]["sigma"]
            else:
                buffer = dx
            data.append(
                ds[var].sel(
                    time=slice(input_json["start"], input_json["end"]),
                    x=slice(x_sw - buffer, x_sw + grid.xsize + buffer),
                    y=slice(y_sw - buffer, y_sw + grid.ysize + buffer),
                )
            )
        # Set south west corner of DALES as origin
        data[-1] = data[-1].assign_coords(
            {"x": data[-1]["x"].values - x_sw, "y": data[-1]["y"].values - y_sw}
        )
    # Merge into xarray dataset
    data = xr.merge(data, compat="override")
    return data, transform, x_sw, y_sw


def create_xarray_dataset_POLYTOPE(input_json, grid, variables):
    data = []
    # Get time epochs

    # Open data and crop data

    # first get the time and transform data....
    # ["ua", "va", "wa", "ta", "hus", "clw", "ps", "tas", "huss"]
    x_sw, y_sw = grid.x0, grid.y0
    var = variables[0]
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*formula_terms")
        with xr.open_mfdataset(
            "/Users/andrevanginkel/Documents/20_Code/21 Input_Output scripts/21.03_paris_NWP/atmo_new.nc",
            decode_coords="all",
        ) as ds:
            transform, _, _, time = get_transform_time(input_json, var, ds)

            print(f"Got {(x_sw,y_sw)=}")
            print(f"{(np.min(grid.xm), np.max(grid.xm))=}")
            print(f"{(np.min(grid.ym), np.max(grid.ym))=}")
            print(f"{(np.min(ds.x).values, np.max(ds.x).values)=}")
            print(f"{(np.min(ds.y).values, np.max(ds.y).values)=}")
            # exit()
    for var_raw in variables:
        if var_raw in ["huss", "ps", "tas"]:
            filename = "/Users/andrevanginkel/Documents/20_Code/21 Input_Output scripts/21.03_paris_NWP/surface_new.nc"
        else:
            filename = "/Users/andrevanginkel/Documents/20_Code/21 Input_Output scripts/21.03_paris_NWP/atmo_new.nc"
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
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*formula_terms")
            with xr.open_mfdataset(filename, decode_coords="all") as ds:
                # Interpolate fluxes to same time
                if var in ["tauu", "tauv", "hfss"]:
                    ds = ds.interp(
                        time=time,
                        assume_sorted=True,
                        kwargs={"fill_value": "extrapolate"},
                    ).chunk({"time": input_json["tchunk"]})

                # make sure we all have the same dates and times in the files
                ds = ds.sel(time=time)
                # Crop data to time and spatial range, using harmonie spatial resolution or filter as buffer
                dx = ds["x"][1].values - ds["x"][0].values
                dy = ds["y"][1].values - ds["y"][0].values
                if (
                    "filter" in input_json
                ):  # add some extra width for gaussian filtering
                    buffer = 4 * input_json["filter"]["sigma"]
                else:
                    buffer = dx
                data.append(
                    ds[var].sel(  # also step TODO
                        x=slice(int(x_sw - buffer), int(x_sw + grid.xsize + buffer)),
                        y=slice(
                            int(y_sw - buffer), int(y_sw + grid.ysize + buffer)
                        ),  # TODO INT
                    )
                    # .isel(x=slice(0, 10), y=slice(0, 10))
                )
        # DO NOT SET SOUTH WEST CORNER OF DALES AS ORIGIN
        # Set south west corner of DALES as origin
        # data[-1] = data[-1].assign_coords(
        #     {"x": data[-1]["x"].values - x_sw, "y": data[-1]["y"].values - y_sw}
        # )
    # Merge into xarray dataset
    data = xr.merge(data, compat="override")
    return data, transform, x_sw, y_sw


def get_transform_time(input_json, var, ds):
    # Read transform information and transform lat/lon of southwest corner to harmonie x/y
    # proj = "+proj=eqc +ellps=WGS84 +a=6378137.0 +lon_0=0.0 +to_meter=111319.4907932736 +no_defs +type=crs"
    proj = ds.rio.crs.to_proj4()
    transform = Transform({"proj4": proj})
    x_sw, y_sw = transform.latlon_to_xy(input_json["lat_sw"], input_json["lon_sw"])
    # Round to 5 meters to avoid numerical error in coordinates
    x_sw = np.round(x_sw, 0)
    y_sw = np.round(y_sw, 0)
    time = ds["time"].values
    return transform, x_sw, y_sw, time
