import numpy as np
import xarray as xr
import warnings

from datetime import datetime
from helper_scripts.grids import GridDalesOpenBC, nesting_idx
from helper_scripts.LBC.dales_nest.get_all_dales_boundaries import (
    get_all_dales_boundaries,
)


def boundary_fields_fine(
    input_json, grid: GridDalesOpenBC, output_path, grid_indices: nesting_idx
):
    if not grid_indices:
        raise ValueError("No nesting indices provided!")

    ix_west = grid_indices.ix_west
    ix_east = grid_indices.ix_east
    iy_south = grid_indices.iy_south
    iy_north = grid_indices.iy_north

    # Get initial boundary fields from initial fields
    (
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
    ) = get_all_dales_boundaries(input_json, grid, ix_west, ix_east, iy_south, iy_north)
    # Add fields to dataset
    if len(input_json["tracernames"]) > 0:
        openboundaries = xr.merge(
            [
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
            ],
            combine_attrs="drop",
        )
    else:
        openboundaries = xr.merge(
            [
                uwest,
                vwest,
                wwest,
                thlwest,
                qtwest,
                e12west,
                ueast,
                veast,
                weast,
                thleast,
                qteast,
                e12east,
                usouth,
                vsouth,
                wsouth,
                thlsouth,
                qtsouth,
                e12south,
                unorth,
                vnorth,
                wnorth,
                thlnorth,
                qtnorth,
                e12north,
                utop,
                vtop,
                wtop,
                thltop,
                qttop,
                e12top,
            ],
            combine_attrs="drop",
        )
    openboundaries = set_openboundary_attrs(input_json, openboundaries)
    # Add global attributes
    openboundaries = openboundaries.assign_attrs(
        {
            "title": f"openboundaries.inp.{input_json['iexpnr']:03d}.nc",
            "history": f"Created on {datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')} UTC",
            "author": input_json["author"],
            "time0": input_json["time0"],
        }
    )
    openboundaries.to_netcdf(
        path=output_path / "input" / openboundaries.attrs["title"],
        mode="w",
        format="NETCDF4",
    )
    return openboundaries


def set_openboundary_attrs(input_json, openboundaries):
    dts = (
        openboundaries.time.values.astype("datetime64[s]")
        - np.datetime64(input_json["time0"], "s")
    ) / np.timedelta64(1, "s")
    openboundaries = openboundaries.assign_coords({"time": ("time", dts)})
    # # Adjust time variable to seconds since initial field
    # ts = openboundaries['time'].values.astype('datetime64[s]')
    # dts = (ts-np.datetime64(input_json['time0'],'s'))/np.timedelta64(1, 's')
    # openboundaries = openboundaries.assign_coords({'time':('time', dts)})
    # Add variable attributes
    openboundaries["time"] = openboundaries["time"].assign_attrs(
        {"longname": "Time"}
    )  # , 'units': f"seconds since {input_json['time0']}"})
    openboundaries.time.encoding["units"] = f"seconds since {input_json['time0']}"
    openboundaries["xt"] = openboundaries["xt"].assign_attrs(
        {"longname": "West-East displacement of cell centers", "units": "m"}
    )
    openboundaries["xm"] = openboundaries["xm"].assign_attrs(
        {"longname": "West-East displacement of cell edges", "units": "m"}
    )
    openboundaries["yt"] = openboundaries["yt"].assign_attrs(
        {"longname": "South-North displacement of cell centers", "units": "m"}
    )
    openboundaries["ym"] = openboundaries["ym"].assign_attrs(
        {"longname": "South-North displacement of cell edges", "units": "m"}
    )
    openboundaries["zt"] = openboundaries["zt"].assign_attrs(
        {"longname": "Vertical displacement of cell centers", "units": "m"}
    )
    openboundaries["zm"] = openboundaries["zm"].assign_attrs(
        {"longname": "Vertical displacement of cell edges", "units": "m"}
    )
    variables = ["u", "v", "w", "thl", "qt", "e12"]
    if len(input_json["tracernames"]) > 0:
        for tracername in input_json["tracernames"]:
            variables.append(tracername)
    units = ["m/s", "m/s", "m/s", "K", "kg/kg", "m/s"]
    if len(input_json["tracernames"]) > 0:
        for tracername in input_json["tracernames"]:
            if len(input_json["tracernames"]) > 2:
                warnings.warn(
                    f"Unit applied to tracer {tracername} might not be correct!"
                )
            units.append("kg/kg")
    long_names = [
        "West-East velocity at ",
        "South-North velocity at ",
        "Vertical velocity at ",
        "Liquid water potential temperature at ",
        "Total water specific humidity at ",
        "Square root of turbulent kinetic energy at ",
    ]
    if len(input_json["tracernames"]) > 0:
        for tracername in input_json["tracernames"]:
            long_names.append(f"scalar field {tracername} at ")
    for ivar, var in enumerate(variables):
        unit = units[ivar]
        long_name = long_names[ivar]
        for boundary in ["West", "East", "South", "North", "top"]:
            openboundaries[var + boundary.lower()] = openboundaries[
                var + boundary.lower()
            ].assign_attrs(
                {"longname": long_name + boundary + " boundary", "units": unit}
            )

    return openboundaries
