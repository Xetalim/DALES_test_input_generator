import xarray as xr

from helper_scripts.grids import GridDalesOpenBC
from datetime import datetime
import logging
from helper_scripts.logging_wrapper import logwrap

logger = logging.getLogger(__name__)
logger.debug("Entered module: %s", __name__)


@logwrap
def initial_fields_fine(input_json, grid: GridDalesOpenBC, output_path):
    """
    Docstring for initial_fields_fine

    :param input_json: Configuration dict for LBC
    :param grid: Grid object describing the output grid for open boundaries, which is different from the normal DALES grid.
    :type grid: GridDalesOpenBC
    :param output_path: Where the data should be output.

    Gets the initial fields for the current simulation from the host initfields.nc.
    """
    # Load data
    with xr.open_mfdataset(f"{input_json['inpath']}initfields.inp.*.nc") as ds:
        initfields_fine = ds.interp(
            xt=grid.xt,
            xm=grid.xm,
            yt=grid.yt,
            ym=grid.ym,
            zt=grid.zt,
            zm=grid.zm,
            assume_sorted=False,
        )
        initfields_fine = initfields_fine.assign_coords(
            {
                "xt": grid.xt,
                "xm": grid.xm,
                "yt": grid.yt,
                "ym": grid.ym,
            }
        )
        # Adjust transform
        # initfields_fine["transform"].attrs["false_easting"] = initfields_fine[
        #     "transform"
        # ].attrs["false_easting"]
        # initfields_fine["transform"].attrs["false_northing"] = initfields_fine[
        #     "transform"
        # ].attrs["false_northing"]
        # proj4 = ""
        # for param in initfields_fine["transform"].attrs["proj4"][1:].split("+"):
        #     line = "+" + param
        #     if "x_0" in param:
        #         line = f"+x_0={initfields_fine['transform'].attrs['false_easting']} "
        #     if "y_0" in param:
        #         line = f"+y_0={initfields_fine['transform'].attrs['false_northing']} "
        #     proj4 = proj4 + line
        # initfields_fine["transform"].attrs["proj4"] = proj4.rstrip()
        # Add global attributes
        initfields_fine = initfields_fine.assign_attrs(
            {
                "title": f"initfields.inp.{input_json['iexpnr']:03d}.nc",
                "history": f"Created on {datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')} UTC",
                "author": input_json["author"],
                "time0": input_json["time0"],
            }
        )
    initfields_fine.to_netcdf(
        path=output_path / "input" / initfields_fine.attrs["title"],
        mode="w",
        format="NETCDF4",
    )
    return initfields_fine
