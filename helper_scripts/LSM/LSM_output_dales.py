# import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
import logging
from helper_scripts.logging_wrapper import logwrap
import matplotlib.pyplot as plt
import numpy as np
import os

from helper_scripts.grids import generate_dales_domain
from helper_scripts.grids import GridDales

# Custom Python scripts/tools/...
from helper_scripts.LSM.vegetation_properties import ifs_vegetation, top10_to_ifs
import logging
from helper_scripts.logging_wrapper import logwrap

logger = logging.getLogger(__name__)
logger.debug("Entered module: %s", __name__)

# Correction factor for aspect ratio of plots
ASPECT_CORR = 2


# @logwrap


#     return lsm_input, lu_types


@logwrap
def some_plots(lsm_input, plotvars, output_path):
    """
    Generate some standard plots of Land Surface Model input data

    Parameters
    ----------
    lsm_input : LSM_input_DALES
        Class containing Dales input parameters for all LU types.
    plotvars : list
        List of variables to plot.

    Returns
    -------
    None.

    """
    import xarray as xr

    data_vars = dict()
    os.makedirs(output_path / ".." / "covers", exist_ok=True)
    for plotvar in plotvars:
        data = getattr(lsm_input, plotvar)
        data_vars[plotvar] = (["y", "x"], data)

    coords = dict(x=(["x"], lsm_input.x), y=(["y"], lsm_input.y))
    ds_lsm = xr.Dataset(data_vars=data_vars, coords=coords)

    for plotvar in list(ds_lsm.variables):
        if "cover" in plotvar:
            vmax = 1
        else:
            vmax = None
        if plotvar == "x" or plotvar == "y":
            continue
        fig, ax = plt.subplots(1)
        ax.set_aspect(
            abs((lsm_input.y[-1] - lsm_input.y[0]) / (lsm_input.x[-1] - lsm_input.x[0]))
            * ASPECT_CORR
        )
        ds_lsm[plotvar].plot(ax=ax, vmin=0, vmax=vmax)
        plt.tight_layout()

        plt.savefig(output_path / ".." / "covers" / f"var_{plotvar}.png", dpi=300)
        plt.close()
    return


class LSM_output_dales:
    """
    Data structure for the required input for the new LSM
    """

    @logwrap
    def __init__(
        self,
        grid,
        lu_types,
        soil_levels=4,
        debug=False,
    ):

        # land use parameters
        self.parnames = [
            "cover",
            "c_veg",
            "z0m",
            "z0h",
            "lai",
            "ar",
            "br",
            "lambda_s",
            "lambda_us",
            "rs_min",
            "gD",
            "tskin",
            "lutype",
        ]
        dtype_float = np.float64
        dtype_int = np.int32
        dtype_str = object
        self.nlu = len(lu_types)
        self.soil_levels = soil_levels
        self.grid = grid

        self.value_dic = {
            "t_soil": np.zeros((soil_levels, grid.jtot, grid.itot), dtype=dtype_float),
            "index_soil": np.zeros(
                (soil_levels, grid.jtot, grid.itot), dtype=dtype_int
            ),
            "theta_soil": np.zeros(
                (soil_levels, grid.jtot, grid.itot), dtype=dtype_float
            ),
        }

        # Soil temperature, moisture content, and index in van Genuchten lookup table.
        fields = ["index_soil", "t_soil", "theta_soil"]

        # LU types
        self.luname = np.empty(len(lu_types), dtype=dtype_str)
        self.lushort = np.empty(len(lu_types), dtype=dtype_str)
        self.lveg = np.empty(len(lu_types), dtype=dtype_str)
        self.laqu = np.empty(len(lu_types), dtype=dtype_str)
        self.ilu = np.zeros(len(lu_types), dtype=dtype_int)

        # Sub-grid fraction of LU type (-)
        zeros = np.zeros((grid.jtot, grid.itot), dtype=dtype_float)
        for lu in lu_types:
            for parname in self.parnames:
                varname = "_".join([parname, lu])
                self.value_dic[varname] = zeros.copy()
                fields.append(varname)
        # total land use cover
        self.c_tot = np.zeros((grid.jtot, grid.itot), dtype=dtype_float)
        fields.append("cover_tot")
        self.c_veg_tot = np.zeros((grid.jtot, grid.itot), dtype=dtype_float)
        fields.append("c_veg_tot")

        # List of fields which are written to the binary input files for DALES
        self.fields = sorted(fields)

        # Bonus, for offline LSM (not written to DALES input)
        for lu in lu_types:
            ## LU type (-)
            varname = "type_" + lu
            self.value_dic[varname] = zeros.copy()

        if debug:
            # Init all values at NaN
            for field in self.fields:
                data = self.value_dic[field]
                if data.dtype == dtype_int:
                    continue
                data[:] = np.nan

        self.setup_lu_types(lu_types)

    def standard_fill_geometry_modification(self, modify_func=None):
        self.value_dic["index_soil"][:, :] = 2

        # represents average over index_soil=2
        for i in range(4):
            self.value_dic["theta_soil"][i, :, :] = [
                0.36867549,
                0.25300502,
                0.14997292,
                0.16459982,
            ][i]
            self.value_dic["t_soil"][i, :, :] = [
                283.26038083,
                286.79894009,
                290.88998902,
                288.09079126,
            ][i]

        # Python -> Fortran indexing
        self.value_dic["index_soil"] += 1

        if modify_func:
            modify_func(self)

        self.init_lutypes_ifs()

    def setup_lu_types(self, lu_types):
        """
        Set up arrays for luname, lushort, lveg, laqu, ilu. Doesn't set any values yet.

        """

        luname = [value["lu_long"] for key, value in lu_types.items()]
        lushort = [value["lu_short"] for key, value in lu_types.items()]
        lveg = [value["lveg"] for key, value in lu_types.items()]
        laqu = [value["laqu"] for key, value in lu_types.items()]
        ilu = np.arange(len(lu_types)) + 1

        self.luname = luname
        self.lushort = lushort
        self.lveg = lveg
        self.laqu = laqu
        self.ilu = ilu

        # set LU cover for each grid cell
        for lu in lu_types:
            lu_types[lu]["lu_domid"] = (
                np.ones((self.grid.jtot, self.grid.itot)) * lu_types[lu]["lu_ids"][0]
            )
            lu_types[lu]["lu_frac"] = np.ones((self.grid.jtot, self.grid.itot)) * 0
            self.value_dic["c_" + lu] = lu_types[lu]["lu_frac"]
        self.lu_types = lu_types

    @logwrap
    def save_netcdf(
        self,
        output_path,
        exp_id,
    ):
        """
        Save to NetCDF for visualisation et al.
        """
        nc = nc4.Dataset(f"{output_path}/lsm.inp_{exp_id:03d}.nc", "w")

        nc.createDimension("x", self.grid.itot)
        nc.createDimension("y", self.grid.jtot)
        nc.createDimension("z", self.soil_levels)
        nc.createDimension("nlu", self.nlu)
        nc.createDimension("str3", size=3)
        nc.createDimension("str32", size=32)
        nc.createDimension("str1", size=1)

        var_x = nc.createVariable("x", float, "x")
        var_y = nc.createVariable("y", float, "y")

        var_x[:] = self.grid.xt[:]
        var_y[:] = self.grid.yt[:]

        for field in self.fields:
            data = self.value_dic[field]
            dims = ["y", "x"] if data.ndim == 2 else ["z", "y", "x"]
            var = nc.createVariable(field, float, dims)
            var[:] = data[:]

        luname = nc4.stringtochar(np.array(self.luname, "S32"))
        var_lun = nc.createVariable(
            "luname", datatype="S1", dimensions=("nlu", "str32")
        )
        var_lun[:, :] = luname

        lushort = nc4.stringtochar(np.array(self.lushort, "S3"))
        var_lus = nc.createVariable(
            "lushort", datatype="S1", dimensions=("nlu", "str3")
        )
        var_lus[:, :] = lushort

        lveg = nc4.stringtochar(np.array([str(b)[0] for b in self.lveg], "S1"))
        var_lveg = nc.createVariable("lveg", datatype="S1", dimensions=("nlu", "str1"))
        var_lveg[:, :] = lveg

        laqu = nc4.stringtochar(np.array([str(b)[0] for b in self.laqu], "S1"))
        var_laqu = nc.createVariable("laqu", datatype="S1", dimensions=("nlu", "str1"))
        var_laqu[:, :] = laqu

        ilu = np.array(self.ilu, dtype=int)
        var_ilu = nc.createVariable("ilu", int, ("nlu",))
        var_ilu[:] = ilu

        nc.close()

        return

    @logwrap
    def init_lutypes_ifs(self):
        """
        Given some pre-existing land covers, gets properties from the IFS lookup table... TODO UNDERSTAND


        """
        shape = (self.grid.jtot, self.grid.itot)
        for lu in self.lu_types.keys():
            for parname in self.parnames:
                if parname == "cover" or parname == "c_veg":
                    parfield = self.lu_types[lu]["lu_frac"].copy()
                else:
                    parfield = np.full(shape, 0.0)

                for vt in self.lu_types[lu]["lu_ids"]:
                    iv = top10_to_ifs[vt]  # Index in ECMWF lookup table

                    if parname == "lutype":
                        parfield[:] = (
                            iv  # LG: Only apply mask to cover and c_veg (DALES crashes when zeros or nans are in
                        )
                        # the array)
                    elif parname == "tskin":
                        parfield[:] = (
                            273.15  # LG: Only apply mask to cover and c_veg (DALES crashes when zeros or nans
                        )
                        # are in the array)
                    elif parname in ["cover", "c_veg"]:
                        continue
                    else:
                        if parname == "ar":
                            parname_ifs = "a_r"
                        elif parname == "br":
                            parname_ifs = "b_r"
                        else:
                            parname_ifs = parname
                        # parfield[mask] = getattr(ifs_vegetation, parname_ifs) [iv]
                        parfield[:] = getattr(ifs_vegetation, parname_ifs)[
                            iv
                        ]  # LG: Only apply mask to cover and c_veg
                self.value_dic[f"{parname}_{lu}"] = parfield

        totcover = self.calc_totcover("cover")
        self.value_dic["cover_tot"] = totcover

        totcveg = self.calc_totcover("c_veg")
        self.value_dic["c_veg_tot"] = totcveg

        # TODO: more consistent way to check for LU type with bare soil
        bs_name = [
            k
            for k in self.lu_types.keys()
            if "bar" in self.lu_types[k]["lu_long"].lower()
        ][0]

        cover = self.value_dic["cover_tot"]
        cover_bs0 = self.value_dic[f"cover_{bs_name}"]
        cover_bs = np.round(1 - cover + cover_bs0, 6)

        self.value_dic["cover_" + bs_name] = cover_bs

        # recalculate
        totcover = self.calc_totcover("cover")
        self.value_dic["cover_tot"] = totcover

        totcveg = self.calc_totcover("c_veg")
        self.value_dic["c_veg_tot"] = totcveg

    @logwrap
    def calc_totcover(self, ctype):
        """
        Calculate sum over cover of individual LU types to check if it sums up to 1

        Parameters
        ----------
        lsm_input : LSM_input_DALES
            Class containing Dales input parameters for all LU types.
        lu_types : dict
            LU type properties.
        ctype : str
            LU cover type to be summed.

        Returns
        -------
        totcover : np.array
            Total LU cover.

        """
        covers = [ctype + "_" + s for s in self.lu_types.keys()]
        totcover = np.zeros([self.grid.jtot, self.grid.itot])
        for c in covers:
            totcover += self.value_dic[c]

        return totcover
