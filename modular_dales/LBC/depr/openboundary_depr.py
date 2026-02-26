from modular_dales.LBC.nesting_idx import NestingIndices
import modular_dales.LBC.nest_dales_in_HARMONIE.prep_harmonie as prep_harmonie
import modular_dales.LBC.nest_dales_in_HARMONIE.boundary as boundary
import modular_dales.LBC.nest_dales_in_HARMONIE.initfields as initfields
import modular_dales.LBC.nest_dales_in_HARMONIE.initfields as initfields
import modular_dales.LBC.nest_dales_in_dales.boundary_fields_fine as boundary_fields_fine
import modular_dales.LBC.nest_dales_in_dales.get_timesteps0 as get_timesteps0
import modular_dales.LBC.nest_dales_in_dales.initial_fields_fine as initial_fields_fine
from modular_dales.Geometry.GridDales import GridDalesOpenBC, GridDales
import json
import yaml
import numpy as np
import dask
from dask.diagnostics import ProgressBar
import logging
from modular_dales.logging_wrapper import logwrap

logger = logging.getLogger(__name__)
logger.debug("Entered module: %s", __name__)


@logwrap
class do_openboundary:
    def __init__(self, grid: GridDalesOpenBC):
        self.grid = grid
        if not type(grid) is GridDalesOpenBC:
            logger.warning(
                "Passing wrong grid type to openboundary, converting to GridDalesOpenBC with as_openBC."
            )
            self.openBCgrid = grid.as_openbc()
        else:
            self.openBCgrid = grid

    def setup(self, config, output_path, nml, exp_id):
        if "subnest" in config["openboundary"]:
            # we want to nest a model in this, we can provide some convenience function for the boundary indices.
            self.subnest_subgrid = GridDalesOpenBC(config["openboundary"]["subnest"])
            self.subnest_supergrid = self.grid

            ixwest, ixeast = list(self.subnest_supergrid.xm).index(
                np.min(self.subnest_subgrid.xm)
            ), list(self.subnest_supergrid.xm).index(np.max(self.subnest_subgrid.xm))

            iysouth, iynorth = list(self.subnest_supergrid.ym).index(
                np.min(self.subnest_subgrid.ym)
            ), list(self.subnest_supergrid.ym).index(np.max(self.subnest_subgrid.ym))

            iztop = list(self.subnest_supergrid.zt).index(
                np.max(self.subnest_subgrid.zt)
            )

            # integer :: crossplane(100) !< Location of the xz crosssection
            # integer :: crossortho(100) !< Location of the yz crosssection

            nml["namcrosssection"]["crossheight"].append(iztop + 1)
            nml["namcrosssection"]["crossplane"].append(iysouth + 1)
            nml["namcrosssection"]["crossplane"].append(iynorth + 1)
            nml["namcrosssection"]["crossortho"].append(ixwest + 1)
            nml["namcrosssection"]["crossortho"].append(ixeast + 1)
        if "supernest" in config["openboundary"]:
            self.supernest_subgrid = self.grid.as_openbc()
            self.supernest_supergrid = GridDales(config["openboundary"]["supernest"])

            ixwest, ixeast = list(self.supernest_supergrid.xm).index(
                np.min(self.supernest_subgrid.xm)
            ), list(self.supernest_supergrid.xm).index(
                np.max(self.supernest_subgrid.xm)
            )

            iysouth, iynorth = list(self.supernest_supergrid.ym).index(
                np.min(self.supernest_subgrid.ym)
            ), list(self.supernest_supergrid.ym).index(
                np.max(self.supernest_subgrid.ym)
            )

            iztop = list(self.supernest_supergrid.zt).index(
                np.max(self.supernest_subgrid.zt)
            )

            self.indices = NestingIndices(
                ix_west=ixwest,
                ix_east=ixeast,
                iy_south=iysouth,
                iy_north=iynorth,
                subgrid_x0=self.supernest_subgrid.x0,
                subgrid_y0=self.supernest_subgrid.y0,
                supergrid_x0=self.supernest_supergrid.x0,
                supergrid_y0=self.supernest_supergrid.y0,
            )
        else:
            self.indices = None

    def prepare_calculation(self, config, output_path, nml, exp_id):

        if config["openboundary"]["source"] == "HARMONIE":
            self.harmonieprepper = prep_harmonie.harmoniePrepper(
                config["openboundary"], self.openBCgrid
            )
            self.harmonieprepper.load_data()
            self.data, self.transform = self.harmonieprepper.prep_harmonie()
            # we need to use the right surface pressure as calculated from the input data
            logger.info(
                "Setting namelist NAMSURFACE:ps to %f", config["openboundary"]["ps"]
            )
            nml["NAMSURFACE"]["ps"] = config["openboundary"]["ps"]
            # we need to use the right skin liq. pot. temperature from the input data
            logger.info(
                "Setting namelist NAMSURFACE:thls to %f", config["openboundary"]["thls"]
            )
            nml["NAMSURFACE"]["thls"] = config["openboundary"]["thls"]
            (self.data,) = dask.optimize(self.data)

            # if True:
            #     logger.debug(f"Saving processed dataset")
            #     self.data = self.data.compute()
            #     self.data.to_netcdf("/ec/res4/scratch/nld4411/dales_nest_harmonie/netcdfs_newnew4/data_processed.nc",engine="netcdf4")
            #     del self.data

            # self.data = xr.open_dataset("/ec/res4/scratch/nld4411/dales_nest_harmonie/netcdfs_newnew4/data_processed.nc",engine="netcdf4",chunks={"time":input_json["tchunk"]})
            # logger.debug(f"Read in saved processed dataset")

            logger.debug("Setting up boundary fields")
            self.boundaries = boundary.boundary_fields(
                config["openboundary"],
                self.openBCgrid,
                self.data,
                output_path=output_path,
            )
            logger.debug("Setting up initial fields")
            self.initfields = initfields.initial_fields(
                config["openboundary"],
                self.openBCgrid,
                self.data,
                self.transform,
                output_path=output_path,
            )
            logger.debug("Setup all openBC fields, optimizing fields now..")
            self.boundaries, self.initfields = dask.optimize(
                self.boundaries, self.initfields
            )
            logger.debug("Optimized fields")

        elif config["openboundary"]["source"] == "DALES":
            self.data = initial_fields_fine.initial_fields_fine(
                config["openboundary"], grid=self.openBCgrid, output_path=output_path
            )

            boundary_fields_fine.boundary_fields_fine(
                config["openboundary"],
                grid=self.openBCgrid,
                output_path=output_path,
                grid_indices=self.indices,
            )
        else:
            raise ValueError("Unknown source for open boundary conditions")

    def write_openbcs(self, output_path):
        # Save data
        logger.debug("Optimizing initial fields")
        initfields_writer = self.initfields.to_netcdf(
            path=output_path / "input" / self.initfields.attrs["title"],
            mode="w",
            format="netcdf4",
            compute=False,
        )
        (initfields_writer,) = dask.optimize(initfields_writer)
        logger.debug("Writing initial fields")

        initfields_writer.compute()

        logger.debug("Optimizing boundaries")
        openboundaries_writer = self.boundaries.to_netcdf(
            path=output_path / "input" / self.boundaries.attrs["title"],
            mode="w",
            format="netcdf4",
            compute=False,
        )
        (openboundaries_writer,) = dask.optimize(openboundaries_writer)
        logger.debug("Writing boundaries")

        openboundaries_writer.compute()
