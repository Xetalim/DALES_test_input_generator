import helper_scripts.LBC.prep_harmonie as prep_harmonie
import helper_scripts.LBC.boundary as boundary
import helper_scripts.LBC.initfields as initfields
import helper_scripts.LBC.initfields as initfields
import helper_scripts.LBC.dales_nest.boundary_fields_fine as boundary_fields_fine
import helper_scripts.LBC.dales_nest.get_timesteps0 as get_timesteps0
import helper_scripts.LBC.dales_nest.initial_fields_fine as initial_fields_fine
from helper_scripts.grids import GridDalesOpenBC, GridDales, nesting_idx
import json
import yaml
import numpy as np
import dask
from dask.diagnostics import ProgressBar
import logging
from helper_scripts.logging_wrapper import logwrap

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

            self.indices = nesting_idx(
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
            self.boundaries = boundary.boundary_fields(
                config["openboundary"],
                self.openBCgrid,
                self.data,
                output_path=output_path,
            )
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
        logger.debug("Writing initial fields")
        self.initfields.to_netcdf(
            path=output_path / "input" / self.initfields.attrs["title"],
            mode="w",
            format="NETCDF4",
            compute=True
        )
        logger.debug("Writing boundaries")
        self.boundaries.to_netcdf(
            path=output_path / "input" / self.boundaries.attrs["title"],
            mode="w",
            format="NETCDF4",
            compute=True
        )
