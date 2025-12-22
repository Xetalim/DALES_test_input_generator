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
import logging
from helper_scripts.logging_wrapper import logwrap

logger = logging.getLogger(__name__)
logger.debug("Entered module: %s", __name__)


@logwrap
def do_openboundary(grid: GridDalesOpenBC, config, output_path, nml, exp_id):
    if not type(grid) is GridDalesOpenBC:
        logger.warning(
            "Passing wrong grid type to openboundary, converting to GridDalesOpenBC with as_openBC."
        )
        openBCgrid = grid.as_openbc()
    else:
        openBCgrid = grid

    if "subnest" in config["openboundary"]:
        # we want to nest a model in this, we can provide some convenience function for the boundary indices.
        subgrid = GridDalesOpenBC(config["openboundary"]["subnest"])
        supergrid = grid

        ixwest, ixeast = list(supergrid.xm).index(np.min(subgrid.xm)), list(
            supergrid.xm
        ).index(np.max(subgrid.xm))

        iysouth, iynorth = list(supergrid.ym).index(np.min(subgrid.ym)), list(
            supergrid.ym
        ).index(np.max(subgrid.ym))

        iztop = list(supergrid.zt).index(np.max(subgrid.zt))

        # integer :: crossplane(100) !< Location of the xz crosssection
        # integer :: crossortho(100) !< Location of the yz crosssection

        nml["namcrosssection"]["crossheight"].append(iztop + 1)
        nml["namcrosssection"]["crossplane"].append(iysouth + 1)
        nml["namcrosssection"]["crossplane"].append(iynorth + 1)
        nml["namcrosssection"]["crossortho"].append(ixwest + 1)
        nml["namcrosssection"]["crossortho"].append(ixeast + 1)
    if "supernest" in config["openboundary"]:
        subgrid = grid.as_openbc()
        supergrid = GridDales(config["openboundary"]["supernest"])

        ixwest, ixeast = list(supergrid.xm).index(np.min(subgrid.xm)), list(
            supergrid.xm
        ).index(np.max(subgrid.xm))

        iysouth, iynorth = list(supergrid.ym).index(np.min(subgrid.ym)), list(
            supergrid.ym
        ).index(np.max(subgrid.ym))

        iztop = list(supergrid.zt).index(np.max(subgrid.zt))

        indices = nesting_idx(
            ix_west=ixwest,
            ix_east=ixeast,
            iy_south=iysouth,
            iy_north=iynorth,
            subgrid_x0=subgrid.x0,
            subgrid_y0=subgrid.y0,
            supergrid_x0=supergrid.x0,
            supergrid_y0=supergrid.y0,
        )
    else:
        indices = None

    if config["openboundary"]["source"] == "HARMONIE":
        data, transform = prep_harmonie.prep_harmonie(
            config["openboundary"], openBCgrid
        )
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
        data, = dask.optimize(data)
        boundaries_writer = boundary.boundary_fields(
            config["openboundary"], openBCgrid, data, output_path=output_path
        )
        initfields_writer = initfields.initial_fields(
            config["openboundary"], openBCgrid, data, transform, output_path=output_path
        )
        logger.debug("Setup all openBC writers, optimizing writers now..")
        boundaries_writer, initfields_writer = dask.optimize(boundaries_writer, initfields_writer)
        logger.debug("Optimized writers")
        logger.debug("Writing initial fields")
        initfields_writer.compute()
        logger.debug("Writing boundaries")
        boundaries_writer.compute()
    elif config["openboundary"]["source"] == "DALES":
        data = initial_fields_fine.initial_fields_fine(
            config["openboundary"], grid=openBCgrid, output_path=output_path
        )

        boundary_fields_fine.boundary_fields_fine(
            config["openboundary"],
            grid=openBCgrid,
            output_path=output_path,
            grid_indices=indices,
        )
    else:
        raise ValueError("Unknown source for open boundary conditions")
