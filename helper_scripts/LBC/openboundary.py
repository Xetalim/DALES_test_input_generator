import helper_scripts.LBC.prep_harmonie as prep_harmonie
import helper_scripts.LBC.boundary as boundary
import helper_scripts.LBC.initfields as initfields
import helper_scripts.LBC.initfields as initfields
import helper_scripts.LBC.dales_nest.boundary_fields_fine as boundary_fields_fine
import helper_scripts.LBC.dales_nest.get_timesteps0 as get_timesteps0
import helper_scripts.LBC.dales_nest.initial_fields_fine as initial_fields_fine
import helper_scripts.grids as grids
import json
import yaml
import logging

logger = logging.getLogger(__name__)


def do_openboundary(grid, config, output_path, nml, exp_id):
    if config["openboundary"]["source"] == "HARMONIE":
        data, transform = prep_harmonie.prep_harmonie(config["openboundary"], grid)
        # we need to use the right surface pressure as calculated from the input data
        logger.info("Setting namelist PHYSICS:ps to %f", config["openboundary"]["ps"])
        nml["PHYSICS"]["ps"] = config["openboundary"]["ps"]
        # we need to use the right skin liq. pot. temperature from the input data
        logger.info(
            "Setting namelist PHYSICS:thls to %f", config["openboundary"]["thls"]
        )
        nml["PHYSICS"]["thls"] = config["openboundary"]["thls"]
        boundary.boundary_fields(
            config["openboundary"], grid, data, output_path=output_path
        )
        initfields.initial_fields(
            config["openboundary"], grid, data, transform, output_path=output_path
        )
    elif config["openboundary"]["source"] == "DALES":
        data = initial_fields_fine.initial_fields_fine(
            config["openboundary"], grid=grid, output_path=output_path
        )
    else:
        raise ValueError("Unknown source for open boundary conditions")
