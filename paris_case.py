"""Example demonstrating complete DALES simulation module usage.

This example shows how to:
1. Create a simulation instance with configuration
2. Add all available modules
3. Configure, validate, prepare, and write output

The modular system uses dataclasses where:
- Each module is a dataclass with optional `sim` parent reference
- Modules expose properties to access sim.config, sim.grid, sim.output_path, sim.nml
- Namelist values are declared via field metadata: metadata={"nml": "SECTION", "key": "name"}
- Runtime-only fields use: field(default=None, init=False, repr=False)
"""

import logging
import os
import subprocess
import pytest
import yaml

from modular_dales.Atmosphere import (
    AtmosphereModule,
    AtmosphericProfile,
    InterpolatedProfile,
)
from modular_dales.Configuration import (
    DefaultNamelistModule,
    EasyOutputModule,
    TimeModule,
    SamplingModule,
)
from modular_dales.Configuration.output_modules import CheckSimulationModule
from modular_dales.Geometry import GridDales, AllGeometry
from modular_dales.LBC import Nest_in_Dales, NestingTopology, do_openboundary
from modular_dales.Radiation.radiation import RadiationModule
import modular_dales.Surface as Surface
import modular_dales.Surface.LSM as LSM
from modular_dales.Surface.LSM.SLuRB.slurb import SLURBModification, SLURBModule
from modular_dales.logging_wrapper import setup_logging
from modular_dales.modular import dales_simulation
from modular_dales.Surface import ConstantSurfaceTemperatureModule
from modular_dales.modular.simulation_module import set_nml_section
from modular_dales.vars import *

setup_logging("logging.yaml")
# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# logger.info("Added TimeModule")

# # Execute workflow
# logger.info("\n--- Initializing simulation ---")
# sim.init_output_folder()
# sim.setup_module_links()
# logger.info("--- Configuring modules ---")
# sim.do_config()

# logger.info("--- Checking settings ---")
# sim.check_settings()

# logger.info("--- Preparing calculations ---")
# sim.prepare_all_calculations()
# sim.apply_job_configuration()
# logger.info("--- Writing files ---")
# sim.write_module_files()
# sim.write_simulation_files()

# logger.info("\nBasic simulation completed successfully!")


if __name__ == "__main__":
    # pytest.skip("unsupported case")
    with open("machine_conf.yaml", "r") as file:
        machine_conf = yaml.safe_load(file)
    """Create a basic DALES simulation with minimal configuration.

    This example demonstrates:
    - Creating a simulation with config and machine settings
    - Adding default and grid modules
    - Executing the simulation workflow
    """
    logger.info("=" * 70)
    logger.info("EXAMPLE 1: Basic Simulation Setup")
    logger.info("=" * 70)

    # Create minimal configuration with just case name and output directory
    case_name = "011_test"
    output_directory = None

    supergrid = GridDales(
        itot=512,
        jtot=512,
        kmax=80,
        xsize=25600.0,
        ysize=25600.0,
        kmax_soil=4,
        xlat=48.85510,
        xlon=2.34953,
        x0=-11296,
        y0=-5979,
        alpha=1.0,
        dz0=5,
        # proj4="+proj=lcc +lat_1=48.85 +lat_0=48.85 +lon_0=2.35 +k_0=1 +R=6371229 +units=m +no_defs +type=crs",  # optional
        wkt="""PROJCRS["unnamed",
    BASEGEOGCRS["unknown",
        DATUM["unnamed",
            ELLIPSOID["Sphere",6371229,0,
                LENGTHUNIT["metre",1,
                    ID["EPSG",9001]]]],
        PRIMEM["Greenwich",0,
            ANGLEUNIT["degree",0.0174532925199433,
                ID["EPSG",9122]]]],
    CONVERSION["unnamed",
        METHOD["Lambert Conic Conformal (1SP)",
            ID["EPSG",9801]],
        PARAMETER["Latitude of natural origin",48.85,
            ANGLEUNIT["degree",0.0174532925199433],
            ID["EPSG",8801]],
        PARAMETER["Longitude of natural origin",2.35,
            ANGLEUNIT["degree",0.0174532925199433],
            ID["EPSG",8802]],
        PARAMETER["Scale factor at natural origin",1,
            SCALEUNIT["unity",1],
            ID["EPSG",8805]],
        PARAMETER["False easting",0.0,
            LENGTHUNIT["metre",1],
            ID["EPSG",8806]],
        PARAMETER["False northing",0.0,
            LENGTHUNIT["metre",1],
            ID["EPSG",8807]],
        PARAMETER["standard_parallel_1",48.85,
            ANGLEUNIT["degree",0.0174532925199433,
                ID["EPSG",9122]]]],
    CS[Cartesian,2],
        AXIS["easting",east,
            ORDER[1],
            LENGTHUNIT["metre",1,
                ID["EPSG",9001]]],
        AXIS["northing",north,
            ORDER[2],
            LENGTHUNIT["metre",1,
                ID["EPSG",9001]]]]""",
    )
    # this is the grid that is inside us
    # subgrid = GridDales(
    #     itot=512,
    #     jtot=512,
    #     kmax=70,
    #     xsize=6400.0,
    #     ysize=6400.0,
    #     kmax_soil=4,
    #     xlat=52.25,
    #     xlon=5.45,
    #     x0=160.0,
    #     y0=160.0,
    #     alpha=1.0,
    #     dz0=5,
    #     # "proj4": "+proj=..."  # optional
    # )

    nesting = NestingTopology()
    # this is us
    nesting += supergrid
    nesting.my_idx = nesting.nestings.index(supergrid)

    # nesting += subgrid

    # Machine configuration (would normally come from machine_conf.yaml)
    # Create simulation instance
    sim = dales_simulation(case_name, machine_conf)
    logger.info(f"Created simulation with case: {case_name}")

    # Add modules - this demonstrates the += operator for module registration
    sim += DefaultNamelistModule()
    logger.info("Added DefaultNamelistModule")

    sim += supergrid
    logger.info("Added GridModule")

    atmo = AtmosphereModule()
    atmo += AtmosphericProfile(variable=ua, shape="lin", params=dict(surf_val=0, ddz=0))
    atmo += AtmosphericProfile(variable=ug, shape="lin", params=dict(surf_val=0, ddz=0))
    atmo += AtmosphericProfile(
        variable=va, shape="lin", params=dict(surf_val=3, ddz=1e-3)
    )
    atmo += AtmosphericProfile(
        variable=vg, shape="lin", params=dict(surf_val=3, ddz=1e-3)
    )
    # atmo += AtmosphericProfile(
    #     variable=thetal, shape="lin", params=dict(surf_val=293.15, ddz=1e-2)
    # )
    atmo += InterpolatedProfile(
        variable=thetal,
        z=[0, 200, 210, 500],
        points=[293.15, 293.15, 298.15, 301.15],
    )
    atmo += AtmosphericProfile(
        variable=qt, shape="lin", params=dict(surf_val=0.0018, ddz=0)
    )
    atmo += AtmosphericProfile(variable=wa, shape="lin", params=dict(surf_val=0, ddz=0))
    atmo += InterpolatedProfile(
        variable=tke,
        z=[0, 4000, 5000],
        points=[1, 1e-8, 1e-8],
    )

    sim += atmo
    sim += TimeModule(
        xday=1,
        xtime=12.0,
        xyear=2023,
        runtime=1300,
        startyear=2023,
        startmonth=1,
        startday=1,
    )

    sim += RadiationModule(
        iradiation=4,  # 0=off, 1=simple, 4=RRTMGP
    )
    lsm = LSM.LSMModule(
        ps=100000,
        z0mav=0.0001,
        z0hav=0.0001,
        iinterp_t=1,
        iinterp_theta=1,
        dz_soil=[1.89, 0.72, 0.21, 0.07],
        albedoav=0.22,
    )

    lsm += LSM.LandUseModification(geometry=AllGeometry(), type="grs")
    soil_levels = [0.01, 0.04, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.5, 2, 3, 5, 8, 12]
    lsm += LSM.SoilTemperatureMoistureFromHarmonie(
        harmonie_soil_file="/Users/andrevanginkel/tmp/soil_split/diditwork.nc",
        harmonie_soil_valid_time="2023-08-20T12:00:00",
        harmonie_soil_height_levels=soil_levels,
        use_as_tskin=True,
    )
    lsm += LSM.FromLCZ()
    sim += lsm

    slurb = SLURBModule(deep_soil_temperature=293)
    # slurb += SLURBModification(
    #     geometry=AllGeometry(), vars=[{"varname": "albedo_av", "value": 10}]
    # )
    sim += slurb
    sim += nesting

    sim += EasyOutputModule(
        output_interval=60,
        enable_output=True,
    )
    sim += CheckSimulationModule(
        check_interval=60, stop_on_invalid=False, check_tendencies=False
    )
    set_nml_section(sim.nml, sim.nml_docs, "no_coriolis", "physics", "lcoriol", False)
    sim.sim_preprocessing_pipeline()
