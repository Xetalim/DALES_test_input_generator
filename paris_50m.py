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
import numpy as np
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
from modular_dales.Geometry import GridDales
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
    case_name = "paris_case_50m"
    output_directory = None

    grid_50m = GridDales(
        itot=1024,
        jtot=1024,
        kmax=128,
        xsize=51200.0,
        ysize=51200.0,
        kmax_soil=4,
        xlat=48.85510,
        xlon=2.34953,
        x0=626701.73,  # 626701.483,6835702.782
        y0=6835702.73,
        alpha=1.009,
        dz0=10,
        # proj4="+proj=lcc +lat_1=48.849991 +lat_0=48.849991 +lon_0=2.349999 +k_0=1 +R=6371229 +units=m +no_defs",
        proj4="EPSG:2154",  # RGF93 / Lambert-93, used by IGN for France
        # proj4="+proj=lcc +lat_1=48.85 +lat_0=48.85 +lon_0=2.35 +k_0=1 +R=6371229 +units=m +no_defs +type=crs",  # optional
    )

    grid_200m = GridDales(
        itot=1024,
        jtot=1024,
        kmax=128,
        xsize=204_800.0,
        ysize=204_800.0,
        kmax_soil=4,
        xlat=48.85510,
        xlon=2.34953,
        x0=549901.56,
        y0=6758902.73,
        alpha=1.012,
        dz0=10,
        # proj4="+proj=lcc +lat_1=48.849991 +lat_0=48.849991 +lon_0=2.349999 +k_0=1 +R=6371229 +units=m +no_defs",
        proj4="EPSG:2154",  # RGF93 / Lambert-93, used by IGN for France
        # proj4="+proj=lcc +lat_1=48.85 +lat_0=48.85 +lon_0=2.35 +k_0=1 +R=6371229 +units=m +no_defs +type=crs",  # optional
    )
    dz0 = 10
    grid = grid_200m
    print(f"{grid.zt=}")
    print(f"{grid.zm=}")
    print(f"{grid.zsize=}")
    # this is the grid that is inside us
    # Machine configuration (would normally come from machine_conf.yaml)
    # Create simulation instance
    sim = dales_simulation(case_name, machine_conf)
    logger.info(f"Created simulation with case: {case_name}")

    # Add modules - this demonstrates the += operator for module registration
    sim += DefaultNamelistModule()
    logger.info("Added DefaultNamelistModule")

    sim += grid
    print(f"{grid.zt=}")
    print(f"{grid.zm=}")
    print(f"{grid.zsize=}")
    logger.info("Added GridModule")

    atmo = AtmosphereModule()
    atmo += AtmosphericProfile(variable=ua, shape="lin", params=dict(surf_val=3, ddz=0))
    atmo += AtmosphericProfile(variable=va, shape="lin", params=dict(surf_val=3, ddz=0))
    atmo += AtmosphericProfile(variable=ug, shape="lin", params=dict(surf_val=3, ddz=0))
    atmo += AtmosphericProfile(variable=vg, shape="lin", params=dict(surf_val=3, ddz=0))
    # atmo += AtmosphericProfile(
    #     variable=thetal, shape="lin", params=dict(surf_val=293.15, ddz=1e-2)
    # )
    atmo += InterpolatedProfile(
        variable=thetal,
        z=[0, 400, 410, 3000, 4000, 5000],
        points=[293.15, 293.15, 298.15, 301.15, 301.15, 301.15],
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
    atmo += AtmosphericProfile(
        variable=w,
        shape="lin",
        params=dict(surf_val=0.0, ddz=0.0),
    )
    sim += atmo
    sim += TimeModule(
        xday=232,
        xtime=0.0,
        xyear=2023,
        runtime=3600,
        startyear=2023,
        startmonth=8,
        startday=20,
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

    lsm += LSM.LandUseModification(geometry="all", type="grs", params={})
    soil_levels = [0.01, 0.04, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.5, 2, 3, 5, 8, 12]
    lsm += LSM.UniformSkinTemperature(
        skin_temperature=293.15,  # 20°C in Kelvin
    )
    lsm += LSM.UniformSoilMoisture([0.2, 0.2, 0.2, 0.2])
    lsm += LSM.UniformSoilTemperature([293.15, 293.15, 293.15, 293.15])
    # lsm += LSM.SoilTemperatureMoistureFromHarmonie(
    #     harmonie_soil_file="/ec/res4/scratch/nld4411/dales_nest_harmonie/paris_20_april/data/GNATU.nc",
    #     harmonie_soil_valid_time="2023-08-20T00:00:00",
    #     harmonie_soil_height_levels=soil_levels,
    #     use_as_tskin=True,
    # )
    lsm += LSM.FromLCZ()
    sim += lsm

    slurb = SLURBModule(deep_soil_temperature=283)
    # slurb += SLURBModification(
    #     geometry="all", vars=[{"varname": "albedo_av", "value": 10}], params={}
    # )
    sim += slurb
    # sim += nesting

    # in_harm = Nest_in_Harmonie(ml_glob="/ec/res4/scratch/nld4411/dales_nest_harmonie/paris_20_april/data/nc_out/ml*.nc",sfc_glob="/ec/res4/scratch/nld4411/dales_nest_harmonie/paris_20_april/data/nc_out/sfc*.nc")

    # openbc = do_openboundary(
    #     time0="2023-08-20T00:00:00",
    #     start="2023-08-20T0:00:00",
    #     end="2023-08-021T23:59:00",
    #     e12=1,
    #     tauh=0,
    #     taum=100,
    #     dxint=supergrid.xsize / supergrid.itot * 4,
    #     dyint=supergrid.ysize / supergrid.jtot * 4,
    # )

    # openbc += in_harm

    # sim += openbc

    sim += EasyOutputModule(
        output_interval=120,
        enable_output=True,
    )
    sim += CheckSimulationModule(
        check_interval=120, stop_on_invalid=False, check_tendencies=False
    )
    set_nml_section(sim.nml, sim.nml_docs, "no_coriolis", "physics", "lcoriol", True)
    sim.sim_preprocessing_pipeline()
