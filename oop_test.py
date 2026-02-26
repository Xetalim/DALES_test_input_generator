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
)
from modular_dales.Geometry import GridDales
from modular_dales.LBC import Nest_in_Dales, NestingTopology, do_openboundary
from modular_dales.logging_wrapper import setup_logging
from modular_dales.modular import dales_simulation
from modular_dales.Surface import ConstantSurfaceTemperatureModule
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
    raise NotImplementedError(
        "This is a test file for development, not meant to be run as a test or example"
    )
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
        itot=64,
        jtot=64,
        kmax=80,
        xsize=640.0,
        ysize=640.0,
        kmax_soil=4,
        xlat=52.25,
        xlon=5.45,
        x0=0,
        y0=0,
        alpha=1.0,
        dz0=5,
        # "proj4": "+proj=..."  # optional
    )

    # this is the grid that is inside us
    subgrid = GridDales(
        itot=64,
        jtot=64,
        kmax=70,
        xsize=320.0,
        ysize=320.0,
        kmax_soil=4,
        xlat=52.25,
        xlon=5.45,
        x0=160.0,
        y0=160.0,
        alpha=1.0,
        dz0=5,
        # "proj4": "+proj=..."  # optional
    )

    nesting = NestingTopology()
    # this is us
    nesting += supergrid
    nesting.my_idx = nesting.nestings.index(supergrid)

    nesting += subgrid

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
    atmo += AtmosphericProfile(
        variable=ua, shape="lin", params=dict(surf_val=0.5, ddz=0)
    )
    atmo += AtmosphericProfile(
        variable=va, shape="lin", params=dict(surf_val=3, ddz=1e-3)
    )
    atmo += AtmosphericProfile(
        variable=thetal, shape="lin", params=dict(surf_val=293.15, ddz=1e-2)
    )
    atmo += AtmosphericProfile(variable=qt, shape="lin", params=dict(surf_val=0, ddz=0))
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
        runtime=600,
        startyear=2023,
        startmonth=1,
        startday=1,
    )

    # sim += RadiationModule(
    #     iradiation=4,  # 0=off, 1=simple, 4=RRTMGP
    # )
    sim += ConstantSurfaceTemperatureModule(
        thls=293.15, z0mav=0.0001, z0hav=0.0001, ps=100000
    )
    # LSM = LSMModule(
    #     ps=100000,
    #     z0mav=0.0001,
    #     z0hav=0.0001,
    #     iinterp_t=1,
    #     iinterp_theta=1,
    #     dz_soil=[1.89, 0.72, 0.21, 0.07],
    #     albedoav=0.22,
    # )
    # LSM += UniformSkinTemperature(293)
    # LSM += UniformSoilTemperature([293, 293, 293, 293])
    # LSM += UniformSoilMoisture([0.3, 0.3, 0.3, 0.3])
    # LSM += LandUseModification(geometry="all", type="grs", params={})
    # sim += LSM
    #
    # IBM = IBMModule()
    # IBM += IBMModification("all", height=0, params={})
    # IBM += IBMModification(
    #     "rectangle_idx",
    #     height=20,
    #     params={"minx": 10, "maxx": 20, "miny": 10, "maxy": 20},
    # )
    # sim += IBM

    sim += nesting

    sim += EasyOutputModule(
        output_interval=30,
        enable_output=True,
    )

    if sim.nml.get("namchecksim") is None:
        sim.nml["namchecksim"] = {}
    sim.nml["namchecksim"]["tcheck"] = 60

    if sim.nml.get("thermodynamics") is None:
        sim.nml["thermodynamics"] = {}
    sim.nml["thermodynamics"]["lconstexner"] = True

    sim.sim_preprocessing_pipeline()
    exit()

    wd = os.getcwd()
    os.chdir(sim.output_path.as_posix())
    subprocess.run("./job.001", check=False)
    run_result = subprocess.run(
        ["/Users/andrevanginkel/bin/combine.sh", "run_001"], check=False
    )
    # print(run_result.stderr)
    os.chdir(wd)

    # Machine configuration (would normally come from machine_conf.yaml)
    # Create simulation instance
    case_name = "012_test_subnest"
    sim2 = dales_simulation(case_name, machine_conf)
    logger.info(f"Created simulation with case: {case_name}")

    # Add modules - this demonstrates the += operator for module registration
    sim2 += DefaultNamelistModule()
    logger.info("Added DefaultNamelistModule")

    sim2 += subgrid
    logger.info("Added GridModule")

    atmo = AtmosphereModule()
    atmo += AtmosphericProfile(
        variable=ua, shape="lin", params=dict(surf_val=0.5, ddz=0)
    )
    atmo += AtmosphericProfile(
        variable=va, shape="lin", params=dict(surf_val=3, ddz=1e-3)
    )
    atmo += AtmosphericProfile(
        variable=thetal, shape="lin", params=dict(surf_val=293.15, ddz=1e-2)
    )
    atmo += AtmosphericProfile(variable=qt, shape="lin", params=dict(surf_val=0, ddz=0))
    atmo += AtmosphericProfile(variable=wa, shape="lin", params=dict(surf_val=0, ddz=0))
    atmo += InterpolatedProfile(
        variable=tke,
        z=[0, 4000, 5000],
        points=[1, 1e-8, 1e-8],
    )

    sim2 += atmo
    sim2 += TimeModule(
        xday=1,
        xtime=12.0,
        xyear=2023,
        runtime=600,
        startyear=2023,
        startmonth=1,
        startday=1,
    )

    # sim += RadiationModule(
    #     iradiation=4,  # 0=off, 1=simple, 4=RRTMGP
    # )
    sim2 += ConstantSurfaceTemperatureModule(
        thls=293.15, z0mav=0.0001, z0hav=0.0001, ps=100000
    )

    if sim2.nml.get("namnetcdfstats") is None:
        sim2.nml["namnetcdfstats"] = {}
    sim2.nml["namnetcdfstats"]["lsync"] = True

    nesting = NestingTopology()
    # this is us
    nesting += supergrid

    nesting += subgrid
    nesting.my_idx = nesting.nestings.index(subgrid)
    sim2 += nesting

    openbc = do_openboundary(
        # tracernames=["tracer1", "tracer2"],
        time0="2023-01-01T12:00:10",
        start="2023-01-01T12:00:00",
        end="2023-01-02T12:00:00",
        e12=0.01,
        dxint=subgrid.xsize / subgrid.itot,
        dyint=subgrid.ysize / subgrid.jtot,
        tauh=20,
        taum=0,
        # linithetero=False,
    )
    openbc += Nest_in_Dales(
        # inpath=sim.output_path / "input/",
        inpath_coarse=sim.output_path / "input/",
        outpath_coarse=sim.output_path / "run_001/",
        outpath_coarse_old=sim.output_path / "run_001/",
    )
    sim2 += openbc
    sim2 += EasyOutputModule(
        output_interval=5,
        enable_output=True,
    )

    if sim2.nml.get("namchecksim") is None:
        sim2.nml["namchecksim"] = {}
    sim2.nml["namchecksim"]["tcheck"] = 60
    if sim2.nml.get("thermodynamics") is None:
        sim2.nml["thermodynamics"] = {}
    sim2.nml["thermodynamics"]["lconstexner"] = True

    sim2.sim_preprocessing_pipeline()

    wd = os.getcwd()
    os.chdir(sim2.output_path.as_posix())
    subprocess.run("./job.001", check=True)
    run_result = subprocess.run(
        ["/Users/andrevanginkel/bin/combine.sh", "run_001"], check=True
    )
    # print(run_result.stderr)
    os.chdir(wd)
    # txt = sim.save_sim_to_yaml()
    # with open("test_sim.yaml", "w") as f:
    #     f.write(txt)
    # sim2 = dales_simulation.load_sim_from_yaml(txt, machine_conf=machine_conf)
    # print(sim2)

    # yaml.dump(sim2, open("test_sim2_dump.yaml", "w"))
    # yaml.dump(sim, open("test_sim_dump.yaml", "w"))
