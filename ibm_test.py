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
import pytest
import yaml

from modular_dales import dales_simulation
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
from modular_dales.IBM import FromAHN, IBMModule
from modular_dales.logging_wrapper import setup_logging
from modular_dales.Surface import ConstantSurfaceTemperatureModule
from modular_dales.vars import ua, va, thetal, qt, wa, tke

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
    pytest.skip("unsupported case")
    with open("machine_conf.yaml", "r") as file:
        machine_conf = yaml.safe_load(file)
    logger.info("=" * 70)
    logger.info("EXAMPLE 1: Basic Simulation Setup")
    logger.info("=" * 70)

    # Create minimal configuration with just case name and output directory
    case_name = "011_test"
    output_directory = None

    # Simple but valid domain configuration
    domain_info = GridDales(
        itot=128,
        jtot=128,
        kmax=40,
        xsize=512.0,
        ysize=512.0,
        kmax_soil=4,
        xlat=52.25,
        xlon=5.45,
        x0=137067.966091,
        y0=454451.510764,
        alpha=1.0,
        dz0=5.0,
        proj4="+proj=sterea +lat_0=52.1561605555556 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.4171,50.3319,465.5524,1.9342,-1.6677,9.1019,4.0725 +units=m +no_defs +type=crs",
    )

    # Machine configuration (would normally come from machine_conf.yaml)
    # Create simulation instance
    sim = dales_simulation(case_name, machine_conf)
    logger.info("Created simulation with case: %s", case_name)

    # Add modules - this demonstrates the += operator for module registration
    sim += DefaultNamelistModule()
    logger.info("Added DefaultNamelistModule")

    sim += domain_info
    logger.info("Added GridModule")

    atmo = AtmosphereModule()
    atmo += AtmosphericProfile(variable=ua, shape="lin", params=dict(surf_val=3, ddz=0))
    atmo += AtmosphericProfile(
        variable=va, shape="lin", params=dict(surf_val=0, ddz=1e-3)
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
        runtime=3600,
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
    IBM = IBMModule()
    IBM += FromAHN()
    # IBM += IBMModification("all", height=0, params={})
    # IBM += IBMModification(
    #     "rectangle_idx",
    #     height=20,
    #     params={"minx": 10, "maxx": 20, "miny": 10, "maxy": 20},
    # )
    sim += IBM

    # sim += nesting

    sim += EasyOutputModule(
        output_interval=5,
        enable_output=True,
    )

    if sim.nml.get("namchecksim") is None:
        sim.nml["namchecksim"] = {}
    sim.nml["namchecksim"]["tcheck"] = 60

    # if sim.nml.get("thermodynamics") is None:
    #     sim.nml["thermodynamics"] = {}
    # sim.nml["thermodynamics"]["lconstexner"] = True

    sim.init_output_folder()
    sim.setup_module_links()
    sim.do_config()
    sim.check_settings()
    sim.prepare_all_calculations()
    sim.apply_job_configuration()
    sim.write_module_files()
    sim.write_simulation_files()
    # get_process_ahn(sim.grid).to_netcdf(sim.output_path / "input" / "ibm.inp.001.nc")

    # exit()

    # wd = os.getcwd()
    # os.chdir(sim.output_path.as_posix())
    # subprocess.run("./job.001")
    # run_result = subprocess.run(
    #     ["/Users/andrevanginkel/bin/combine.sh", "run_001"], check=False
    # )
    # # print(run_result.stderr)
    # os.chdir(wd)
