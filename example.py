"""Compact KNMI -> DALES open-boundary case setup."""

import logging
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
from modular_dales.Geometry import GridDales
from modular_dales.IBM.IBM import FromAHN, IBMModule
from modular_dales.LBC import Nest_in_KNMI, NestingTopology, do_openboundary
from modular_dales.logging_wrapper import setup_logging
from modular_dales.modular import dales_simulation
from modular_dales.Surface import ConstantSurfaceTemperatureModule
from modular_dales.modular.simulation_module import set_nml_section
from modular_dales.vars import qt, thetal, tke, ua, ug, va, vg, wa

setup_logging("logging.yaml")
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

if __name__ == "__main__":
    with open("machine_conf.yaml", "r", encoding="utf-8") as f:
        machine_conf = yaml.safe_load(f)

    logger.info("=" * 70)
    logger.info("DALES nested in KNMI HARMONIE")
    logger.info("=" * 70)

    case_name = "knmi_nest_test"
    grid = GridDales(
        itot=128,
        jtot=128,
        kmax=80,
        xsize=512,  # m
        ysize=512,  # m
        kmax_soil=4,
        xlat=51.9990,
        xlon=4.3731,
        x0=85225,
        y0=445930,
        alpha=1.01,
        dz0=2,
        proj4=(
            "+proj=sterea +lat_0=52.1561605555556 +lon_0=5.38763888888889 "
            "+k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel "
            "+towgs84=565.4171,50.3319,465.5524,-0.398957,0.343988,-1.87740,4.0725 "
            "+units=m +no_defs"
        ),
    )

    sim = dales_simulation(case_name, machine_conf)
    logger.info("Created simulation: %s", case_name)
    sim += DefaultNamelistModule()
    sim += grid

    atmo = AtmosphereModule()
    for variable in (ua, ug, va, vg, wa):
        atmo += AtmosphericProfile(
            variable=variable, shape="lin", params={"surf_val": 3, "ddz": 0}
        )
    atmo += AtmosphericProfile(
        variable=thetal, shape="lin", params={"surf_val": 290.0, "ddz": 3.5e-3}
    )
    atmo += AtmosphericProfile(
        variable=qt, shape="lin", params={"surf_val": 0.008, "ddz": -2.5e-6}
    )
    atmo += InterpolatedProfile(
        variable=tke,
        z=[0, 4000, 5000],
        points=[1, 1e-8, 1e-8],
    )
    sim += atmo

    sim += TimeModule(
        xday=26,
        xtime=18.0,
        xyear=2026,
        runtime=21600,
        startyear=2026,
        startmonth=3,
        startday=26,
    )
    sim += ConstantSurfaceTemperatureModule(
        thls=290.0, z0mav=0.01, z0hav=0.01, ps=101325.0
    )

    # openbc = do_openboundary(
    #     time0="2026-03-26T18:00:00",
    #     start="2026-03-26T18:00:00",
    #     end="2026-03-27T00:00:00",
    #     e12=0.01,
    #     tchunk=1,
    #     dxint=2.0 * 4,
    #     dyint=2.0 * 4,
    #     dxturb=2.0,
    #     dyturb=2.0,
    #     linithetero=True,
    # )
    # openbc += Nest_in_KNMI(
    #     ml_glob="/Users/andrevanginkel/Documents/30_Data/32_campus_case/32.01_harm_2026032618/HARM43_V1_P5_2026032618/HA43_N55ML.nc",
    #     sfc_glob="/Users/andrevanginkel/Documents/30_Data/32_campus_case/32.01_harm_2026032618/HARM43_V1_P3_2026032618/HA43_N55.nc",
    #     w_from_continuity=True,
    # )
    # sim += openbc

    IBM = IBMModule()
    IBM += FromAHN()
    sim += IBM

    sim += EasyOutputModule(output_interval=1, enable_output=True)
    sim += SamplingModule(output_interval=1, enable_output=True)

    for section, key, value in [
        ("namnetcdfstats", "lsync", True),
        ("physics", "lcoriol", False),
        ("namchecksim", "tcheck", 1),
        ("run", "nprocx", 0),
        ("run", "nprocy", 0),
    ]:
        set_nml_section(sim.nml, sim.nml_docs, "user_defined", section, key, value)

    logger.info("Running preprocessing pipeline...")
    sim.sim_preprocessing_pipeline()
    logger.info("Done! Case written to: %s", sim.output_path)
