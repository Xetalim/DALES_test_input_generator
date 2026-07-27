"""Example: single DALES case nested in KNMI operational HARMONIE data.

Usage:
    # Set the paths to your pre-converted KNMI NetCDF files, then run:
    python openbc_knmi.py

    # The script will create a DALES case directory with initial and
    # boundary conditions derived from the KNMI HARMONIE output.
"""

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
from modular_dales.Radiation.radiation import RadiationModule
from modular_dales.Surface.LSM.SLuRB.slurb import SLURBModule
from modular_dales.logging_wrapper import setup_logging
from modular_dales.modular import dales_simulation
from modular_dales.Surface import LSM, ConstantSurfaceTemperatureModule
from modular_dales.modular.simulation_module import set_nml_section
from modular_dales.vars import *

setup_logging("logging.yaml")
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# ============================================================================
# USER CONFIGURATION — edit these variables to match your data
# ============================================================================

# Glob patterns pointing to pre-converted KNMI HARMONIE NetCDF files.
# N55ML: 3D fields on hybrid levels (u, v, t, q).
# N55 or N20: surface fields (pres, t(2m), r(2m), etc.).
KNMI_ML_GLOB = "/Users/andrevanginkel/Documents/60_fun/various_knmi_apis/data/HARM43_V1_P5_2026062721/out.nc"
KNMI_SFC_GLOB = "/Users/andrevanginkel/Documents/60_fun/various_knmi_apis/data/HARM43_V1_P3_2026062721/out.nc"

# Domain centre — TU Delft campus
# EPSG:28992 (RD New): x ≈ 83000, y ≈ 447000
CENTRE_LAT = 51.9990
CENTRE_LON = 4.3731

# Domain origin in EPSG:28992 (meters).  The domain is placed so that
# the TU Delft campus is roughly centred.
RD_X0 = 85225  # - 6400.0  # 78960 m — southwest corner x
RD_Y0 = 345930  # - 6400.0  # 439710 m — southwest corner y

# EPSG:28992 proj4 string
PROJ4_RD = (
    "+proj=sterea +lat_0=52.1561605555556 +lon_0=5.38763888888889 "
    "+k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel "
    "+towgs84=565.4171,50.3319,465.5524,-0.398957,0.343988,-1.87740,4.0725 "
    "+units=m +no_defs"
)

# Domain size in grid points and meters
ITOT, JTOT = 512, 512
DX, DY = 200.0, 200.0  # grid spacing in meters
XSIZE = ITOT * DX
YSIZE = JTOT * DY

# Vertical levels
KMAX = 128
DZ0 = 20  # lowest level spacing (stretched grid)

# Time settings
START = "2026-06-28T02:00:00"
END = "2026-06-28T03:00:00"
RUNTIME = 21600  # seconds (6 hours)
START_YEAR, START_MONTH, START_DAY = 2026, 6, 28
START_HOUR = 2.0

# Surface parameters (simple constant-temperature surface)
THLS = 294.0  # surface liquid pot. temperature (K)
PS = 101325.0  # surface pressure (Pa)

# Derive w from horizontal wind divergence? (False → w=0)
W_FROM_CONTINUITY = True

# ============================================================================

if __name__ == "__main__":
    with open("machine_conf.yaml", "r") as f:
        machine_conf = yaml.safe_load(f)

    logger.info("=" * 70)
    logger.info("DALES nested in KNMI HARMONIE")
    logger.info("=" * 70)

    case_name = "knmi_nest_test"

    # --- Grid ---------------------------------------------------------------
    # x0/y0 in EPSG:28992 (RD New) meters, centred on TU Delft campus.
    grid = GridDales(
        itot=ITOT,
        jtot=JTOT,
        kmax=KMAX,
        xsize=XSIZE,
        ysize=YSIZE,
        kmax_soil=4,
        xlat=CENTRE_LAT,
        xlon=CENTRE_LON,
        x0=RD_X0,
        y0=RD_Y0,
        alpha=1.025,
        dz0=DZ0,
        proj4=PROJ4_RD,
    )

    # --- Nesting topology (single domain, no children) ----------------------
    nesting = NestingTopology()
    nesting += grid
    nesting.my_idx = nesting.nestings.index(grid)

    # --- Simulation ---------------------------------------------------------
    sim = dales_simulation(case_name, machine_conf)
    logger.info("Created simulation: %s", case_name)

    sim += DefaultNamelistModule()
    sim += grid

    # --- Atmosphere (placeholder profiles; overridden by open BC) -----------
    atmo = AtmosphereModule()
    atmo += AtmosphericProfile(variable=ua, shape="lin", params=dict(surf_val=0, ddz=0))
    atmo += AtmosphericProfile(variable=ug, shape="lin", params=dict(surf_val=0, ddz=0))
    atmo += AtmosphericProfile(variable=va, shape="lin", params=dict(surf_val=0, ddz=0))
    atmo += AtmosphericProfile(variable=vg, shape="lin", params=dict(surf_val=0, ddz=0))
    atmo += AtmosphericProfile(
        variable=thetal, shape="lin", params=dict(surf_val=THLS, ddz=3.5e-3)
    )
    atmo += AtmosphericProfile(
        variable=qt, shape="lin", params=dict(surf_val=0.008, ddz=-2.5e-6)
    )
    atmo += AtmosphericProfile(variable=wa, shape="lin", params=dict(surf_val=0, ddz=0))
    atmo += InterpolatedProfile(
        variable=tke,
        z=[0, 4000, 5000],
        points=[1, 1e-8, 1e-8],
    )
    sim += atmo

    # --- Time ---------------------------------------------------------------
    sim += TimeModule(
        xday=START_DAY,
        xtime=START_HOUR,
        xyear=START_YEAR,
        runtime=RUNTIME,
        startyear=START_YEAR,
        startmonth=START_MONTH,
        startday=START_DAY,
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
    # --- Nesting topology ---------------------------------------------------
    sim += nesting

    # --- Open boundary conditions from KNMI HARMONIE output -----------------
    openbc = do_openboundary(
        time0=START,
        start=START,
        end=END,
        e12=0.01,
        tchunk=1,
        dxint=DX * 4,
        dyint=DY * 4,
        # tauh=600,
        # taum=600,
        dxturb=DX,
        dyturb=DY,
        linithetero=True,
    )
    openbc += Nest_in_KNMI(
        ml_glob=KNMI_ML_GLOB,
        sfc_glob=KNMI_SFC_GLOB,
        w_from_continuity=W_FROM_CONTINUITY,
    )
    sim += openbc

    # IBM = IBMModule()
    # IBM += FromAHN()
    # sim += IBM

    # --- Output -------------------------------------------------------------
    sim += EasyOutputModule(output_interval=60, enable_output=True)
    # sim += SamplingModule(output_interval=1, enable_output=True)

    # --- Namelist tweaks ----------------------------------------------------

    set_nml_section(
        sim.nml, sim.nml_docs, "user_defined", "namnetcdfstats", "lsync", True
    )
    set_nml_section(sim.nml, sim.nml_docs, "user_defined", "physics", "lcoriol", False)
    set_nml_section(sim.nml, sim.nml_docs, "user_defined", "namchecksim", "tcheck", 60)
    set_nml_section(sim.nml, sim.nml_docs, "user_defined", "run", "nprocx", 0)
    set_nml_section(sim.nml, sim.nml_docs, "user_defined", "run", "nprocy", 0)
    # set_nml_section(
    #     sim.nml, sim.nml_docs, "user_defined", "NAMGREYZONE", "lgreyzone", True
    # )
    # set_nml_section(sim.nml, sim.nml_docs, "user_defined", "NAMGREYZONE", "ledmf", True)
    # set_nml_section(
    #     sim.nml, sim.nml_docs, "user_defined", "NAMGREYZONE", "lscaleaware", True
    # )

    # --- Run the preprocessing pipeline -------------------------------------
    logger.info("Running preprocessing pipeline...")
    sim.sim_preprocessing_pipeline()
    logger.info("Done! Case written to: %s", sim.output_path)
