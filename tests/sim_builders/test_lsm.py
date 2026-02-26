from __future__ import annotations

import numpy as np

from modular_dales import dales_simulation
from modular_dales.Atmosphere import (
    AtmosphereModule,
    AtmosphericProfile,
    InterpolatedProfile,
)
from modular_dales.Configuration.defaultnamelist import DefaultNamelistModule
from modular_dales.Configuration.run_and_time import TimeModule
from modular_dales.Geometry.GridDales import GridDales
from modular_dales.Radiation.radiation import RadiationModule
from modular_dales.Surface.LSM.LSM import LSMModule, LandUseModification
from modular_dales.Surface.LSM.modular_temps_moisture import (
    UniformSkinTemperature,
    UniformSoilMoisture,
    UniformSoilTemperature,
)
from modular_dales.vars import *  # noqa: F401,F403


def simple_LSM_case(machine_conf: dict) -> dales_simulation:
    """Construct a simple simulation including an LSM configuration."""
    sim = dales_simulation("LSM_case", machine_conf)

    sim += DefaultNamelistModule()

    domain_info = GridDales(
        itot=16,
        jtot=16,
        kmax=20,
        xsize=160.0,
        ysize=160.0,
        kmax_soil=4,
        xlat=52.25,
        xlon=5.45,
        x0=0,
        y0=0,
        alpha=1.0,
        dz0=5.0,
    )
    sim += domain_info

    atmo = AtmosphereModule()
    atmo += AtmosphericProfile(
        variable=ua, shape="lin", params=dict(surf_val=3, ddz=1e-3)
    )
    atmo += AtmosphericProfile(variable=va, shape="lin", params=dict(surf_val=0, ddz=0))
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
    lsm = LSMModule(
        ps=100000,
        z0mav=0.0001,
        z0hav=0.0001,
        iinterp_t=1,
        iinterp_theta=1,
        dz_soil=[1.89, 0.72, 0.21, 0.07],
        albedoav=0.22,
    )
    lsm += UniformSkinTemperature(293)
    lsm += UniformSoilTemperature([293, 293, 293, 293])
    lsm += UniformSoilMoisture([0.3, 0.3, 0.3, 0.3])
    lsm += LandUseModification(geometry="all", type="grs", params={})
    sim += lsm

    sim += RadiationModule(
        iradiation=4,
    )
    time_module = TimeModule(xtime=0.0, xday=1, xyear=2025, runtime=60)
    sim += time_module
    return sim


def assert_lsm_files_written(machine_conf) -> None:
    """Test that LSM files are written correctly for simple_LSM_case."""

    sim = simple_LSM_case(machine_conf("lsm_files"))
    sim.sim_preprocessing_pipeline()

    lsm_files = list(sim.output_path.glob("input/lsm.inp_*.nc"))
    assert len(lsm_files) == 1, f"Expected 1 LSM file, found {len(lsm_files)}"
    for lsm_file in lsm_files:
        assert lsm_file.is_file(), f"Expected LSM file {lsm_file} to exist"
        assert (
            lsm_file.stat().st_size > 0
        ), f"Expected LSM file {lsm_file} to be non-empty"
