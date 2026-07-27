from __future__ import annotations

import yaml

import numpy as np

from modular_dales import (
    AtmosphereModule,
    AtmosphericProfile,
    DefaultNamelistModule,
    EasyOutputModule,
    GridDales,
    InterpolatedProfile,
    RadiationModule,
    TimeModule,
    dales_simulation,
)
from modular_dales.Configuration.output_modules import CheckSimulationModule
from modular_dales.Emission.emission import (
    EmissionModule,
    EmissionPointSource,
    EmissionTracer,
)
from modular_dales.Surface.LSM.LSM import LSMModule, LandUseModification
from modular_dales.Surface.LSM.modular_temps_moisture import (
    UniformSkinTemperature,
    UniformSoilMoisture,
    UniformSoilTemperature,
)
from modular_dales.Surface import (
    ConstantSurfaceTemperatureModule,
)
from modular_dales.modular.simulation_module import set_nml_section
from modular_dales.modular.time_dependent import TimedependentModule
from modular_dales.modular.time_dependent_scalars import TimeDependentScalar
from modular_dales.vars import (
    VariableDefinition,
    get_var_by_name,
    qt,
    register_var,
    tke,
    thetal,
    ua,
    ug,
    va,
    vg,
    wa,
)

if __name__ == "__main__":
    with open("machine_conf.yaml", "r", encoding="utf-8") as f:
        machine_conf = yaml.safe_load(f)

    sim = dales_simulation("004_daan", machine_conf)
    sim += DefaultNamelistModule()

    sim += GridDales(
        itot=64 // 2,
        jtot=64 // 2,
        kmax=96,
        xsize=3200.0 / 2,
        ysize=3200.0 / 2,
        kmax_soil=4,
        xlat=52.25,
        xlon=5.45,
        x0=0.0,
        y0=0.0,
        alpha=1.02,
        dz0=4.0,
    )

    atmo = AtmosphereModule()
    atmo += AtmosphericProfile(
        variable=ua,
        shape="lin",
        params=dict(surf_val=1.0, ddz=0.0),
    )
    atmo += AtmosphericProfile(
        variable=va,
        shape="lin",
        params=dict(surf_val=0.0, ddz=0.0),
    )
    atmo += AtmosphericProfile(
        variable=ug,
        shape="lin",
        params=dict(surf_val=1.0, ddz=0.0),
    )
    atmo += AtmosphericProfile(
        variable=vg,
        shape="lin",
        params=dict(surf_val=0.0, ddz=0.0),
    )
    atmo += AtmosphericProfile(
        variable=wa,
        shape="lin",
        params=dict(surf_val=0.0, ddz=0.0),
    )

    # Neutral boundary layer up to 1 km, with stable stratification above.
    def func(z, lapse_rate, ml_val, offset_val=1.25, z_ml=500):
        u = np.zeros(z.shape, dtype=float)
        u[z <= z_ml] = ml_val
        u[z > z_ml] = ml_val + offset_val + (z[z > z_ml] - z_ml) * lapse_rate
        return u

    z = np.linspace(0, 5000, 1000, dtype=float)
    atmo += InterpolatedProfile(
        variable=thetal,
        z=z.tolist(),
        points=func(z, lapse_rate=0.006, ml_val=295.6, offset_val=5, z_ml=500).tolist(),
    )
    atmo += AtmosphericProfile(
        variable=qt,
        shape="lin",
        params=dict(surf_val=0.000, ddz=0.0),
    )
    atmo += InterpolatedProfile(
        variable=tke,
        z=[0.0, 500.0, 2500.0, 5000.0],
        points=[0.2, 0.1, 1.0e-3, 1.0e-4],
    )
    sim += atmo

    times = [0, 1800, 3600 * 6, 3600 * 6 + 1]
    values = [301.0, 301, 304.5, 305.5]
    sim += TimedependentModule(timesteps=times)
    sim += ConstantSurfaceTemperatureModule(
        thls=TimeDependentScalar(times, values), z0mav=0.1, z0hav=0.01, ps=101325.0
    )

    sim += EasyOutputModule(output_interval=120, enable_output=True)
    sim += TimeModule(
        xtime=12.0,
        xday=180,
        xyear=2025,
        runtime=3600 * 6,
        startyear=2025,
        startmonth=6,
        startday=29,
    )
    sim += CheckSimulationModule(check_interval=120)
    set_nml_section(
        sim.nml, sim.nml_docs, "user_defined", "namnetcdfstats", "lsync", True
    )
    set_nml_section(
        sim.nml, sim.nml_docs, "user_defined", "thermodynamics", "lmoist", False
    )
    # set_nml_section(
    #     sim.nml, sim.nml_docs, "user_defined", "namsubgrid", "lanisotropic", True
    # )
    sim.sim_preprocessing_pipeline()
    print(f"Case generated at: {sim.output_path}")
