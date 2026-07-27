from __future__ import annotations

import yaml

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
from modular_dales.modular.simulation_module import set_nml_section
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
        itot=64,
        jtot=64,
        kmax=128,
        xsize=3200.0,
        ysize=3200.0,
        kmax_soil=4,
        xlat=52.25,
        xlon=5.45,
        x0=0.0,
        y0=0.0,
        alpha=1.02,
        dz0=10.0,
    )

    if "co2" not in get_var_by_name():
        register_var(
            VariableDefinition(
                name="co2",
                long_name="Carbon Dioxide (CO2)",
                unit="ppm",
                can_be_time_dependent=True,
            )
        )
    co2 = get_var_by_name()["co2"]

    atmo = AtmosphereModule()
    atmo += AtmosphericProfile(
        variable=ua,
        shape="lin",
        params=dict(surf_val=3.0, ddz=0.0),
    )
    atmo += AtmosphericProfile(
        variable=va,
        shape="lin",
        params=dict(surf_val=3.0, ddz=0.0),
    )
    atmo += AtmosphericProfile(
        variable=ug,
        shape="lin",
        params=dict(surf_val=3.0, ddz=0.0),
    )
    atmo += AtmosphericProfile(
        variable=vg,
        shape="lin",
        params=dict(surf_val=3.0, ddz=0.0),
    )
    atmo += AtmosphericProfile(
        variable=wa,
        shape="lin",
        params=dict(surf_val=0.0, ddz=0.0),
    )

    # Neutral boundary layer up to 1 km, with stable stratification above.
    atmo += InterpolatedProfile(
        variable=thetal,
        z=[0.0, 1000.0, 3000.0, 5000.0],
        points=[289.0, 289.0, 297.0, 305.0],
    )
    atmo += AtmosphericProfile(
        variable=qt,
        shape="lin",
        params=dict(surf_val=0.007, ddz=0.0),
    )
    atmo += InterpolatedProfile(
        variable=tke,
        z=[0.0, 500.0, 2500.0, 5000.0],
        points=[0.2, 0.1, 1.0e-3, 1.0e-4],
    )
    atmo += AtmosphericProfile(
        variable=co2,
        shape="lin",
        params=dict(surf_val=0, ddz=0.0),
    )
    sim += atmo

    lsm = LSMModule(
        ps=101325.0,
        z0mav=1.0e-4,
        z0hav=1.0e-4,
        iinterp_t=1,
        iinterp_theta=1,
        dz_soil=[1.89, 0.72, 0.21, 0.07],
        albedoav=0.08,
    )
    lsm += UniformSkinTemperature(289.0)
    lsm += UniformSoilTemperature([289.0, 289.0, 289.0, 289.0])
    lsm += UniformSoilMoisture([0.98, 0.98, 0.98, 0.98])
    lsm += LandUseModification(geometry="all", type="aqu", params={})
    sim += lsm

    emis = EmissionModule()
    emis += EmissionTracer(
        name="co2",
        long_name="Carbon Dioxide (CO2)",
        unit="ppm",
        molar_mass=44.009,
        lemis=True,
    )
    emis += EmissionPointSource(
        tracer_name="co2",
        x_idx=4,
        y_idx=32,
        height=20.0,
        temperature=289.0,
        volume=1.0,
        emission=10.0,
        stack_exit_area=1.0,
    )
    sim += emis

    sim += RadiationModule(iradiation=4)
    sim += EasyOutputModule(output_interval=60, enable_output=True)
    sim += TimeModule(
        xtime=12.0,
        xday=180,
        xyear=2025,
        runtime=36000,
        startyear=2025,
        startmonth=6,
        startday=29,
    )
    sim += CheckSimulationModule(check_interval=60)
    set_nml_section(
        sim.nml, sim.nml_docs, "user_defined", "namsubgrid", "lanisotrop", True
    )  # anisotropic SGS speeds up the simulation a lot, and we don't really care about the small inconsistencies, so we set to True
    sim.sim_preprocessing_pipeline()
    print(f"Case generated at: {sim.output_path}")
