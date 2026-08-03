"""Build a 64x64 periodic-precursor + recycled-openboundary LS2D setup.

Workflow:
1. Create and run a periodic 64x64 precursor with LS2D nudging.
2. Force precursor cross-sections (edge/middle in xz/yz and one xy top layer).
3. Create a 64x64 open-boundary child nested in LS2D atmosphere profiles.
4. Add periodic turbulence perturbations from precursor cross-sections.
"""

from __future__ import annotations

import datetime
import logging
import os
from pathlib import Path
import subprocess

import dask
import xarray as xr
import yaml

from modular_dales import (
    DefaultNamelistModule,
    GridDales,
    LSMModule,
    Nest_in_AtmosphereProfiles,
    Periodic_Dales_Turbulence_Perturbations,
    PeriodicPrecursorCrossSections,
    UniformSkinTemperature,
    UniformSoilMoisture,
    UniformSoilTemperature,
    dales_simulation,
    do_openboundary,
)
from modular_dales.Atmosphere.ls2d_atmosphere import LS2DAtmosphereModule
from modular_dales.logging_wrapper import setup_logging

from triple_nest import _attach_common_physics, make_integer_stretched_grid

setup_logging("logging.yaml")
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

_GRID64 = dict(
    itot=64,
    jtot=64,
    kmax=72,
    xsize=6400.0 * 1.5,
    ysize=6400.0 * 1.5,
    kmax_soil=4,
    xlat=52.25,
    xlon=5.45,
    x0=136173.0 - ((6400.0 * 1.5) / 2.0),
    y0=455912.0 - ((6400.0 * 1.5) / 2.0),
    proj4="EPSG:28992",
    alpha=1.02,
    dz0=20.0,
)


def _new_grid(spec: dict) -> GridDales:
    grid = GridDales(**spec)
    zt, zm, dz = make_integer_stretched_grid(grid.kmax, grid.dz0, grid.alpha)
    grid.zt = zt
    grid.zm = zm
    grid.dz = dz
    grid.zsize = float(grid.zm[-1])
    return grid


def _run_job_direct(case_dir: Path, stage_name: str) -> None:
    logger.info("Running %s in %s", stage_name, case_dir)
    subprocess.run(["./job.001"], cwd=case_dir, check=True)


def _build_periodic_precursor_64(machine_conf: dict) -> dales_simulation:
    sim = dales_simulation("periodic_recycle_precursor_64", machine_conf)
    sim += DefaultNamelistModule()

    grid = _new_grid(_GRID64)
    sim += grid

    _attach_common_physics(
        sim,
        grid,
        center_vertical_crosssections=False,
        use_ls2d=True,
        is_nested=False,
    )

    # Force old-style DALES cross-sections at edges/middle + selected top layer.
    sim += PeriodicPrecursorCrossSections(top_layer_index=-1)

    return sim


def _copy_surface_state_from_parent(
    parent_sim: dales_simulation,
    child_sim: dales_simulation,
) -> None:
    sim_lsm = child_sim.retrieve_module(LSMModule)

    with xr.open_dataset(parent_sim.output_path / "run_001" / "surfcross.001.nc") as ds:
        sim_lsm += UniformSkinTemperature(ds["tskin"].isel(time=0).mean().item())

    with xr.open_dataset(parent_sim.output_path / "input" / "lsm.inp_001.nc") as ds:
        sim_lsm += UniformSoilTemperature(
            ds["t_soil"].mean(dim=("x", "y")).values.tolist()
        )
        sim_lsm += UniformSoilMoisture(
            ds["theta_soil"].mean(dim=("x", "y")).values.tolist()
        )


def _build_openboundary_child_64(
    machine_conf: dict,
    precursor_sim: dales_simulation,
) -> dales_simulation:
    sim = dales_simulation("periodic_recycle_openbc_64", machine_conf)
    sim += DefaultNamelistModule()

    grid = _new_grid(_GRID64)
    sim += grid

    _attach_common_physics(
        sim,
        grid,
        center_vertical_crosssections=False,
        use_ls2d=True,
        is_nested=True,
    )

    ls2d_atmo = sim.retrieve_module(LS2DAtmosphereModule)

    openbc = do_openboundary(
        time0="2026-05-01T12:00:00",
        start="2026-05-01T12:00:00",
        end="2026-05-01T13:00:00",
        e12=0.01,
        dxint=grid.xsize / grid.itot,
        dyint=grid.ysize / grid.jtot,
        tchunk=100,
    )
    openbc += Nest_in_AtmosphereProfiles(atmosphere_module=ls2d_atmo)
    openbc += Periodic_Dales_Turbulence_Perturbations(
        periodic_outpath=(precursor_sim.output_path / "run_001").as_posix(),
        filter_scale_m=2000.0,
        top_layer_index=-1,
        perturbation_variables=["u", "v", "w", "thl", "qt"],
        tau=100.0,
    )
    sim += openbc

    _copy_surface_state_from_parent(precursor_sim, sim)
    return sim


def build_periodic_recycle_case(
    machine_conf: dict,
) -> tuple[dales_simulation, dales_simulation]:
    precursor = _build_periodic_precursor_64(machine_conf)
    precursor.sim_preprocessing_pipeline()
    # _run_job_direct(precursor.output_path, "periodic_precursor_64")

    child = _build_openboundary_child_64(machine_conf, precursor)
    child.sim_preprocessing_pipeline()
    return precursor, child


if __name__ == "__main__":
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    # dask.config.set(scheduler="threads")
    # dask.config.set(scheduler="single-threaded")
    import webbrowser

    from dask.distributed import LocalCluster, Client

    n_cores = os.cpu_count() or 2

    cluster = LocalCluster(
        n_workers=0,
        threads_per_worker=1,
        memory_limit="8GB",
        dashboard_address=":8787",
    )

    cluster.adapt(
        minimum=1,
        maximum=max(1, n_cores - 1),
        wait_count=3,
        interval="1s",
    )

    client = Client(cluster)

    webbrowser.open(cluster.dashboard_link)

    print(f"Dask dashboard: {cluster.dashboard_link}")
    machine_conf_path = Path(__file__).resolve().parent / "machine_conf.yaml"
    with machine_conf_path.open("r", encoding="utf-8") as handle:
        machine_conf = yaml.safe_load(handle)

    logger.info("Building periodic recycle setup")
    precursor_sim, child_sim = build_periodic_recycle_case(machine_conf)

    _run_job_direct(child_sim.output_path, "periodic_recycle_openbc_64")

    logger.info("Precursor output path: %s", precursor_sim.output_path)
    logger.info("Open-boundary child output path: %s", child_sim.output_path)
