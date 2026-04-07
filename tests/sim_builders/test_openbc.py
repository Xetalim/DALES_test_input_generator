import subprocess
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from modular_dales import (
    AtmosphereModule,
    AtmosphericProfile,
    ConstantSurfaceTemperatureModule,
    DefaultNamelistModule,
    GridDales,
    TimeModule,
    dales_simulation,
)
from modular_dales.Configuration.output_modules import EasyOutputModule
from modular_dales.LBC import NestingTopology
from modular_dales.vars import *  # noqa: F401,F403

from tests.openbc_test_input.openbc_fixtures import (
    atmo_netcdf_file,
    surface_netcdf_file,
)


@pytest.mark.skip(reason="openBCModule not fully implemented yet")
def test_openbc(atmo_netcdf_file, surface_netcdf_file) -> None:
    """Smoke test for OpenBC module using generated NetCDF inputs.

    This is preserved from the original test file and remains skipped until
    openBCModule is fully implemented.
    """

    print("Testing openBCModule with atmospheric profile from netCDF file")

    file1 = atmo_netcdf_file(xlen=50, ylen=50, levlen=21, timesteps=2)
    subprocess.run(["ncview", str(file1)], check=True)
    print(file1)
    file2 = surface_netcdf_file(xlen=50, ylen=50, levlen=21, timesteps=2)
    subprocess.run(["ncview", str(file2)], check=True)
    print(file2)


def _build_nested_sim_with_easyoutput(
    machine_conf: dict, casename: str
) -> dales_simulation:
    """Create a small nested DALES sim with EasyOutput and cross-sections.

    This mirrors the structure of the example in oop_test, but is reduced in
    size so it is suitable for use in unit tests. The simulation consists of a
    coarse "supergrid" with a smaller "subgrid" nested inside via
    NestingTopology. EasyOutputModule enables both fielddumps and
    cross-sections so that DALES will write ``fielddump.nc`` and
    cross-section NetCDF files (e.g. ``crossxz.nc`` or ``crossyz.nc``) after
    running ``job.001``.
    """

    sim = dales_simulation(casename, machine_conf)

    # Basic namelist + coarse grid
    sim += DefaultNamelistModule()

    supergrid = GridDales(
        itot=16,
        jtot=16,
        kmax=40,
        xsize=640.0,
        ysize=640.0,
        kmax_soil=4,
        xlat=52.25,
        xlon=5.45,
        x0=0.0,
        y0=0.0,
        alpha=1.0,
        dz0=10.0,
    )

    # Subgrid fully contained inside the supergrid domain
    subgrid = GridDales(
        itot=64,
        jtot=64,
        kmax=30,
        xsize=320.0,
        ysize=320.0,
        kmax_soil=4,
        xlat=52.25,
        xlon=5.45,
        x0=160.0,
        y0=160.0,
        alpha=1.0,
        dz0=10.0,
    )

    sim += supergrid

    # Simple atmosphere with linear profiles for basic variables
    atmo = AtmosphereModule()
    atmo += AtmosphericProfile(
        variable=ua,
        shape="lin",
        params=dict(surf_val=3.0, ddz=1e-3),
    )
    atmo += AtmosphericProfile(
        variable=va,
        shape="lin",
        params=dict(surf_val=0.0, ddz=0.0),
    )
    atmo += AtmosphericProfile(
        variable=w,
        shape="lin",
        params=dict(surf_val=0.0, ddz=0.0),
    )
    atmo += AtmosphericProfile(
        variable=thetal,
        shape="lin",
        params=dict(surf_val=293.15, ddz=1e-2),
    )
    atmo += AtmosphericProfile(
        variable=qt,
        shape="lin",
        params=dict(surf_val=0.01, ddz=0.0),
    )

    sim += atmo

    # Short runtime so the produced files stay small
    sim += TimeModule(xtime=0.0, xday=1, xyear=2025, runtime=60)

    # Simple surface module
    sim += ConstantSurfaceTemperatureModule(
        thls=293.15,
        z0mav=0.0001,
        z0hav=0.0001,
        ps=100000,
        albedoav=0.22,
    )

    # Nesting topology: supergrid hosting subgrid; this also configures
    # cross-section indices in namcrosssection so DALES writes cross* files.
    nesting = NestingTopology()
    nesting += supergrid
    nesting.my_idx = nesting.nestings.index(supergrid)
    nesting += subgrid
    sim += nesting

    # Enable fielddumps and cross-sections via EasyOutputModule
    sim += EasyOutputModule(
        output_interval=1,
        enable_output=True,
    )

    # Minimal additional namelist tweaks to keep parity with other examples
    if sim.nml.get("namchecksim") is None:
        sim.nml["namchecksim"] = {}
    sim.nml["namchecksim"]["tcheck"] = 60

    if sim.nml.get("thermodynamics") is None:
        sim.nml["thermodynamics"] = {}
    sim.nml["thermodynamics"]["lconstexner"] = True

    return sim


def _assert_crosssection_matches_fielddump(
    fielddump_path: Path,
    crosssection_path: Path,
    var: str = "thl",
    fd_coord_dim: str = "xt",
    cs_coord_dim: str = "xt",
    cross_index: int = 0,
    time_index: int = 0,
    rtol: float = 1e-6,
    atol: float = 1e-8,
) -> None:
    """Compare a 2D slice between two DALES NetCDF outputs at a grid index.

    Both NetCDF files are expected to contain a variable ``var`` on a grid
    with coordinate ``fd_coord_dim`` (e.g. ``xt`` or ``yt``) in the first
    dataset and ``cs_coord_dim`` in the second. The function selects the
    coordinate value from the first dataset at ``index`` along
    ``fd_coord_dim``, locates the nearest matching coordinate in the second
    dataset along ``cs_coord_dim``, and asserts the resulting 2D slices are
    numerically equal within the provided tolerances.
    """

    with xr.open_dataset(fielddump_path) as ds_ref, xr.open_dataset(
        crosssection_path
    ) as ds_other:
        if var not in ds_ref or var not in ds_other:
            raise AssertionError(f"Variable '{var}' not found in both datasets")

        if fd_coord_dim not in ds_ref.coords:
            raise AssertionError(
                f"Coordinate '{fd_coord_dim}' not found in fielddump coords"
            )
        if cs_coord_dim not in ds_other.coords:
            raise AssertionError(
                f"Coordinate '{cs_coord_dim}' not found in crosssection coords"
            )

        coord_vals_ref = ds_ref[fd_coord_dim]
        n_ref = coord_vals_ref.sizes[fd_coord_dim]

        def _get_slices(idx_ref: int):
            """Return aligned slices for a given reference index.

            Uses the reference dataset coordinate at ``idx_ref`` to find the
            nearest matching coordinate in the other dataset and returns the
            corresponding 2D slices at ``time_index``.
            """

            coord_val = coord_vals_ref.isel({fd_coord_dim: idx_ref}).item()
            coord_vals_other = ds_other[cs_coord_dim].values
            other_index = int(np.abs(coord_vals_other - coord_val).argmin())

            field_slice_local = ds_ref[var].isel(
                {"time": time_index, fd_coord_dim: idx_ref}
            )
            cross_slice_local = ds_other[var].isel(
                {"time": time_index, cs_coord_dim: other_index}
            )
            return field_slice_local, cross_slice_local, other_index

        coord_value = ds_other.coords[cs_coord_dim].values[cross_index]
        index = list(ds_ref[fd_coord_dim].values).index(coord_value)

        # First try the requested index directly.
        base_field_slice, base_cross_slice, base_other_index = _get_slices(index)

        if np.allclose(
            base_field_slice.values, base_cross_slice.values, rtol=rtol, atol=atol
        ):
            # Slices match at the requested index: test passes.
            return
        else:
            # If that fails, try small shifts in the reference index to detect
            # a possible off-by-N error between fielddump and cross-section.
            max_shift = 50
            for shift in range(-max_shift, max_shift + 1):
                if shift == 0:
                    continue
                idx_shift = index + shift
                if idx_shift < 0 or idx_shift >= n_ref:
                    continue

                field_slice_s, cross_slice_s, other_idx_s = _get_slices(idx_shift)
                if np.allclose(
                    field_slice_s.values,
                    cross_slice_s.values,
                    rtol=rtol,
                    atol=atol,
                ):
                    # Slices only match with a shifted index: report and fail
                    raise AssertionError(
                        f"{fielddump_path.name}/{crosssection_path.name} slices do not match at index "
                        f"{index} along '{fd_coord_dim}', but do match when the "
                        f"fielddump index is shifted by {shift} (to {idx_shift}) "
                        f"and crosssection index {other_idx_s} is used along "
                        f"'{cs_coord_dim}'. This suggests an index offset of "
                        f"{shift} between fielddump and crosssection coordinates."
                    )
                else:
                    print(
                        shift,
                        np.mean(
                            np.abs(field_slice_s.values - cross_slice_s.values) ** 2
                        ),
                    )
                    # This shift did not help; try the next one.
                    continue

            # No shift produced a match: re-raise the original assertion.
            np.testing.assert_allclose(
                base_field_slice.values,
                base_cross_slice.values,
                rtol=rtol,
                atol=atol,
            )


@pytest.fixture(
    params=[
        pytest.param(1, id="1_cores"),
        pytest.param(2, id="2_cores"),
        pytest.param(4, id="4_cores"),
        pytest.param(8, id="8_cores"),
    ]
)
def core_changer(request):
    return request.param


def test_crosssection_matches_fielddump_at_grid_index(
    machine_conf, core_changer
) -> None:
    """Ensure cross-section values equal fielddump values at same grid index.

    This is an end-to-end test: it builds a small nested simulation with
    EasyOutputModule, runs ``job.001`` and ``combine.sh``, then compares a
    cross-section NetCDF file against ``fielddump.nc`` at a coordinate index
    taken from the cross-section grid.
    """

    conf = machine_conf("openbc_crosssection_fielddump")
    conf["job_conf"]["numcores"] = core_changer
    sim = _build_nested_sim_with_easyoutput(
        conf, casename="openbc_crosssection_fielddump"
    )
    sim.sim_preprocessing_pipeline()

    outdir = sim.output_path

    # Run the DALES job and combine script to generate diagnostic files.
    subprocess.run(["./job.001"], check=True, cwd=outdir.as_posix())
    subprocess.run(["combine.sh", "run_001"], check=True, cwd=outdir.as_posix())

    fielddump_file = outdir / "fielddump.nc"
    if not fielddump_file.is_file():
        pytest.skip("fielddump.nc not found; DALES run did not produce fielddump")

    # Prefer a yz cross-section at fixed x if available, fall back to others.
    candidates = [
        outdir / "crossyz.nc",
        outdir / "crossxz.nc",
        outdir / "crossxy.nc",
    ]
    var = "thl"
    for crosssection_file in candidates:
        # Choose which coordinate dimension is held constant in the cross-section.
        if crosssection_file.name == "crossyz.nc":
            if var == "u":
                fd_coord_dim = "xm"
                cs_coord_dim = "xm"
            else:
                fd_coord_dim = "xt"
                cs_coord_dim = "xt"
        elif crosssection_file.name == "crossxz.nc":
            if var == "v":
                fd_coord_dim = "ym"
                cs_coord_dim = "ym"
            else:
                fd_coord_dim = "yt"
                cs_coord_dim = "yt"
        else:  # crossxy.nc
            fd_coord_dim = "zt"
            cs_coord_dim = "zt"

        # Use the cross-section grid as reference and compare against fielddump.
        _assert_crosssection_matches_fielddump(
            fielddump_file,
            crosssection_file,
            var=var,
            fd_coord_dim=fd_coord_dim,
            cs_coord_dim=cs_coord_dim,
            cross_index=0,
            time_index=5,
        )
