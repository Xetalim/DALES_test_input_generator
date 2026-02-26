from __future__ import annotations

import pytest

from .generic_testers import (
    assert_roundtrip_simulation_outputs_identical,
    run_simulation_and_check_job,
    assert_files_written,
)
from .sim_builders.test_basic import _build_basic_sim
from .sim_builders.test_emission import emission_case
from .sim_builders.test_ibm import (
    IBM_case,
    assert_ibm_files_written,
)
from .sim_builders.test_lsm_harmonie_soil import (
    LSM_with_HARM_soil,
    assert_lsm_HARM_files_written,
)
from .sim_builders.test_lsm import simple_LSM_case, assert_lsm_files_written
from .sim_builders.test_sim2 import sim_2_case
from .sim_builders.test_slurb import (
    assert_slurb_files_written,
    simple_SLURB_case,
)
from .sim_builders.test_slurb_lcz import (
    SLURB_LCZ_case,
    assert_slurb_lcz_files_written,
)
from .sim_builders.test_timedep_atmosphere import (
    timedep_atmosphere_case,
    assert_timedep_atmosphere_files_written,
)
from .sim_builders.test_timedep_atmosphere_params import (
    timedep_atmosphere_with_params_case,
    assert_timedep_atmosphere_with_params_files_written,
)
from .sim_builders.test_timedep_atmosphere_params_long import (
    timedep_atmosphere_with_params_case_long,
)
from .sim_builders.test_openbc_atmosphere_profiles import (
    openbc_atmosphere_profiles_case,
    openbc_atmosphere_profiles_timedep_case,
    assert_openbc_atmosphere_profiles_files_written,
    assert_openbc_atmosphere_profiles_timedep_files_written,
)


@pytest.fixture(
    params=[
        pytest.param(_build_basic_sim, id="_build_basic_sim"),
        pytest.param(timedep_atmosphere_case, id="timedep_atmosphere_case"),
        pytest.param(
            timedep_atmosphere_with_params_case,
            id="timedep_atmosphere_with_params_case",
        ),
        pytest.param(sim_2_case, id="sim_2_case"),
        pytest.param(IBM_case, id="IBM_case"),
        pytest.param(simple_LSM_case, id="simple_LSM_case"),
        pytest.param(LSM_with_HARM_soil, id="LSM_with_HARM_soil"),
        pytest.param(SLURB_LCZ_case, id="SLURB_LCZ_case"),
        pytest.param(simple_SLURB_case, id="simple_SLURB_case"),
        pytest.param(emission_case, id="emission_case"),
        pytest.param(
            openbc_atmosphere_profiles_case,
            id="openbc_atmosphere_profiles_case",
        ),
        pytest.param(
            openbc_atmosphere_profiles_timedep_case,
            id="openbc_atmosphere_profiles_timedep_case",
        ),
    ]
)
def sim_builder(request):
    return request.param


@pytest.fixture(
    params=[
        pytest.param(_build_basic_sim, id="_build_basic_sim"),
        pytest.param(timedep_atmosphere_case, id="timedep_atmosphere_case"),
        pytest.param(
            timedep_atmosphere_with_params_case,
            id="timedep_atmosphere_with_params_case",
        ),
        pytest.param(sim_2_case, id="sim_2_case"),
        pytest.param(IBM_case, id="IBM_case"),
        pytest.param(simple_LSM_case, id="simple_LSM_case"),
        pytest.param(LSM_with_HARM_soil, id="LSM_with_HARM_soil"),
        pytest.param(SLURB_LCZ_case, id="SLURB_LCZ_case"),
        pytest.param(simple_SLURB_case, id="simple_SLURB_case"),
        pytest.param(emission_case, id="emission_case"),
        pytest.param(
            openbc_atmosphere_profiles_case,
            id="openbc_atmosphere_profiles_case",
        ),
        pytest.param(
            openbc_atmosphere_profiles_timedep_case,
            id="openbc_atmosphere_profiles_timedep_case",
        ),
    ]
)
def sim_builder_to_run(request):
    return request.param


@pytest.fixture(
    params=[
        pytest.param(
            assert_timedep_atmosphere_files_written,
            id="timedep_atmosphere_case",
        ),
        pytest.param(
            assert_timedep_atmosphere_with_params_files_written,
            id="timedep_atmosphere_with_params_case",
        ),
        pytest.param(sim_2_case, id="sim_2_case"),
        pytest.param(assert_ibm_files_written, id="IBM_case"),
        pytest.param(assert_slurb_lcz_files_written, id="SLURB_LCZ_case"),
        pytest.param(assert_slurb_files_written, id="SLURB_case"),
        pytest.param(
            assert_openbc_atmosphere_profiles_files_written,
            id="openbc_atmosphere_profiles_case",
        ),
        pytest.param(
            assert_openbc_atmosphere_profiles_timedep_files_written,
            id="openbc_atmosphere_profiles_timedep_case",
        ),
    ]
)
def sim_builder_for_write_test(request):
    return request.param


def test_roundtrip_simulation_outputs_identical(machine_conf, sim_builder) -> None:
    """End-to-end roundtrip test, parametrized over all simulation builders."""

    assert_roundtrip_simulation_outputs_identical(sim_builder, machine_conf)


def test_simulation_runs(machine_conf, sim_builder_to_run) -> None:
    """Simulation run + job.001 test, parametrized over all simulation builders."""

    run_simulation_and_check_job(sim_builder_to_run, machine_conf)


def test_files_written(machine_conf, sim_builder_for_write_test) -> None:
    """Tests that the expected files are written, parametrized over all simulation builders."""

    sim_builder_for_write_test(machine_conf)
