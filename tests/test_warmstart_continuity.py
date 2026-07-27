from __future__ import annotations

import copy
import shutil
import subprocess
from pathlib import Path
from typing import Any

import f90nml
import numpy as np
import pytest
import xarray as xr

from modular_dales.Configuration.output_modules import (
    BulkMicrophysicsStatisticsOutputModule,
    ColumnStatisticsOutputModule,
    IndependentOutputModule,
    NetCDFStatisticsSyncModule,
    ParticlesOutputModule,
    QuadrantStatisticsOutputModule,
    SamplingTendencyOutputModule,
    StatTendencyOutputModule,
    StressStatisticsOutputModule,
    TiltStatisticsOutputModule,
    VariableBudgetOutputModule,
)
from modular_dales.Configuration.physics_modules import (
    BulkMicrophysicsSettingsModule,
    TracerSettingsModule,
)

from .sim_builders.test_warmstart import build_warmstart_base_sim


def _dales_exec_exists(machine_conf: dict[str, Any]) -> bool:
    exe = machine_conf.get("job_conf", {}).get("dales_exec")
    return bool(exe and Path(exe).exists())


def _run_job(case_dir: Path, timeout_seconds: int, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["./job.001", *args],
        cwd=str(case_dir),
        text=True,
        capture_output=True,
        timeout=timeout_seconds,
        check=False,
    )


def _likely_cause_hint(msg: str) -> str:
    lower = msg.lower()
    if "time" in lower and "mismatch" in lower:
        return "Scheduler drift around dtav/timeav/tnextwrite; verify cadence and restart boundary alignment."
    if "missing" in lower or "extra" in lower:
        return "Incomplete restart state or conditional allocation path differs across split and continuous runs."
    if "variable" in lower:
        return "Module accumulator/counter state likely not fully persisted into restart files."
    return "Check module namelist toggles, restart files, and optional tracer dependencies."


def _time_coord_name(ds: xr.Dataset) -> str | None:
    for candidate in ("time", "t", "Time", "TIME"):
        if candidate in ds.coords:
            return candidate
    for dim in ds.dims:
        if str(dim).lower() == "time":
            return str(dim)
    return None


def _relative_diff(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    denom = np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-30)
    return np.abs(a - b) / denom


def _compare_netcdf_file(
    cont_path: Path,
    split_path: Path,
    restart_time: float,
    abs_tol: float,
    rel_tol: float,
    include_restart_time: bool,
) -> list[str]:
    with xr.open_dataset(cont_path, engine="netcdf4") as dcont, xr.open_dataset(split_path, engine="netcdf4") as dsplit:
        cont = dcont.load()
        split = dsplit.load()

    failures: list[str] = []
    tcoord = _time_coord_name(cont)
    if tcoord and tcoord in split.coords:
        cont_times = np.asarray(cont[tcoord].values)
        split_times = np.asarray(split[tcoord].values)
        mask_cont = cont_times >= restart_time if include_restart_time else cont_times > restart_time
        mask_split = split_times >= restart_time if include_restart_time else split_times > restart_time
        cont_post = cont_times[mask_cont]
        split_post = split_times[mask_split]
        if cont_post.shape != split_post.shape or not np.array_equal(cont_post, split_post):
            failures.append(
                f"{cont_path.name}: time grid mismatch post-restart "
                f"(continuous={cont_post.shape[0]} records, split={split_post.shape[0]} records). "
                f"hint={_likely_cause_hint('time mismatch')}"
            )
            return failures

        cont = cont.sel({tcoord: cont_post})
        split = split.sel({tcoord: split_post})

    common_vars = sorted(set(cont.data_vars).intersection(set(split.data_vars)))
    if not common_vars:
        failures.append(f"{cont_path.name}: no overlapping data variables to compare")
        return failures

    for var_name in common_vars:
        a = np.asarray(cont[var_name].values)
        b = np.asarray(split[var_name].values)
        if a.shape != b.shape:
            failures.append(
                f"{cont_path.name}: variable={var_name} shape mismatch continuous={a.shape} split={b.shape}. "
                f"hint={_likely_cause_hint('variable shape mismatch')}"
            )
            continue

        exact = np.array_equal(a, b)
        close = np.allclose(a, b, atol=abs_tol, rtol=rel_tol, equal_nan=True)
        if exact or close:
            continue

        abs_diff = np.abs(a - b)
        rel_diff = _relative_diff(a, b)
        max_abs = float(np.nanmax(abs_diff))
        max_rel = float(np.nanmax(rel_diff))

        first_time_msg = "first_time=unknown"
        if tcoord and tcoord in cont[var_name].dims:
            axis = cont[var_name].dims.index(tcoord)
            times = np.asarray(cont[tcoord].values)
            reduced_axes = tuple(i for i in range(a.ndim) if i != axis)
            exceeded = np.any((abs_diff > abs_tol) & (rel_diff > rel_tol), axis=reduced_axes)
            idx = np.where(exceeded)[0]
            if idx.size > 0:
                first_time_msg = f"first_time={float(times[idx[0]])}"

        failures.append(
            f"{cont_path.name}: variable={var_name} mismatch (exact={exact}, within_tol={close}, "
            f"max_abs={max_abs:.6e}, max_rel={max_rel:.6e}, {first_time_msg}). "
            f"hint={_likely_cause_hint('variable mismatch')}"
        )

    return failures


def _compare_outputs(
    cont_run_dir: Path,
    split_run_dir: Path,
    patterns: list[str],
    restart_time: float,
    abs_tol: float,
    rel_tol: float,
    include_restart_time: bool,
) -> list[str]:
    failures: list[str] = []

    cont_files: set[str] = set()
    split_files: set[str] = set()
    for pattern in patterns:
        cont_files.update(fp.name for fp in cont_run_dir.glob(pattern) if fp.is_file())
        split_files.update(fp.name for fp in split_run_dir.glob(pattern) if fp.is_file())

    missing = sorted(cont_files - split_files)
    extra = sorted(split_files - cont_files)
    if missing:
        failures.append(
            f"Missing files in split run: {missing}. hint={_likely_cause_hint('missing file')}"
        )
    if extra:
        failures.append(
            f"Extra files in split run: {extra}. hint={_likely_cause_hint('extra file')}"
        )

    for rel in sorted(cont_files.intersection(split_files)):
        cpath = cont_run_dir / rel
        spath = split_run_dir / rel
        if cpath.suffix.lower() != ".nc":
            continue
        failures.extend(
            _compare_netcdf_file(
                cpath,
                spath,
                restart_time=restart_time,
                abs_tol=abs_tol,
                rel_tol=rel_tol,
                include_restart_time=include_restart_time,
            )
        )

    return failures


DEFAULTS = {
    "domain": {"itot": 16, "jtot": 16, "kmax": 48, "xsize": 160.0, "ysize": 160.0},
    "timing": {"tfinal": 360, "trestart": 180, "dtav": 60, "timeav": 60},
    "tolerances": {"abs": 0.0, "rel": 0.0},
    "include_restart_time": True,
}


def _build_case_modules(module_name: str, timing: dict[str, Any], domain: dict[str, Any]) -> list[Any]:
    dtav = float(timing["dtav"])
    timeav = float(timing["timeav"])

    modules: list[Any] = [
        TracerSettingsModule(nsv=0),
        NetCDFStatisticsSyncModule(lsync=True),
        IndependentOutputModule(
            stats_enabled=True,
            stats_dtav=dtav,
            stats_timeav=timeav,
            timestat_enabled=True,
            timestat_dtav=dtav,
            budget_enabled=True,
            budget_dtav=dtav,
            budget_timeav=timeav,
        ),
    ]

    if module_name == "modvarbudget":
        modules.extend(
            [
                TracerSettingsModule(nsv=2),
                VariableBudgetOutputModule(enabled=True, dtav=dtav, timeav=timeav),
            ]
        )
    elif module_name == "modsampling":
        modules.append(SamplingTendencyOutputModule(dtav=dtav, timeav=timeav))
    elif module_name == "modquadrant":
        modules.append(
            QuadrantStatisticsOutputModule(
                enabled=True,
                dtav=dtav,
                timeav=timeav,
                hole=0.0,
                iwind=1,
                klow=2,
                khigh=int(domain["kmax"]),
            )
        )
    elif module_name == "modsamptend":
        modules.extend(
            [
                BulkMicrophysicsSettingsModule(imicro=2, l_sb=True, nsv=2),
                SamplingTendencyOutputModule(dtav=dtav, timeav=timeav),
            ]
        )
    elif module_name == "modstattend":
        modules.append(StatTendencyOutputModule(enabled=True, dtav=dtav, timeav=timeav))
    elif module_name == "modcolstat":
        modules.append(
            ColumnStatisticsOutputModule(
                enabled=True,
                npoints=2,
                x_idx=[3, 7],
                y_idx=[3, 7],
            )
        )
    elif module_name == "modbulkmicrostat3":
        modules.extend(
            [
                BulkMicrophysicsSettingsModule(imicro=2, l_sb=True, nsv=2),
                BulkMicrophysicsStatisticsOutputModule(enabled=True, dtav=dtav, timeav=timeav),
            ]
        )
    elif module_name == "modtilt":
        modules.append(TiltStatisticsOutputModule(enabled=True, dtav=dtav, timeav=timeav))
    elif module_name == "modstress":
        modules.append(StressStatisticsOutputModule(enabled=True, dtav=dtav, timeav=timeav))
    elif module_name == "modparticles":
        modules.append(ParticlesOutputModule(enabled=True, dtav=dtav, timeav=timeav))
    else:
        raise ValueError(f"Unknown warmstart module case: {module_name}")

    return modules


MODULE_CASES = [
    {
        "name": "modvarbudget",
        "patterns": ["*varbudget*", "*budget*"],
    },
    {
        "name": "modsampling",
        "mpi_layouts": [(1, 1), (2, 1)],
        "patterns": ["*sampling*", "*samptend*"],
    },
    {
        "name": "modquadrant",
        "patterns": ["*quadrant*"],
    },
    {
        "name": "modsamptend",
        "patterns": ["*samptend*", "*sampling*"],
    },
    {
        "name": "modstattend",
        "patterns": ["*stattend*"],
    },
    {
        "name": "modcolstat",
        "patterns": ["*colstat*"],
    },
    {
        "name": "modbulkmicrostat3",
        "patterns": ["*bulkmicro*", "*microstat*"],
    },
    {
        "name": "modtilt",
        "patterns": ["*tilt*"],
    },
    {
        "name": "modstress",
        "patterns": ["*stress*"],
    },
    {
        "name": "modparticles",
        "patterns": ["*particle*"],
    },
]


def _run_single_layout(
    *,
    module_case: dict[str, Any],
    base_machine_conf: dict[str, Any],
    tmp_path: Path,
    domain: dict[str, Any],
    timing: dict[str, Any],
    tolerances: dict[str, Any],
    include_restart_time: bool,
    layout: tuple[int, int],
) -> list[str]:
    nprocx, nprocy = int(layout[0]), int(layout[1])
    np_ranks = max(1, nprocx * nprocy)

    module_name = str(module_case["name"])
    case_base_name = f"warmstart_{module_name}_np{np_ranks}"
    case_conf = copy.deepcopy(base_machine_conf)
    case_conf["job_conf"]["numcores"] = np_ranks
    case_conf["case_conf"]["BASE_OUTPUT_PATH"] = str(tmp_path)

    tfinal = int(timing["tfinal"])
    trestart = int(timing["trestart"])

    cont_case_dir = tmp_path / f"{case_base_name}_continuous"
    split_case_dir = tmp_path / f"{case_base_name}_split"
    for case_dir in (cont_case_dir, split_case_dir):
        if case_dir.exists():
            shutil.rmtree(case_dir)

    cont_sim = build_warmstart_base_sim(
        case_conf,
        cont_case_dir.name,
        domain=domain,
        tfinal=tfinal,
        trestart=trestart,
        lwarmstart=False,
        extra_modules=_build_case_modules(module_name, timing, domain),
    )
    cont_sim.sim_preprocessing_pipeline()
    cont_case_dir = Path(cont_sim.output_path)

    split_runtime = trestart + 1
    split_sim = build_warmstart_base_sim(
        case_conf,
        split_case_dir.name,
        domain=domain,
        tfinal=split_runtime,
        trestart=trestart,
        lwarmstart=False,
        extra_modules=_build_case_modules(module_name, timing, domain),
    )
    split_sim.sim_preprocessing_pipeline()
    split_case_dir = Path(split_sim.output_path)

    timeout_seconds = 1200

    # 1) Continuous reference run (0 -> Tfinal)
    cont_res = _run_job(cont_case_dir, timeout_seconds)
    if cont_res.returncode != 0:
        return [
            f"{module_name} np={np_ranks}: continuous run failed rc={cont_res.returncode}\n"
            f"stdout:\n{cont_res.stdout[-4000:]}\n"
            f"stderr:\n{cont_res.stderr[-4000:]}"
        ]

    # 2) First split segment must be cold start (0 -> Trestart)
    split_input_nml = split_case_dir / "input" / "namoptions.001"
    split_input_run = f90nml.read(str(split_input_nml)).get("RUN", {})
    if bool(split_input_run.get("lwarmstart", False)):
        return [
            f"{module_name} np={np_ranks}: split segment-1 namelist unexpectedly has lwarmstart=True"
        ]

    split_res_1 = _run_job(split_case_dir, timeout_seconds)
    if split_res_1.returncode != 0:
        return [
            f"{module_name} np={np_ranks}: split segment-1 failed rc={split_res_1.returncode}\n"
            f"stdout:\n{split_res_1.stdout[-4000:]}\n"
            f"stderr:\n{split_res_1.stderr[-4000:]}"
        ]

    split_run_dir = split_case_dir / "run_001"
    if not any(split_run_dir.glob("initd*")):
        return [
            f"{module_name} np={np_ranks}: restart artifact missing after split segment-1"
        ]

    # 3) Switch the existing split run namelist in-place to warmstart mode.
    split_run_nml = split_run_dir / "namoptions.001"
    split_run_nml_data = f90nml.read(str(split_run_nml))
    if "RUN" not in split_run_nml_data:
        split_run_nml_data["RUN"] = {}
    split_run_nml_data["RUN"]["runtime"] = float(tfinal)
    split_run_nml_data["RUN"]["lwarmstart"] = True
    split_run_nml_data.write(str(split_run_nml), force=True)

    split_run_cfg = f90nml.read(str(split_run_nml)).get("RUN", {})
    if not bool(split_run_cfg.get("lwarmstart", False)):
        return [
            f"{module_name} np={np_ranks}: split segment-2 namelist missing lwarmstart=True"
        ]

    # 4) Second split segment warmstart (Trestart -> Tfinal)
    split_res_2 = _run_job(split_case_dir, timeout_seconds, "WARMSTART")
    if split_res_2.returncode != 0:
        return [
            f"{module_name} np={np_ranks}: split segment-2 warmstart failed rc={split_res_2.returncode}\n"
            f"stdout:\n{split_res_2.stdout[-4000:]}\n"
            f"stderr:\n{split_res_2.stderr[-4000:]}"
        ]

    cont_run_dir = cont_case_dir / "run_001"
    patterns = list(module_case.get("patterns", ["*.nc"]))
    failures = _compare_outputs(
        cont_run_dir=cont_run_dir,
        split_run_dir=split_run_dir,
        patterns=patterns,
        restart_time=float(trestart),
        abs_tol=float(tolerances["abs"]),
        rel_tol=float(tolerances["rel"]),
        include_restart_time=include_restart_time,
    )

    return [f"{module_name} np={np_ranks}: {msg}" for msg in failures]


@pytest.mark.parametrize("module_case", MODULE_CASES, ids=[c["name"] for c in MODULE_CASES])
def test_warmstart_continuity_module(module_case: dict[str, Any], machine_conf, tmp_path: Path) -> None:
    conf = machine_conf(f"warmstart_{module_case['name']}")
    if not _dales_exec_exists(conf):
        pytest.skip("DALES executable in machine_conf.yaml is not available")

    domain = copy.deepcopy(DEFAULTS["domain"])
    domain.update(module_case.get("domain", {}))
    assert int(domain["itot"]) <= 32 and int(domain["jtot"]) <= 32 and int(domain["kmax"]) <= 96
    timing = copy.deepcopy(DEFAULTS["timing"])
    timing.update(module_case.get("timing", {}))
    tolerances = copy.deepcopy(DEFAULTS["tolerances"])
    tolerances.update(module_case.get("tolerances", {}))
    include_restart_time = bool(module_case.get("include_restart_time", DEFAULTS["include_restart_time"]))

    assert 0 < int(timing["trestart"]) < int(timing["tfinal"])

    layouts = module_case.get("mpi_layouts", [(1, 1)])
    all_failures: list[str] = []
    for layout in layouts:
        all_failures.extend(
            _run_single_layout(
                module_case=module_case,
                base_machine_conf=conf,
                tmp_path=tmp_path,
                domain=domain,
                timing=timing,
                tolerances=tolerances,
                include_restart_time=include_restart_time,
                layout=(int(layout[0]), int(layout[1])),
            )
        )

    assert not all_failures, "\n".join(all_failures)