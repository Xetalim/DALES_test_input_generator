from __future__ import annotations

import json
import logging
import re
import shutil
import stat
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import f90nml
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import yaml

LOGGER = logging.getLogger("run_standard_cases")


@dataclass
class RunnerSettings:
    # Input/output paths
    cases_file: Path = Path("standard_cases.txt")
    machine_conf: Path = Path("machine_conf.yaml")
    source_root: Path | None = None
    output_root: Path | None = None

    # Run behavior
    # One invocation creates one run folder.
    run_index: int | None = None  # None -> auto detect next run_XXX
    overwrite_existing_run: bool = False
    recreate_figures_only: bool = False
    expnr: int = 1
    numcores: int | None = None
    run_mode: str = "mpi"  # "mpi", "singlecore", or "norun"
    case_filter: str | None = None
    log_level: str = "INFO"

    # Forced baseline output settings
    genstat_dtav: float = 60.0
    genstat_timeav: float = 60.0

    # Dynamic namelist overrides applied on top of staged namoptions
    # Example:
    # {
    #   "run": {"runtime": 1800, "nprocx": 0, "nprocy": 0},
    #   "namchecksim": {"tcheck": 120},
    # }
    namelist_overrides: dict[str, dict[str, Any]] = field(default_factory=dict)

    # NetCDF comparison settings (floating-point tolerant)
    compare_run_subdir: str = "run_001"
    compare_to_baseline_run: int = 1
    atol: float = 1e-8
    rtol: float = 1e-6

    # Plot settings
    plot_variables: list[str] | None = None  # None means auto (all common numeric vars)
    max_variables_per_file: int = 8
    profile_dim_candidates: tuple[str, ...] = (
        "zt",
        "zm",
        "z",
        "lev",
        "level",
        "height",
    )
    time_dim_candidates: tuple[str, ...] = ("time", "t")


# ================================================================
# USER SETTINGS: edit here
# ================================================================
SETTINGS = RunnerSettings(
    run_mode="mpi",
    recreate_figures_only=False,
    # source_root=Path("/Users/andrevanginkel/Documents/20_Code/24_dales_source/dales/cases"),
    # output_root=Path("/Users/andrevanginkel/Documents/40_Input_and_Runs/42_Dales_Cases/42.01_generated_cases/standard_case_regression"),
    namelist_overrides={
        "run": {"nprocx": 0, "nprocy": 0},
        "thermodynamics": {"lbaseexner": True},
    },
    plot_variables=None,
    atol=1e-7,
    rtol=1e-5,
    run_index=3,
)


def enforce_case_size_and_runtime_limits(nml: f90nml.Namelist) -> None:
    """Clamp heavy case settings so regression runs stay bounded.

    Limits enforced:
    - DOMAIN.itot <= 64
    - DOMAIN.jtot <= 64
    - RUN.runtime <= 7200 seconds (2 hours)

    If itot/jtot are reduced, xsize/ysize are scaled proportionally so
    horizontal grid spacing remains unchanged.
    """
    if "DOMAIN" not in nml:
        nml["DOMAIN"] = {}
    if "RUN" not in nml:
        nml["RUN"] = {}

    domain = nml["DOMAIN"]
    run = nml["RUN"]

    def _cap_int_with_scaled_size(
        count_key: str,
        size_key: str,
        max_count: int,
    ) -> None:
        original_count = domain.get(count_key)
        if original_count is None:
            return

        original_count_int = int(original_count)
        if original_count_int <= max_count:
            return

        original_size = domain.get(size_key)
        if original_size is not None:
            # Preserve grid spacing: new_size = old_size * new_count / old_count.
            domain[size_key] = (
                float(original_size) * float(max_count) / float(original_count_int)
            )

        domain[count_key] = max_count
        LOGGER.info(
            "Capped DOMAIN.%s from %s to %s and adjusted DOMAIN.%s",
            count_key,
            original_count_int,
            max_count,
            size_key,
        )

    _cap_int_with_scaled_size("itot", "xsize", max_count=32)
    _cap_int_with_scaled_size("jtot", "ysize", max_count=32)

    runtime = run.get("runtime")
    if runtime is not None and float(runtime) > 7200.0:
        run["runtime"] = 7200
        LOGGER.info("Capped RUN.runtime from %s to 7200", runtime)


def remove_timer_namelist(nml: f90nml.Namelist) -> None:
    """Remove TIMER namelist group if present.

    Some standard cases ship with:
    &timer
        ltimer = .true.
        ltimer_print = .true.
        ltimer_write = .true.
    /
    This function strips that group from staged namoptions.
    """
    timer_keys = [key for key in list(nml.keys()) if str(key).lower() == "timer"]
    for key in timer_keys:
        del nml[key]
        LOGGER.info("Removed namelist group '&%s' from staged namoptions", key)
    timer_keys = [
        key for key in list(nml.keys()) if str(key).lower() == "namcrosssection"
    ]
    for key in timer_keys:
        del nml[key]
        LOGGER.info("Removed namelist group '&%s' from staged namoptions", key)


def load_machine_conf(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def parse_standard_cases(cases_file: Path) -> tuple[Path | None, list[str]]:
    source_root: Path | None = None
    cases: list[str] = []

    lines = cases_file.read_text(encoding="utf-8").splitlines()
    root_re = re.compile(r"^/.*")
    tree_top_re = re.compile(r"^[├└]──\s+(.+)$")
    plain_name_re = re.compile(r"^[A-Za-z0-9_.-]+$")

    for raw_line in lines:
        line = raw_line.rstrip()
        if not line:
            continue

        if source_root is None and root_re.match(line):
            candidate = Path(line)
            if candidate.exists():
                source_root = candidate
                continue

        tree_match = tree_top_re.match(line)
        if tree_match is not None:
            name = tree_match.group(1).strip()
            if name:
                cases.append(name)
            continue

        stripped = line.strip()
        if stripped.startswith("#"):
            continue
        if plain_name_re.match(stripped):
            cases.append(stripped)

    unique_cases = list(dict.fromkeys(cases))
    return source_root, unique_cases


def render_job_file(
    destination_case_dir: Path,
    job_template_path: Path,
    machine_conf: dict,
    numcores_override: int | None,
) -> None:
    content = job_template_path.read_text(encoding="utf-8")
    job_conf = machine_conf.get("job_conf", {})

    input_dir = destination_case_dir / "input"

    required_filetransfers = ""
    for static_input_file in ("rrtmg_lw.nc", "rrtmg_sw.nc"):
        src = input_dir / static_input_file
        if src.exists():
            required_filetransfers += (
                f'ln -sf "$(pwd)/input/{static_input_file}" '
                f'"$RUNDIR/{static_input_file}" || exit\n'
            )

    replacements = {
        "required_folders": "",
        "required_filetransfers": required_filetransfers,
        "required_folder_rsyncs": "",
        "dales_exec": str(job_conf.get("dales_exec", "dales")),
        "debugger": str(job_conf.get("debugger", "lldb")),
        "numcores": str(
            int(numcores_override)
            if numcores_override is not None
            else int(job_conf.get("numcores", 1))
        ),
    }

    for key, value in replacements.items():
        content = content.replace(f"{{{{{key}}}}}", value)

    job_path = destination_case_dir / "job.001"
    job_path.write_text(content, encoding="utf-8")
    mode = job_path.stat().st_mode
    job_path.chmod(mode | stat.S_IXUSR)


def preferred_namoptions_path(input_dir: Path, expnr: int) -> Path:
    # Accept alternative case names like "namoptions-1536.001" as well.
    # Important: prefer those case-specific files over plain namoptions.001,
    # because the latter may come from the generic input template.
    exact_variant_candidates: list[Path] = []
    exact_plain_candidates: list[Path] = []
    exact_expnr_suffix: list[Path] = []
    general_candidates: list[Path] = []

    for path in sorted(input_dir.glob("namoptions*")):
        if not path.is_file():
            continue

        # Require a trailing .NNN experiment suffix.
        if re.match(r"^namoptions.*\.\d{3}$", path.name) is None:
            continue

        general_candidates.append(path)
        if path.name.endswith(f".{expnr:03d}"):
            exact_expnr_suffix.append(path)
            if path.name.startswith("namoptions-"):
                exact_variant_candidates.append(path)
            elif path.name == f"namoptions.{expnr:03d}":
                exact_plain_candidates.append(path)

    # Prefer case-specific files such as namoptions-1536.001.
    if exact_variant_candidates:
        return exact_variant_candidates[0]

    # Then plain namoptions.001.
    if exact_plain_candidates:
        return exact_plain_candidates[0]

    if exact_expnr_suffix:
        return exact_expnr_suffix[0]
    if general_candidates:
        return general_candidates[0]

    raise FileNotFoundError(
        f"No namoptions file found in {input_dir}. "
        "Expected e.g. namoptions.001 or namoptions-XXXX.001"
    )


def apply_namelist_overrides(
    nml: f90nml.Namelist,
    overrides: dict[str, dict[str, Any]],
) -> None:
    for section, section_overrides in overrides.items():
        if section not in nml:
            nml[section] = {}
        for key, value in section_overrides.items():
            nml[section][key] = value


def enforce_namgenstat_and_overrides(input_dir: Path, settings: RunnerSettings) -> None:
    source_namoptions = preferred_namoptions_path(input_dir, settings.expnr)
    nml = f90nml.read(source_namoptions.as_posix())

    remove_timer_namelist(nml)

    if "namgenstat" not in nml:
        nml["namgenstat"] = {}
    nml["namgenstat"]["lstat"] = True
    nml["namgenstat"]["dtav"] = settings.genstat_dtav
    nml["namgenstat"]["timeav"] = settings.genstat_timeav

    if "namnetcdfstats" not in nml:
        nml["namnetcdfstats"] = {}
    nml["namnetcdfstats"]["lsync"] = True

    apply_namelist_overrides(nml, settings.namelist_overrides)
    enforce_case_size_and_runtime_limits(nml)

    target = input_dir / f"namoptions.{settings.expnr:03d}"
    nml.write(target.as_posix(), force=True)


def copy_case_inputs(
    source_case_dir: Path,
    destination_case_dir: Path,
    input_template_dir: Path,
) -> None:
    input_dir = destination_case_dir / "input"
    input_dir.mkdir(parents=True, exist_ok=True)

    if input_template_dir.exists():
        shutil.copytree(input_template_dir, input_dir, dirs_exist_ok=True)

    for entry in source_case_dir.iterdir():
        if entry.name.startswith("results_"):
            continue

        target = input_dir / entry.name
        if entry.is_dir():
            shutil.copytree(entry, target, dirs_exist_ok=True, symlinks=True)
        else:
            # Dereference file symlinks so staged cases remain self-contained.
            if entry.is_symlink():
                try:
                    resolved = entry.resolve(strict=True)
                    if resolved.is_file():
                        shutil.copy2(resolved, target)
                        continue
                except FileNotFoundError:
                    LOGGER.warning("Broken symlink in source case: %s", entry)
            shutil.copy2(entry, target, follow_symlinks=False)


def run_case(case_dir: Path, run_mode: str) -> None:
    cmd = ["./job.001"]
    if run_mode == "norun":
        cmd.append("NORUN")
    elif run_mode == "singlecore":
        cmd.append("singlecore")

    log_file = case_dir / "run.log"
    LOGGER.info("Executing %s in %s", " ".join(cmd), case_dir)
    with log_file.open("w", encoding="utf-8") as f:
        process = subprocess.run(
            cmd,
            cwd=case_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        f.write(process.stdout)
        if process.returncode != 0:
            raise RuntimeError(
                f"Case run failed for {case_dir.name} with exit code {process.returncode}. "
                f"Inspect {log_file}."
            )


def find_relative_netcdf_files(root: Path) -> set[str]:
    return {
        p.relative_to(root).as_posix()
        for p in root.rglob("*.nc")
        if p.is_file() and not p.is_symlink()
    }


def choose_profile_dim(dims: tuple[str, ...], settings: RunnerSettings) -> str | None:
    for candidate in settings.profile_dim_candidates:
        if candidate in dims:
            return candidate
    return None


def choose_time_dim(dims: tuple[str, ...], settings: RunnerSettings) -> str | None:
    for candidate in settings.time_dim_candidates:
        if candidate in dims:
            return candidate
    return None


def _safe_plot_stem(rel_path: str, var_name: str) -> str:
    base = rel_path.replace("/", "__").replace(".", "_")
    return f"{base}__{var_name}"


def _plot_profiles_heatmap_triplet(
    ref_2d: xr.DataArray,
    cur_2d: xr.DataArray,
    rel_path: str,
    var_name: str,
    out_path: Path,
) -> None:
    ref_values = np.asarray(ref_2d.values, dtype=np.float64)
    cur_values = np.asarray(cur_2d.values, dtype=np.float64)
    diff_values = cur_values - ref_values

    # Detrend by removing the time-mean at each vertical level.
    # This highlights relative anomalies and reduces mean-offset dominance.
    ref_detrended = ref_values - np.nanmean(ref_values, axis=0, keepdims=True)
    cur_detrended = cur_values - np.nanmean(cur_values, axis=0, keepdims=True)
    diff_detrended = cur_detrended - ref_detrended

    fig, axes = plt.subplots(
        nrows=3, ncols=2, figsize=(16, 10), constrained_layout=True
    )

    # LEFT COLUMN (kept as before): baseline, current, difference
    # Share color scale between baseline and current.
    finite_ref_cur = np.concatenate(
        [
            ref_values[np.isfinite(ref_values)],
            cur_values[np.isfinite(cur_values)],
        ]
    )
    if finite_ref_cur.size == 0:
        vmin, vmax = -1.0, 1.0
    else:
        vmin = float(np.nanmin(finite_ref_cur))
        vmax = float(np.nanmax(finite_ref_cur))

    im0 = axes[0, 0].imshow(
        ref_values, aspect="auto", origin="lower", vmin=vmin, vmax=vmax
    )
    axes[0, 0].set_title(f"Baseline: {var_name} | File: {rel_path}")
    axes[0, 0].set_ylabel("z-index")
    fig.colorbar(im0, ax=axes[0, 0], orientation="vertical")

    im1 = axes[1, 0].imshow(
        cur_values, aspect="auto", origin="lower", vmin=vmin, vmax=vmax
    )
    axes[1, 0].set_title(f"Current: {var_name} | File: {rel_path}")
    axes[1, 0].set_ylabel("z-index")
    fig.colorbar(im1, ax=axes[1, 0], orientation="vertical")

    finite_diff = diff_values[np.isfinite(diff_values)]
    if finite_diff.size == 0:
        dmax = 1.0
    else:
        dmax = float(np.nanmax(np.abs(finite_diff)))
        if dmax == 0.0:
            dmax = 1.0
    cmap = "RdBu_r"
    vmin = -dmax
    vmax = dmax
    if np.nanmax(diff_values) <= 0.0:
        cmap = "Blues_r"
        vmax = np.nanmax(diff_values)
        vmin = np.nanmin(diff_values)
    if np.nanmin(diff_values) >= 0.0:
        cmap = "Reds"
        vmax = np.nanmax(diff_values)
        vmin = np.nanmin(diff_values)
    im2 = axes[2, 0].imshow(
        diff_values,
        aspect="auto",
        origin="lower",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )
    axes[2, 0].set_title(
        f"Difference (Current - Baseline): {var_name} | File: {rel_path}"
    )
    axes[2, 0].set_ylabel("z-index")
    axes[2, 0].set_xlabel("time-index")
    fig.colorbar(im2, ax=axes[2, 0], orientation="vertical")

    # RIGHT COLUMN: detrended baseline, detrended current, detrended difference
    finite_det = np.concatenate(
        [
            ref_detrended[np.isfinite(ref_detrended)],
            cur_detrended[np.isfinite(cur_detrended)],
        ]
    )
    if finite_det.size == 0:
        det_vmin, det_vmax = -1.0, 1.0
    else:
        det_vmin = float(np.nanmin(finite_det))
        det_vmax = float(np.nanmax(finite_det))

    im3 = axes[0, 1].imshow(
        ref_detrended,
        aspect="auto",
        origin="lower",
        vmin=det_vmin,
        vmax=det_vmax,
    )
    axes[0, 1].set_title(f"Detrended Baseline: {var_name} | File: {rel_path}")
    axes[0, 1].set_ylabel("z-index")
    fig.colorbar(im3, ax=axes[0, 1], orientation="vertical")

    im4 = axes[1, 1].imshow(
        cur_detrended,
        aspect="auto",
        origin="lower",
        vmin=det_vmin,
        vmax=det_vmax,
    )
    axes[1, 1].set_title(f"Detrended Current: {var_name} | File: {rel_path}")
    axes[1, 1].set_ylabel("z-index")
    fig.colorbar(im4, ax=axes[1, 1], orientation="vertical")

    finite_diff_det = diff_detrended[np.isfinite(diff_detrended)]
    if finite_diff_det.size == 0:
        ddet_max = 1.0
    else:
        ddet_max = float(np.nanmax(np.abs(finite_diff_det)))
        if ddet_max == 0.0:
            ddet_max = 1.0

    im5 = axes[2, 1].imshow(
        diff_detrended,
        aspect="auto",
        origin="lower",
        cmap="RdBu_r",
        vmin=-ddet_max,
        vmax=ddet_max,
    )
    axes[2, 1].set_title(
        f"Detrended Difference (Current - Baseline): {var_name} | File: {rel_path}"
    )
    axes[2, 1].set_ylabel("z-index")
    axes[2, 1].set_xlabel("time-index")
    fig.colorbar(im5, ax=axes[2, 1], orientation="vertical")

    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def generate_profiles_heatmaps(
    baseline_profiles_file: Path,
    current_profiles_file: Path,
    rel_path: str,
    heatmap_dir: Path,
    settings: RunnerSettings,
) -> dict[str, Any]:
    heatmap_dir.mkdir(parents=True, exist_ok=True)

    outputs: list[str] = []
    skipped: dict[str, str] = {}

    with xr.open_dataset(
        baseline_profiles_file, decode_times=False
    ) as ds_ref, xr.open_dataset(current_profiles_file, decode_times=False) as ds_cur:
        common_vars = sorted(set(ds_ref.data_vars) & set(ds_cur.data_vars))

        for var_name in common_vars:
            if (
                settings.plot_variables is not None
                and var_name not in settings.plot_variables
            ):
                continue
            if not np.issubdtype(ds_ref[var_name].dtype, np.number):
                continue
            if not np.issubdtype(ds_cur[var_name].dtype, np.number):
                continue

            ref = ds_ref[var_name]
            cur = ds_cur[var_name]

            time_dim = choose_time_dim(tuple(ref.dims), settings)
            profile_dim = choose_profile_dim(tuple(ref.dims), settings)
            if time_dim is None or profile_dim is None:
                skipped[var_name] = "missing_time_or_profile_dim"
                continue
            if time_dim not in cur.dims or profile_dim not in cur.dims:
                skipped[var_name] = "dims_not_shared_between_runs"
                continue

            ref_2d = ref.transpose(time_dim, profile_dim)
            cur_2d = cur.transpose(time_dim, profile_dim)
            if ref_2d.shape != cur_2d.shape:
                skipped[var_name] = "shape_mismatch"
                continue

            out_path = (
                heatmap_dir / f"{_safe_plot_stem(rel_path, var_name)}__heatmap.png"
            )
            _plot_profiles_heatmap_triplet(ref_2d, cur_2d, rel_path, var_name, out_path)
            outputs.append(out_path.as_posix())

    return {
        "profiles_file": rel_path,
        "heatmap_outputs": outputs,
        "skipped_variables": skipped,
    }


def _plot_profile(
    ref_profile: xr.DataArray,
    cur_profile: xr.DataArray,
    profile_dim: str,
    rel_path: str,
    var_name: str,
    out_path: Path,
) -> None:
    plt.figure(figsize=(6, 5))
    plt.plot(ref_profile.values, ref_profile[profile_dim].values, label="baseline")
    plt.plot(cur_profile.values, cur_profile[profile_dim].values, label="current")
    plt.xlabel("value")
    plt.ylabel(profile_dim)
    plt.title(f"Profile: {var_name}\nFile: {rel_path}")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=140)
    plt.close()


def _plot_average(
    ref: xr.DataArray,
    cur: xr.DataArray,
    settings: RunnerSettings,
    rel_path: str,
    var_name: str,
    out_path: Path,
) -> None:
    time_dim = choose_time_dim(tuple(ref.dims), settings)
    plt.figure(figsize=(7, 4))

    if time_dim is not None:
        ref_mean_t = ref.mean(dim=[d for d in ref.dims if d != time_dim], skipna=True)
        cur_mean_t = cur.mean(dim=[d for d in cur.dims if d != time_dim], skipna=True)
        plt.plot(ref_mean_t[time_dim].values, ref_mean_t.values, label="baseline")
        plt.plot(cur_mean_t[time_dim].values, cur_mean_t.values, label="current")
        plt.xlabel(time_dim)
        plt.ylabel("domain mean")
        plt.title(f"Domain Mean Time Series: {var_name}\nFile: {rel_path}")
    else:
        ref_val = float(ref.mean(skipna=True).values)
        cur_val = float(cur.mean(skipna=True).values)
        plt.bar([0, 1], [ref_val, cur_val], tick_label=["baseline", "current"])
        plt.ylabel("domain mean")
        plt.title(f"Domain Mean Comparison: {var_name}\nFile: {rel_path}")

    plt.grid(True, alpha=0.3)
    if time_dim is not None:
        plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=140)
    plt.close()


def compare_netcdf_file(
    baseline_file: Path,
    current_file: Path,
    rel_path: str,
    plots_dir: Path,
    settings: RunnerSettings,
) -> dict[str, Any]:
    with xr.open_dataset(baseline_file, decode_times=False) as ds_ref, xr.open_dataset(
        current_file, decode_times=False
    ) as ds_cur:
        common_vars = sorted(set(ds_ref.data_vars) & set(ds_cur.data_vars))

        numeric_vars: list[str] = []
        for var in common_vars:
            if np.issubdtype(ds_ref[var].dtype, np.number) and np.issubdtype(
                ds_cur[var].dtype, np.number
            ):
                numeric_vars.append(var)

        if settings.plot_variables is not None:
            numeric_vars = [v for v in numeric_vars if v in settings.plot_variables]

        numeric_vars = numeric_vars[: settings.max_variables_per_file]

        var_reports: dict[str, Any] = {}
        changed_vars: list[str] = []
        close_vars: list[str] = []

        for var_name in numeric_vars:
            ref = ds_ref[var_name]
            cur = ds_cur[var_name]

            if ref.shape != cur.shape:
                var_reports[var_name] = {
                    "status": "shape_mismatch",
                    "baseline_shape": list(ref.shape),
                    "current_shape": list(cur.shape),
                }
                changed_vars.append(var_name)
                continue

            ref_vals = np.asarray(ref.values, dtype=np.float64)
            cur_vals = np.asarray(cur.values, dtype=np.float64)
            diff = cur_vals - ref_vals

            mean_abs_diff = float(np.nanmean(np.abs(diff)))
            max_abs_diff = float(np.nanmax(np.abs(diff)))
            rms_diff = float(np.sqrt(np.nanmean(diff**2)))
            baseline_mean = float(np.nanmean(ref_vals))
            current_mean = float(np.nanmean(cur_vals))
            baseline_abs_max = float(np.nanmax(np.abs(ref_vals)))
            tol = settings.atol + settings.rtol * baseline_abs_max
            is_close = bool(max_abs_diff <= tol)

            profile_plot_path: str | None = None
            average_plot_path: str | None = None

            stem = _safe_plot_stem(rel_path, var_name)
            avg_plot = plots_dir / f"{stem}__avg.png"
            _plot_average(ref, cur, settings, rel_path, var_name, avg_plot)
            average_plot_path = avg_plot.as_posix()

            profile_dim = choose_profile_dim(tuple(ref.dims), settings)
            if profile_dim is not None:
                reduce_dims = [d for d in ref.dims if d != profile_dim]
                ref_profile = ref.mean(dim=reduce_dims, skipna=True)
                cur_profile = cur.mean(dim=reduce_dims, skipna=True)

                if ref_profile.shape == cur_profile.shape:
                    prof_plot = plots_dir / f"{stem}__profile.png"
                    _plot_profile(
                        ref_profile,
                        cur_profile,
                        profile_dim,
                        rel_path,
                        var_name,
                        prof_plot,
                    )
                    profile_plot_path = prof_plot.as_posix()

            if is_close:
                close_vars.append(var_name)
            else:
                changed_vars.append(var_name)

            var_reports[var_name] = {
                "status": "compared",
                "is_close": is_close,
                "mean_abs_diff": mean_abs_diff,
                "max_abs_diff": max_abs_diff,
                "rms_diff": rms_diff,
                "baseline_mean": baseline_mean,
                "current_mean": current_mean,
                "baseline_abs_max": baseline_abs_max,
                "tolerance": tol,
                "average_plot": average_plot_path,
                "profile_plot": profile_plot_path,
            }

        return {
            "relative_path": rel_path,
            "baseline_file": baseline_file.as_posix(),
            "current_file": current_file.as_posix(),
            "num_compared_variables": len(var_reports),
            "close_variables": close_vars,
            "changed_variables": changed_vars,
            "all_variables_close": len(changed_vars) == 0,
            "variables": var_reports,
        }


def compare_netcdf_sets(
    baseline_run_dir: Path,
    current_run_dir: Path,
    case_plot_dir: Path,
    settings: RunnerSettings,
) -> dict[str, Any]:
    baseline_files = find_relative_netcdf_files(baseline_run_dir)
    current_files = find_relative_netcdf_files(current_run_dir)

    missing = sorted(baseline_files - current_files)
    extra = sorted(current_files - baseline_files)
    common = sorted(baseline_files & current_files)

    file_reports: dict[str, Any] = {}
    changed_files: list[str] = []
    close_files: list[str] = []
    profiles_heatmaps: dict[str, Any] = {}

    for rel_path in common:
        report = compare_netcdf_file(
            baseline_file=baseline_run_dir / rel_path,
            current_file=current_run_dir / rel_path,
            rel_path=rel_path,
            plots_dir=case_plot_dir,
            settings=settings,
        )
        file_reports[rel_path] = report
        if report["all_variables_close"]:
            close_files.append(rel_path)
        else:
            changed_files.append(rel_path)

        if Path(rel_path).name == "profiles.001.nc":
            heatmap_dir = case_plot_dir / "profiles_heatmaps"
            profiles_heatmaps = generate_profiles_heatmaps(
                baseline_profiles_file=baseline_run_dir / rel_path,
                current_profiles_file=current_run_dir / rel_path,
                rel_path=rel_path,
                heatmap_dir=heatmap_dir,
                settings=settings,
            )

    return {
        "baseline_run_dir": baseline_run_dir.as_posix(),
        "current_run_dir": current_run_dir.as_posix(),
        "counts": {
            "baseline_total_nc": len(baseline_files),
            "current_total_nc": len(current_files),
            "missing_in_current": len(missing),
            "extra_in_current": len(extra),
            "common": len(common),
            "changed_files": len(changed_files),
            "close_files": len(close_files),
        },
        "missing_in_current": missing,
        "extra_in_current": extra,
        "changed_files": changed_files,
        "close_files": close_files,
        "files": file_reports,
        "profiles_heatmaps": profiles_heatmaps,
    }


def filter_cases(cases: list[str], case_filter: str | None) -> list[str]:
    if case_filter is None:
        return cases

    regex = re.compile(case_filter)
    selected = [case for case in cases if regex.search(case)]
    if not selected:
        raise ValueError(f"No cases matched case_filter='{case_filter}'")
    return selected


def resolve_source_root(
    settings: RunnerSettings,
    machine_conf: dict,
    source_root_from_file: Path | None,
) -> Path:
    source_root = settings.source_root or source_root_from_file
    if source_root is None:
        source_root = (
            Path(machine_conf["case_conf"]["SOURCE_PATH"]).expanduser().resolve()
            / "cases"
        )
    source_root = source_root.expanduser().resolve()
    if not source_root.exists():
        raise FileNotFoundError(f"Source root does not exist: {source_root}")
    return source_root


def resolve_output_root(settings: RunnerSettings, machine_conf: dict) -> Path:
    if settings.output_root is not None:
        return settings.output_root.expanduser().resolve()

    return (
        Path(machine_conf["case_conf"]["BASE_OUTPUT_PATH"]).expanduser().resolve()
        / "standard_case_regression_new"
    )


def list_existing_run_indices(output_root: Path) -> list[int]:
    indices: list[int] = []
    for entry in output_root.glob("run_*"):
        if not entry.is_dir():
            continue
        match = re.match(r"run_(\d+)$", entry.name)
        if match is None:
            continue
        indices.append(int(match.group(1)))
    return sorted(indices)


def _report_path(run_dir: Path) -> Path:
    return run_dir / "comparison_report.json"


def load_existing_case_results(run_dir: Path) -> dict[str, Any]:
    report_path = _report_path(run_dir)
    if not report_path.exists():
        return {}

    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        LOGGER.warning(
            "Could not parse existing report %s; treating as empty", report_path
        )
        return {}

    cases = report.get("cases", {})
    if not isinstance(cases, dict):
        return {}
    return cases


def is_run_complete(run_dir: Path, selected_cases: list[str]) -> bool:
    case_results = load_existing_case_results(run_dir)
    if not case_results:
        return False

    for case_name in selected_cases:
        entry = case_results.get(case_name)
        if not isinstance(entry, dict):
            return False
        if entry.get("status") not in {"completed", "staged"}:
            return False
    return True


def resolve_run_index(
    settings: RunnerSettings,
    output_root: Path,
    selected_cases: list[str],
) -> int:
    if settings.run_index is not None:
        return int(settings.run_index)

    existing = list_existing_run_indices(output_root)
    if settings.recreate_figures_only:
        if not existing:
            raise FileNotFoundError(
                f"No existing runs found in {output_root} for recreate_figures_only mode"
            )
        return existing[-1]

    if not existing:
        return 1

    latest = existing[-1]
    latest_dir = output_root / f"run_{latest:03d}"
    if not is_run_complete(latest_dir, selected_cases):
        return latest
    return latest + 1


def write_run_report(
    run_dir: Path,
    run_index: int,
    existing_runs_before_invocation: list[int],
    settings: RunnerSettings,
    per_case_results: dict[str, Any],
) -> None:
    report_path = _report_path(run_dir)
    report_path.write_text(
        json.dumps(
            {
                "run": run_index,
                "existing_runs_before_invocation": existing_runs_before_invocation,
                "settings": {
                    "atol": settings.atol,
                    "rtol": settings.rtol,
                    "run_mode": settings.run_mode,
                    "namelist_overrides": settings.namelist_overrides,
                },
                "cases": per_case_results,
            },
            indent=2,
        ),
        encoding="utf-8",
    )


def main() -> None:
    settings = SETTINGS
    logging.basicConfig(
        level=getattr(logging, settings.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    if settings.run_mode not in {"mpi", "singlecore", "norun"}:
        raise ValueError("run_mode must be one of: 'mpi', 'singlecore', 'norun'")

    machine_conf = load_machine_conf(settings.machine_conf)
    source_root_from_file, parsed_cases = parse_standard_cases(settings.cases_file)

    if not parsed_cases:
        raise ValueError(
            f"No cases found in {settings.cases_file}. "
            "Use either one-case-per-line or tree output with top-level case entries."
        )

    source_root = resolve_source_root(settings, machine_conf, source_root_from_file)
    output_root = resolve_output_root(settings, machine_conf)
    selected_cases = filter_cases(parsed_cases, settings.case_filter)

    input_template_dir = Path("input_template") / "input"
    job_template_path = Path("input_template") / "job.001"
    if not job_template_path.exists():
        raise FileNotFoundError(f"Job template not found: {job_template_path}")

    LOGGER.info("Source root: %s", source_root)
    LOGGER.info("Output root: %s", output_root)
    LOGGER.info("Cases selected: %d", len(selected_cases))
    existing_runs = list_existing_run_indices(output_root)
    run_index = resolve_run_index(settings, output_root, selected_cases)
    run_dir = output_root / f"run_{run_index:03d}"

    if settings.recreate_figures_only:
        if not run_dir.exists():
            raise FileNotFoundError(
                f"Run folder not found for recreate_figures_only mode: {run_dir}"
            )
        LOGGER.info(
            "Recreate-figures-only mode enabled: using existing run %03d in %s",
            run_index,
            run_dir,
        )

    resume_existing = False
    run_dir_non_empty = run_dir.exists() and any(run_dir.iterdir())
    if run_dir_non_empty and not settings.recreate_figures_only:
        if is_run_complete(run_dir, selected_cases):
            if not settings.overwrite_existing_run:
                raise FileExistsError(
                    f"Run folder already exists and appears complete: {run_dir}. "
                    "Set overwrite_existing_run=True or choose another run_index."
                )
        else:
            resume_existing = True

    run_dir.mkdir(parents=True, exist_ok=True)
    LOGGER.info("Existing runs before this invocation: %s", existing_runs)
    if resume_existing:
        LOGGER.info("Resuming incomplete run %03d in %s", run_index, run_dir)
    else:
        LOGGER.info("Starting run %03d in %s", run_index, run_dir)

    per_case_results: dict[str, Any] = load_existing_case_results(run_dir)

    for case_name in selected_cases:
        existing_case_report = per_case_results.get(case_name)
        if (
            not settings.recreate_figures_only
            and isinstance(existing_case_report, dict)
            and existing_case_report.get("status") in {"completed", "staged"}
        ):
            LOGGER.info("Skipping case '%s' (already finished in this run)", case_name)
            continue

        source_case_dir = source_root / case_name
        if not source_case_dir.exists():
            raise FileNotFoundError(
                f"Case '{case_name}' not found under source root: {source_case_dir}"
            )

        destination_case_dir = run_dir / case_name
        destination_case_dir.mkdir(parents=True, exist_ok=True)

        if not settings.recreate_figures_only:
            LOGGER.info("Staging case '%s'", case_name)
            copy_case_inputs(source_case_dir, destination_case_dir, input_template_dir)
            enforce_namgenstat_and_overrides(destination_case_dir / "input", settings)
            render_job_file(
                destination_case_dir,
                job_template_path=job_template_path,
                machine_conf=machine_conf,
                numcores_override=settings.numcores,
            )
            run_case(destination_case_dir, run_mode=settings.run_mode)
        else:
            LOGGER.info(
                "Rebuilding figures for case '%s' from existing outputs", case_name
            )

        case_report: dict[str, Any] = {
            "case": case_name,
            "source": source_case_dir.as_posix(),
            "destination": destination_case_dir.as_posix(),
            "status": (
                "postprocessed"
                if settings.recreate_figures_only
                else ("staged" if settings.run_mode == "norun" else "completed")
            ),
        }

        if run_index >= 2 and (
            settings.recreate_figures_only or settings.run_mode != "norun"
        ):
            baseline_case_dir = (
                output_root / f"run_{settings.compare_to_baseline_run:03d}" / case_name
            )
            baseline_run_subdir = baseline_case_dir / settings.compare_run_subdir
            current_run_subdir = destination_case_dir / settings.compare_run_subdir

            if not baseline_run_subdir.exists():
                LOGGER.warning(
                    "Skipping comparison for case '%s': baseline folder missing: %s",
                    case_name,
                    baseline_run_subdir,
                )
            elif not current_run_subdir.exists():
                raise FileNotFoundError(
                    f"Current run output folder missing: {current_run_subdir}"
                )
            else:
                case_plot_dir = run_dir / "comparison_plots" / case_name
                case_plot_dir.mkdir(parents=True, exist_ok=True)

                comparison = compare_netcdf_sets(
                    baseline_run_dir=baseline_run_subdir,
                    current_run_dir=current_run_subdir,
                    case_plot_dir=case_plot_dir,
                    settings=settings,
                )
                case_report["comparison_to_baseline"] = comparison

        per_case_results[case_name] = case_report

        # Checkpoint after each case so interrupted runs can resume.
        write_run_report(
            run_dir=run_dir,
            run_index=run_index,
            existing_runs_before_invocation=existing_runs,
            settings=settings,
            per_case_results=per_case_results,
        )

    write_run_report(
        run_dir=run_dir,
        run_index=run_index,
        existing_runs_before_invocation=existing_runs,
        settings=settings,
        per_case_results=per_case_results,
    )
    LOGGER.info("Wrote run report: %s", _report_path(run_dir))

    LOGGER.info("All requested runs finished successfully.")


if __name__ == "__main__":
    main()
