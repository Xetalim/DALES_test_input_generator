"""Helpers for resolving external DALES data with a fixed upstream source.

External datasets are always sourced from:
- Repository: https://github.com/dalesteam/dales
- Branch: dev

RRTMG assets are initialized via the DALES submodule layout.
"""

from __future__ import annotations

from dataclasses import dataclass
import logging
import pathlib
import shutil
import subprocess

logger = logging.getLogger(__name__)

DALES_REPO_URL = "https://github.com/dalesteam/dales"
DALES_REPO_BRANCH = "dev"
RRTMG_LW_PATH = pathlib.Path("external/RRTMG/RRTMG_LW")
RRTMG_SW_PATH = pathlib.Path("external/RRTMG/RRTMG_SW")
RRTMGP_DATA_PATH = pathlib.Path("external/rrtmgp-data")
SPARSE_PATHS = [
    "data",
    RRTMG_LW_PATH.as_posix(),
    RRTMG_SW_PATH.as_posix(),
    RRTMGP_DATA_PATH.as_posix(),
]


@dataclass(frozen=True)
class ExternalDataPaths:
    """Concrete cached paths to external DALES assets."""

    rrtmg_lw: pathlib.Path
    rrtmg_sw: pathlib.Path
    rrtmgp_data_dir: pathlib.Path
    van_genuchten_parameters: pathlib.Path


@dataclass(frozen=True)
class RrtmgDataPaths:
    """Concrete cached paths to RRTMG/RRTMGP assets."""

    rrtmg_lw: pathlib.Path
    rrtmg_sw: pathlib.Path
    rrtmgp_data_dir: pathlib.Path


def cache_root(sim) -> pathlib.Path:
    """Return the root cache directory for generated helper artifacts."""
    root = pathlib.Path.cwd() / "COG_CACHE" / "external_data"
    root.mkdir(parents=True, exist_ok=True)
    return root


def _run_git(args: list[str], cwd: pathlib.Path) -> None:
    try:
        subprocess.run(
            ["git", *args],
            cwd=cwd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"git {' '.join(args)} failed in {cwd}: {exc.stderr.strip() or exc.stdout.strip()}"
        ) from exc


def _dales_checkout_dir(sim) -> pathlib.Path:
    return cache_root(sim) / "dales_dev"


def _ensure_dales_dev_checkout(sim) -> pathlib.Path:
    repo = _dales_checkout_dir(sim)
    repo.parent.mkdir(parents=True, exist_ok=True)

    if not repo.exists():
        _run_git(
            [
                "clone",
                "--filter=blob:none",
                "--depth",
                "1",
                "--no-checkout",
                DALES_REPO_URL,
                repo.as_posix(),
            ],
            cwd=repo.parent,
        )

    _run_git(["checkout", DALES_REPO_BRANCH], cwd=repo)
    _run_git(["sparse-checkout", "init", "--cone"], cwd=repo)
    _run_git(["sparse-checkout", "set", *SPARSE_PATHS], cwd=repo)
    _run_git(["read-tree", "-mu", "HEAD"], cwd=repo)

    _run_git(["submodule", "sync", "--recursive"], cwd=repo)
    _run_git(
        [
            "submodule",
            "update",
            "--init",
            "--recursive",
            "--depth",
            "1",
            RRTMG_LW_PATH.as_posix(),
            RRTMG_SW_PATH.as_posix(),
        ],
        cwd=repo,
    )
    # Keep full clone for rrtmgp-data as requested.
    _run_git(
        [
            "submodule",
            "update",
            "--init",
            "--recursive",
            RRTMGP_DATA_PATH.as_posix(),
        ],
        cwd=repo,
    )
    _stage_attribution_files(repo, cache_root(sim) / "rrtmg")
    return repo


def _stage_attribution_files(repo: pathlib.Path, target_dir: pathlib.Path) -> None:
    """Copy README/LICEN?E markdown docs from RRTMG sources into cache metadata."""
    target_dir.mkdir(parents=True, exist_ok=True)
    roots = [
        repo / RRTMG_LW_PATH,
        repo / RRTMG_SW_PATH,
        repo / RRTMGP_DATA_PATH,
    ]
    for root in roots:
        if not root.exists():
            continue
        for pattern in ("README*", "LICENSE*"):
            for source in root.glob(pattern):
                if source.exists() and source.is_file():
                    shutil.copy2(source, target_dir / f"{root.name}_{source.name}")


def _require_file(path: pathlib.Path, hint: str) -> pathlib.Path:
    if not path.exists() or not path.is_file():
        raise FileNotFoundError(f"Missing required file: {path}. {hint}")
    return path


def _require_dir(path: pathlib.Path, hint: str) -> pathlib.Path:
    if not path.exists() or not path.is_dir():
        raise FileNotFoundError(f"Missing required directory: {path}. {hint}")
    return path


def resolve_external_data_paths(sim) -> ExternalDataPaths:
    """Resolve radiation and LSM datasets from fixed DALES dev checkout."""
    repo = _ensure_dales_dev_checkout(sim)
    hint = "Expected DALES dev checkout with initialized radiation submodules."
    rrtmg_lw = _require_file(
        repo / "external" / "RRTMG" / "RRTMG_LW" / "data" / "rrtmg_lw.nc",
        hint,
    )
    rrtmg_sw = _require_file(
        repo / "external" / "RRTMG" / "RRTMG_SW" / "data" / "rrtmg_sw.nc",
        hint,
    )
    rrtmgp_data_dir = _require_dir(repo / "external" / "rrtmgp-data", hint)
    van_genuchten = _require_file(
        repo / "data" / "van_genuchten_parameters.nc",
        "Ensure the DALES repository includes the tracked data directory.",
    )

    return ExternalDataPaths(
        rrtmg_lw=rrtmg_lw,
        rrtmg_sw=rrtmg_sw,
        rrtmgp_data_dir=rrtmgp_data_dir,
        van_genuchten_parameters=van_genuchten,
    )


def resolve_rrtmg_data_paths(sim) -> RrtmgDataPaths:
    """Resolve only RRTMG/RRTMGP datasets needed by radiation schemes."""
    repo = _ensure_dales_dev_checkout(sim)
    hint = "Expected DALES dev checkout with initialized radiation submodules."
    rrtmg_lw = _require_file(
        repo / "external" / "RRTMG" / "RRTMG_LW" / "data" / "rrtmg_lw.nc",
        hint,
    )
    rrtmg_sw = _require_file(
        repo / "external" / "RRTMG" / "RRTMG_SW" / "data" / "rrtmg_sw.nc",
        hint,
    )
    rrtmgp_data_dir = _require_dir(repo / "external" / "rrtmgp-data", hint)

    return RrtmgDataPaths(
        rrtmg_lw=rrtmg_lw,
        rrtmg_sw=rrtmg_sw,
        rrtmgp_data_dir=rrtmgp_data_dir,
    )


def resolve_van_genuchten_path(sim) -> pathlib.Path:
    """Resolve only the LSM van Genuchten parameter file."""
    repo = _ensure_dales_dev_checkout(sim)
    return _require_file(
        repo / "data" / "van_genuchten_parameters.nc",
        "Ensure the DALES repository includes the tracked data directory.",
    )


def stage_backrad_file(
    sim, source_file: pathlib.Path, target_name: str
) -> pathlib.Path:
    """Copy backrad input into cache and return staged path."""
    root = cache_root(sim) / "backrad"
    root.mkdir(parents=True, exist_ok=True)
    staged = root / target_name
    shutil.copy2(source_file, staged)
    return staged
