from __future__ import annotations

import filecmp
from pathlib import Path

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MACHINE_CONF_PATH = PROJECT_ROOT / "machine_conf.yaml"


def load_machine_conf() -> dict:
    """Load the machine configuration YAML as a dictionary."""
    with MACHINE_CONF_PATH.open("r") as f:
        return yaml.safe_load(f)


def assert_dirs_equal(left: Path, right: Path) -> None:
    """Recursively assert that two directories have identical contents.

    Raises a ValueError describing the differences if any are found.
    """
    cmp = filecmp.dircmp(left, right, shallow=False)

    problems = []
    if cmp.left_only:
        problems.append(f"Only in {left}: {sorted(cmp.left_only)}")
    if cmp.right_only:
        problems.append(f"Only in {right}: {sorted(cmp.right_only)}")
    if cmp.diff_files:
        problems.append(
            f"Differing files between {left} and {right}: {sorted(cmp.diff_files)}"
        )
    if cmp.funny_files:
        problems.append(
            f"Problematic files between {left} and {right}: {sorted(cmp.funny_files)}"
        )

    for subdir in cmp.common_dirs:
        assert_dirs_equal(left / subdir, right / subdir)

    if problems:
        raise ValueError("\n".join(problems))
