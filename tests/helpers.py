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
    problems = []

    def compare_dirs(left_dir: Path, right_dir: Path) -> None:
        left_entries = {entry.name: entry for entry in left_dir.iterdir()}
        right_entries = {entry.name: entry for entry in right_dir.iterdir()}

        left_only = sorted(set(left_entries) - set(right_entries))
        right_only = sorted(set(right_entries) - set(left_entries))
        if left_only:
            problems.append(f"Only in {left_dir}: {left_only}")
        if right_only:
            problems.append(f"Only in {right_dir}: {right_only}")

        common_names = sorted(set(left_entries) & set(right_entries))
        differing_files = []
        problematic_names = []

        for name in common_names:
            left_entry = left_entries[name]
            right_entry = right_entries[name]

            left_is_dir = left_entry.is_dir()
            right_is_dir = right_entry.is_dir()
            left_is_file = left_entry.is_file()
            right_is_file = right_entry.is_file()

            if left_is_dir and right_is_dir:
                compare_dirs(left_entry, right_entry)
                continue

            if left_is_file and right_is_file:
                if not filecmp.cmp(left_entry, right_entry, shallow=False):
                    differing_files.append(name)
                continue

            problematic_names.append(name)

        if differing_files:
            problems.append(
                f"Differing files between {left_dir} and {right_dir}: {differing_files}"
            )
        if problematic_names:
            problems.append(
                f"Problematic files between {left_dir} and {right_dir}: {problematic_names}"
            )

    compare_dirs(left, right)

    if problems:
        raise ValueError("\n".join(problems))
