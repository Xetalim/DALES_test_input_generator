from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest

from .helpers import load_machine_conf


@pytest.fixture
def machine_conf(tmp_path: Path):
    """Provide a machine configuration with a temporary BASE_OUTPUT_PATH.

    Returns a function that takes a case name and returns a machine
    configuration dict with BASE_OUTPUT_PATH pointing to a unique
    temporary directory for that case.
    """

    def _machine_conf(tmp_path: Path, casename: str) -> dict:
        conf = deepcopy(load_machine_conf())
        base = tmp_path / f"case_{casename}"
        base.mkdir(parents=True, exist_ok=True)
        conf["case_conf"]["BASE_OUTPUT_PATH"] = str(base)
        return conf

    return lambda casename: _machine_conf(tmp_path, casename)
