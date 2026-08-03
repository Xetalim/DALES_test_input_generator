from types import SimpleNamespace

import pytest

from modular_dales.LBC.periodic_precursor_crosssections import (
    PeriodicPrecursorCrossSections,
)


class DummyGrid:
    def __init__(self, itot: int, jtot: int, kmax: int):
        self.itot = itot
        self.jtot = jtot
        self.kmax = kmax


@pytest.mark.parametrize(
    "top_layer_index, expected_height",
    [
        (0, [1]),
        (-1, [8]),
        (3, [4]),
    ],
)
def test_periodic_precursor_crosssections_sets_edges_middle_and_top(
    top_layer_index: int,
    expected_height: list[int],
):
    sim = SimpleNamespace(
        grid=DummyGrid(itot=10, jtot=12, kmax=8),
        nml={},
        nml_docs={},
    )

    module = PeriodicPrecursorCrossSections(sim=sim, top_layer_index=top_layer_index)
    module.prepare_calculation()

    assert sim.nml["namcrosssection"]["crossortho"] == [1, 6, 10]
    assert sim.nml["namcrosssection"]["crossplane"] == [1, 7, 12]
    assert sim.nml["namcrosssection"]["crossheight"] == expected_height


def test_periodic_precursor_crosssections_rejects_invalid_top_layer():
    sim = SimpleNamespace(
        grid=DummyGrid(itot=4, jtot=4, kmax=3),
        nml={},
        nml_docs={},
    )

    module = PeriodicPrecursorCrossSections(sim=sim, top_layer_index=3)
    with pytest.raises(IndexError):
        module.prepare_calculation()
