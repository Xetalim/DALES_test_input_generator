from modular_dales.LBC.openbc import (
    do_openboundary,
    Nest_in_Harmonie,
    Nest_in_KNMI,
    Periodic_Dales_Turbulence_Perturbations,
)


def test_periodic_addon_can_be_combined_with_harmonie_source():
    mod = do_openboundary()
    mod += Nest_in_Harmonie(ml_glob="ml", sfc_glob="sfc")
    mod += Periodic_Dales_Turbulence_Perturbations(
        periodic_outpath="/tmp/periodic",
        filter_scale_m=2000.0,
    )

    assert mod.nest_in_harmonie is not None
    assert mod.periodic_dales_turbulence_perturbations is not None
    assert (
        mod.periodic_dales_turbulence_perturbations.periodic_outpath == "/tmp/periodic"
    )


def test_periodic_addon_can_be_combined_with_knmi_source():
    mod = do_openboundary()
    mod += Nest_in_KNMI(ml_glob="ml", sfc_glob="sfc")
    mod += Periodic_Dales_Turbulence_Perturbations(
        periodic_outpath="/tmp/periodic",
        top_layer_index=-1,
    )

    assert mod.nest_in_knmi is not None
    assert mod.periodic_dales_turbulence_perturbations is not None
    assert mod.periodic_dales_turbulence_perturbations.top_layer_index == -1
