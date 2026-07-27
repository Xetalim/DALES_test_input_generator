import numpy as np
import xarray as xr

from modular_dales.Surface.LSM.LCZ.get_from_LCZ import lcz_dict, map_lcz_to_fields


def test_map_lcz_to_fields_maps_values_for_valid_lcz_cells():
    lcz_da = xr.DataArray(
        np.array([[1.0, 1.0], [2.0, 2.0]], dtype=float),
        dims=("y", "x"),
    )
    ifs_da = xr.DataArray(
        np.array([[20.0, 10.0], [20.0, 30.0]], dtype=float),
        dims=("y", "x"),
    )

    fields = map_lcz_to_fields(lcz_da, ifs_da, lcz_dict)
    f_bld = fields["building_surface_fraction"].values
    z0_urb = fields["roughness_length"].values

    assert np.isclose(f_bld[0, 0], lcz_dict[1].building_surface_fraction)
    assert np.isclose(f_bld[1, 0], lcz_dict[2].building_surface_fraction)
    assert np.isclose(f_bld[0, 1], lcz_dict[1].building_surface_fraction)
    assert np.isclose(f_bld[1, 1], lcz_dict[2].building_surface_fraction)

    assert np.isclose(z0_urb[0, 1], lcz_dict[1].roughness_length)
    assert np.isclose(z0_urb[1, 1], lcz_dict[2].roughness_length)
