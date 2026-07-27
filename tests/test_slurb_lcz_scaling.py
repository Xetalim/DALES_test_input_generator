import numpy as np
import xarray as xr

from modular_dales.Geometry.GridDales import GridDales
from modular_dales.Surface.LSM.LSM_output_dales import LSM_output_dales
from modular_dales.Surface.LSM.SLuRB.slurb import slbCreatorClass


def test_apply_slurb_parameters_lcz_scales_f_bld_with_urban_cover():
    grid = GridDales(
        itot=2,
        jtot=2,
        kmax=1,
        xsize=200.0,
        ysize=200.0,
        kmax_soil=1,
        xlat=52.0,
        xlon=5.0,
        x0=0.0,
        y0=0.0,
        alpha=1.0,
        dz0=10.0,
    )

    lu_types = {
        "urb": {
            "lu_long": "Urban",
            "lu_short": "urb",
            "ifs_id": 20,
            "lveg": False,
            "laqu": False,
        }
    }

    lsm = LSM_output_dales(grid=grid, lu_types=lu_types, soil_levels=1)
    lsm.value_dic["cover_urb"][:, :] = np.array([[0.0, 0.2], [0.5, 1.0]])

    lsm.LCZ_ds = xr.Dataset(
        {
            "building_surface_fraction": (("y", "x"), np.array([[0.8, 0.5], [0.4, 0.3]])),
            "aspect_ratio": (("y", "x"), np.array([[1.0, 1.0], [1.0, 1.0]])),
            "height_roughness_element": (("y", "x"), np.array([[10.0, 10.0], [10.0, 10.0]])),
            "roughness_length": (("y", "x"), np.array([[0.3, 0.4], [0.5, 0.6]])),
        }
    )

    slb_generator = slbCreatorClass(grid)
    lsm.apply_slurb_parameters_lcz(slb_generator)

    expected_f_bld = np.array([[0.0, 0.1], [0.2, 0.3]])
    assert np.allclose(slb_generator.f_bld, expected_f_bld)
    assert np.allclose(slb_generator.z0_urb, lsm.LCZ_ds["roughness_length"].values)
