import csv
import logging
from pathlib import Path
from typing import Optional

import numpy as np
import xarray as xr

from modular_dales.Surface.LSM.translation_tables.vegetation_properties import (
    top10_to_ifs,
)

logger = logging.getLogger(__name__)


AGS_PARAM_NAMES = [
    "gm25",
    "Ammax25",
    "f0",
    "alpha0",
    "co2_comp298",
    "T1gm",
    "T1Am",
]


_VEG_LOOKUP = {
    "fce": {
        "gm25": 3.0e-3,
        "Ammax25": 2.2,
        "f0": 0.89,
        "a": 0.0,
        "b": 0.0,
        "T1gm": 278.0,
        "T1Am": 281.0,
        "co2_comp298": 82.2,
        "alpha0": 0.017,
    },
    "fbd": {
        "gm25": 2.0e-3,
        "Ammax25": 1.83,
        "f0": 0.80,
        "a": 0.0,
        "b": 0.0,
        "T1gm": 278.0,
        "T1Am": 281.0,
        "co2_comp298": 42.0,
        "alpha0": 0.0142,
    },
    "ara": {
        "gm25": 1.3e-3,
        "Ammax25": 2.20,
        "f0": 0.85,
        "alpha0": 0.0142,
        "a": 2.381,
        "b": 0.6103,
        "T1gm": 278.0,
        "T1Am": 281.0,
        "co2_comp298": 42.0,
    },
    "crp": {
        "gm25": 1.4e-3,
        "Ammax25": 1.83,
        "f0": 0.92,
        "alpha0": 0.0142,
        "a": 2.381,
        "b": 0.6103,
        "T1gm": 278.0,
        "T1Am": 281.0,
        "co2_comp298": 42.0,
    },
    "sem": {
        "gm25": 1.0e-3,
        "Ammax25": 1.83,
        "f0": 0.80,
        "alpha0": 0.0142,
        "a": 2.381,
        "b": 0.6103,
        "T1gm": 278.0,
        "T1Am": 281.0,
        "co2_comp298": 42.0,
    },
    "grs": {
        "C3": {
            "gm25": 1.5e-3,
            "Ammax25": 2.0,
            "f0": 0.75,
            "alpha0": 0.0142,
            "a": 2.381,
            "b": 0.6103,
            "T1gm": 278.0,
            "T1Am": 281.0,
            "co2_comp298": 42.0,
        },
        "C4": {
            "gm25": 2.3e-3,
            "Ammax25": 1.83,
            "f0": 0.70,
            "alpha0": 0.0117,
            "a": 5.323,
            "b": 0.8923,
            "T1gm": 286.15,
            "T1Am": 286.15,
            "co2_comp298": 2.6,
        },
    },
}


def _sorted_on_xy(ds: xr.Dataset, var_name: str) -> xr.DataArray:
    da = ds[var_name]
    if "x" not in da.coords or "y" not in da.coords:
        raise ValueError(f"Dataset variable '{var_name}' must contain x/y coordinates")
    if da["x"].values[0] > da["x"].values[-1]:
        da = da.sortby("x")
    if da["y"].values[0] > da["y"].values[-1]:
        da = da.sortby("y")
    return da


def _sample_on_grid(da: xr.DataArray, grid) -> np.ndarray:
    xq = xr.DataArray(grid.xt, dims=("x",), coords={"x": grid.xt})
    yq = xr.DataArray(grid.yt, dims=("y",), coords={"y": grid.yt})
    sampled = da.interp(x=xq, y=yq, method="nearest")
    return np.asarray(sampled.values)


def apply_top10_to_lsm_writer(
    lsm_writer,
    top10_path: Path,
    fill_north_sea: bool = True,
) -> None:
    ds_lu = xr.open_dataset(top10_path)
    da_lu = _sorted_on_xy(ds_lu, "land_use")
    top10_codes = _sample_on_grid(da_lu, lsm_writer.grid)
    top10_codes_int = np.where(np.isfinite(top10_codes), top10_codes, -9999).astype(int)

    # Reset fractions before re-assigning from the raster map.
    for lu in lsm_writer.lu_types.keys():
        lsm_writer.lu_types[lu]["lu_frac"][:, :] = 0.0

    for lu, lu_info in lsm_writer.lu_types.items():
        if "lu_ids" in lu_info and lu_info["lu_ids"] is not None:
            lu_ids = np.asarray(lu_info["lu_ids"], dtype=int)
            mask = np.isin(top10_codes_int, lu_ids)
        elif "ifs_id" in lu_info:
            ifs_map = np.vectorize(lambda x: top10_to_ifs.get(int(x), -999))(
                top10_codes_int
            )
            mask = ifs_map == int(lu_info["ifs_id"])
        else:
            mask = np.zeros_like(top10_codes, dtype=bool)

        lsm_writer.lu_types[lu]["lu_frac"][mask] = 1.0
        lsm_writer.value_dic[f"c_{lu}"][:, :] = lsm_writer.lu_types[lu]["lu_frac"]

    if fill_north_sea:
        if not hasattr(lsm_writer, "lat") or not hasattr(lsm_writer, "lon"):
            logger.warning(
                "Top10 north-sea fill requested, but no 2D lat/lon arrays are available; skipping"
            )
            fill_north_sea = False

    if fill_north_sea:
        slope = (50.9 - 51.396) / (1.7 - 3.414)
        sea_mask = (lsm_writer.lon < 3.414) & (
            (lsm_writer.lat - 51.396) > slope * (lsm_writer.lon - 3.414)
        )
        if "aqu" in lsm_writer.lu_types:
            for lu in lsm_writer.lu_types.keys():
                lsm_writer.lu_types[lu]["lu_frac"][sea_mask] = 0.0
            lsm_writer.lu_types["aqu"]["lu_frac"][sea_mask] = 1.0
            lsm_writer.value_dic["c_aqu"][:, :] = lsm_writer.lu_types["aqu"]["lu_frac"]

    lsm_writer.init_lutypes_ifs()
    lsm_writer.recalculate_remaining_cover()


class BOFEKInfo:
    def __init__(self, csv_path: Path):
        lookup = {f"B{i}": i + 5 for i in range(1, 19)}
        lookup.update({f"O{i}": i + 23 for i in range(1, 19)})

        soil_rows: dict[int, list[tuple[float, float, str]]] = {}
        with csv_path.open("r", encoding="utf-8") as f:
            reader = csv.reader(f)
            next(reader)  # header
            for row in reader:
                soil_id = int(float(row[0]))
                z_top = float(row[7])
                z_bot = float(row[8])
                soil_code = row[29].strip()
                soil_rows.setdefault(soil_id, []).append((z_top, z_bot, soil_code))

        self.soil_id = np.array(sorted(soil_rows.keys()), dtype=int)
        self.soil_id_lu = np.full(self.soil_id.max() + 1, -9999, dtype=int)
        for i, sid in enumerate(self.soil_id):
            self.soil_id_lu[sid] = i

        max_layers = max(len(v) for v in soil_rows.values())
        n_types = len(self.soil_id)
        self.lookup_index = np.full((n_types, max_layers), -1, dtype=int)
        self.n_layers = np.zeros(n_types, dtype=int)

        for i, sid in enumerate(self.soil_id):
            rows = soil_rows[sid]
            self.n_layers[i] = len(rows)
            for k, (_z_top, _z_bot, soil_code) in enumerate(rows):
                self.lookup_index[i, k] = lookup.get(soil_code, -1)


def apply_bofek_to_lsm_writer(
    lsm_writer,
    bofek_nc_path: Path,
    bofek_csv_path: Path,
) -> None:
    bf = BOFEKInfo(bofek_csv_path)
    ds_soil = xr.open_dataset(bofek_nc_path)
    da_soil = _sorted_on_xy(ds_soil, "bofek_code")
    bf_codes = _sample_on_grid(da_soil, lsm_writer.grid).astype(int)

    nz, ny, nx = lsm_writer.value_dic["index_soil"].shape
    index_soil = np.full((nz, ny, nx), 2, dtype=int)

    for j in range(ny):
        for i in range(nx):
            code = int(bf_codes[j, i])
            if code <= 0 or code >= bf.soil_id_lu.size:
                continue
            table_idx = bf.soil_id_lu[code]
            if table_idx < 0:
                continue

            n_layers = int(bf.n_layers[table_idx])
            if n_layers <= 0:
                continue
            lu_profile = bf.lookup_index[table_idx, :n_layers]
            lu_profile = lu_profile[lu_profile >= 0]
            if lu_profile.size == 0:
                continue

            for k in range(nz):
                src_k = min(k, lu_profile.size - 1)
                index_soil[k, j, i] = int(lu_profile[src_k])

    # Python -> Fortran indexing as expected by DALES lsm input
    lsm_writer.value_dic["index_soil"][:, :, :] = index_soil + 1


def get_veg_params(lu: str, planttype: Optional[int] = None) -> dict[str, float]:
    lu = lu.lower()
    if lu not in _VEG_LOOKUP:
        logger.warning("Unknown AGS land-use type '%s'; using zero parameters", lu)
        return {k: 0.0 for k in AGS_PARAM_NAMES}

    if lu == "grs":
        key = "C4" if planttype == 4 else "C3"
        params = _VEG_LOOKUP[lu][key].copy()
    else:
        params = _VEG_LOOKUP[lu].copy()

    return {k: float(params.get(k, 0.0)) for k in AGS_PARAM_NAMES}


def apply_ags_parameters_to_lsm_writer(
    lsm_writer,
    grass_planttype: Optional[int] = None,
) -> None:
    shape = (lsm_writer.grid.jtot, lsm_writer.grid.itot)

    for par in AGS_PARAM_NAMES:
        for lu, lu_info in lsm_writer.lu_types.items():
            lu_short = str(lu_info.get("lu_short", lu)).lower()
            planttype = grass_planttype if lu_short == "grs" else None
            veg_params = get_veg_params(lu_short, planttype=planttype)

            field_name = f"{par}_{lu}"
            lsm_writer.value_dic[field_name] = np.full(
                shape, veg_params[par], dtype=float
            )
            if field_name not in lsm_writer.fields:
                lsm_writer.fields.append(field_name)
