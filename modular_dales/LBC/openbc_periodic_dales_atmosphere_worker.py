"""Worker that adds periodic DALES turbulence to atmosphere-based open boundaries."""

from __future__ import annotations

from functools import partial
import glob
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Iterable, Optional, Tuple

import numpy as np
import xarray as xr

from modular_dales.LBC.nest_dales_in_dales.get_all_dales_boundaries import (
    _file_matches_boundary_selection,
    _promote_boundary_stagger_dims,
)
from modular_dales.LBC.openbc_atmosphere_worker import OpenBCAtmosphereWorker

if TYPE_CHECKING:
    from modular_dales.LBC.openbc import do_openboundary

logger = logging.getLogger(__name__)


class OpenBCPeriodicDalesAtmosphereWorker:
    """Compose atmosphere open boundaries with periodic-precursor perturbations."""

    def __init__(self, module: "do_openboundary") -> None:
        self.module = module

    def prepare(self) -> Tuple[xr.Dataset, xr.Dataset]:
        base_worker = OpenBCAtmosphereWorker(self.module)
        boundaries, initfields = base_worker.prepare()

        boundaries = self.apply_perturbations_to_boundaries(boundaries)
        return boundaries, initfields

    def apply_perturbations_to_boundaries(self, boundaries: xr.Dataset) -> xr.Dataset:
        perturbations_by_name, target_time_seconds = (
            self._build_periodic_boundary_perturbations(boundaries)
        )
        return self._add_perturbations(
            boundaries,
            perturbations_by_name,
            target_time_seconds,
        )

    def _periodic_cfg(self):
        cfg = getattr(
            self.module,
            "periodic_dales_turbulence_perturbations",
            None,
        )
        if cfg is not None:
            return cfg
        return getattr(self.module, "nest_in_periodic_dales_and_atmosphere", None)

    def _build_periodic_boundary_perturbations(
        self,
        boundaries: xr.Dataset,
    ) -> Tuple[Dict[str, xr.DataArray], Optional[np.ndarray]]:
        cfg = self._periodic_cfg()
        if cfg is None:
            return {}, None
        outpath = Path(cfg.periodic_outpath)
        perturbation_vars = list(cfg.perturbation_variables or [])
        template_meta = self._collect_template_metadata(boundaries, perturbation_vars)

        yz_files = self._list_crosssection_files(outpath, "crossyz.*.*.nc")
        xz_files = self._list_crosssection_files(outpath, "crossxz.*.*.nc")
        xy_files = self._list_crosssection_files(outpath, "crossxy.*.*.nc")

        # Explicit section selection strategy:
        # - west  boundary <- west edge from crossyz
        # - east  boundary <- center section from crossyz
        # - south boundary <- south edge from crossxz
        # - north boundary <- center section from crossxz
        # - top   boundary <- user-selected top level from crossxy
        #
        # For staggered DALES fields we always select BOTH coordinate families
        # simultaneously (xt/xm, yt/ym, zt/zm). This prevents NaNs for variables
        # living on different staggers (u on xm, v on ym, w on zm).
        yz_west = self._open_selected_cross_dataset(
            yz_files,
            boundary="west",
            sel_index=self._build_lateral_selector(
                yz_files,
                primary="xt",
                stagger="xm",
                mode="edge_low",
            ),
        )
        yz_center = self._open_selected_cross_dataset(
            yz_files,
            boundary="east",
            sel_index=self._build_lateral_selector(
                yz_files,
                primary="xt",
                stagger="xm",
                mode="center",
            ),
        )
        xz_south = self._open_selected_cross_dataset(
            xz_files,
            boundary="south",
            sel_index=self._build_lateral_selector(
                xz_files,
                primary="yt",
                stagger="ym",
                mode="edge_low",
            ),
        )
        xz_center = self._open_selected_cross_dataset(
            xz_files,
            boundary="north",
            sel_index=self._build_lateral_selector(
                xz_files,
                primary="yt",
                stagger="ym",
                mode="center",
            ),
        )
        xy_top = self._open_selected_cross_dataset(
            xy_files,
            boundary="top",
            sel_index=self._build_top_selector(
                xy_files,
                top_layer_index=cfg.top_layer_index,
                primary="zt",
                stagger="zm",
            ),
        )

        source_by_boundary = {
            "west": yz_west,
            "east": yz_center,
            "south": xz_south,
            "north": xz_center,
            "top": xy_top,
        }

        out: Dict[str, xr.DataArray] = {}
        target_time_seconds: Optional[np.ndarray] = None

        for bnd, source_ds in source_by_boundary.items():
            for var in perturbation_vars:
                target_name = f"{var}{bnd}"
                if target_name not in template_meta:
                    continue
                source_var = self._get_source_var(source_ds, var)
                if source_var is None:
                    logger.warning(
                        "Periodic cross-sections do not contain '%s'; skipping perturbations for %s.",
                        var,
                        target_name,
                    )
                    continue

                pert = self._compute_perturbation(source_var)
                pert = self._align_to_template(pert, template_meta[target_name])
                if "time" in pert.coords and target_time_seconds is None:
                    target_time_seconds = self._time_coord_to_seconds(pert["time"])
                out[target_name] = pert.rename(target_name)

        if not out:
            logger.warning(
                "No periodic perturbations were added; check periodic cross-section files and variables."
            )
            return {}, None

        return out, target_time_seconds

    def _cross_chunks(self) -> Optional[dict]:
        if self.module.tchunk is None:
            return None
        return {"time": int(self.module.tchunk)}

    def _list_crosssection_files(self, outpath: Path, pattern: str) -> list[str]:
        files = sorted(glob.glob((outpath / pattern).as_posix()))
        if not files:
            raise FileNotFoundError(
                f"No cross-section files found for pattern '{pattern}' in '{outpath}'."
            )
        return files

    def _open_selected_cross_dataset(
        self,
        files: list[str],
        boundary: str,
        sel_index: Dict[str, float],
    ) -> xr.Dataset:
        selected_files = [
            file_path
            for file_path in files
            if _file_matches_boundary_selection(file_path, sel_index)
        ]

        if not selected_files:
            logger.warning(
                "No cross-section files matched selection %s for boundary %s; falling back to all files.",
                sel_index,
                boundary,
            )
            selected_files = files

        preprocess = partial(_promote_boundary_stagger_dims, boundary=boundary)
        chunks = self._cross_chunks()
        return xr.open_mfdataset(
            selected_files,
            combine="by_coords",
            join="outer",
            chunks=chunks,
            preprocess=preprocess,
        )

    def _collect_template_metadata(
        self,
        boundaries: xr.Dataset,
        perturbation_vars: list[str],
    ) -> Dict[str, Dict[str, object]]:
        metadata: Dict[str, Dict[str, object]] = {}
        for bnd in ["west", "east", "south", "north", "top"]:
            for var in perturbation_vars:
                key = f"{var}{bnd}"
                if key not in boundaries:
                    continue
                da = boundaries[key]
                metadata[key] = {
                    "dims": tuple(da.dims),
                    "coords": {dim: da[dim] for dim in da.dims if dim in da.coords},
                }
        return metadata

    def _coord_unique_sorted(self, ds: xr.Dataset, coord: str) -> np.ndarray:
        if coord not in ds.coords:
            return np.array([])
        values = np.asarray(ds[coord].values, dtype=float).reshape(-1)
        if values.size == 0:
            return np.array([])
        return np.unique(values)

    def _pick_center_value(self, values: np.ndarray) -> float:
        midpoint = 0.5 * (float(values[0]) + float(values[-1]))
        return float(values[np.argmin(np.abs(values - midpoint))])

    def _pick_level_value(
        self, values: np.ndarray, top_layer_index: Optional[int]
    ) -> float:
        if top_layer_index is None:
            return float(values[-1])

        idx = int(top_layer_index)
        if idx < 0:
            idx += values.size
        if idx < 0 or idx >= values.size:
            raise IndexError(
                f"top_layer_index={top_layer_index} is out of bounds for {values.size} levels."
            )
        return float(values[idx])

    def _build_lateral_selector(
        self,
        files: list[str],
        primary: str,
        stagger: str,
        mode: str,
    ) -> Dict[str, float]:
        primary_vals = self._coord_values_from_files(files, primary)
        stagger_vals = self._coord_values_from_files(files, stagger)

        if primary_vals.size == 0 and stagger_vals.size == 0:
            raise ValueError(
                f"Cross-section files are missing both '{primary}' and '{stagger}' coordinates."
            )

        if mode == "edge_low":
            p_target = (
                float(primary_vals[0])
                if primary_vals.size > 0
                else float(stagger_vals[0])
            )
            s_target = float(stagger_vals[0]) if stagger_vals.size > 0 else p_target
        elif mode == "center":
            p_target = (
                self._pick_center_value(primary_vals)
                if primary_vals.size > 0
                else self._pick_center_value(stagger_vals)
            )
            s_target = (
                self._pick_center_value(stagger_vals)
                if stagger_vals.size > 0
                else p_target
            )
        else:
            raise ValueError(f"Unknown section selection mode: {mode}")

        selector: Dict[str, float] = {}
        if primary_vals.size > 0:
            selector[primary] = p_target
        if stagger_vals.size > 0:
            selector[stagger] = s_target
        return selector

    def _build_top_selector(
        self,
        files: list[str],
        top_layer_index: Optional[int],
        primary: str,
        stagger: str,
    ) -> Dict[str, float]:
        z_vals = self._coord_values_from_files(files, primary)
        zm_vals = self._coord_values_from_files(files, stagger)
        if z_vals.size == 0 and zm_vals.size == 0:
            raise ValueError("crossxy files do not expose usable vertical coordinates.")

        z_target = (
            self._pick_level_value(z_vals, top_layer_index)
            if z_vals.size > 0
            else self._pick_level_value(zm_vals, top_layer_index)
        )
        zm_target = (
            self._pick_level_value(zm_vals, top_layer_index)
            if zm_vals.size > 0
            else z_target
        )

        selector: Dict[str, float] = {}
        if z_vals.size > 0:
            selector[primary] = z_target
        if zm_vals.size > 0:
            selector[stagger] = zm_target
        return selector

    def _coord_values_from_files(self, files: list[str], coord: str) -> np.ndarray:
        values: list[np.ndarray] = []
        for file_path in files:
            try:
                with xr.open_dataset(file_path) as ds:
                    if coord not in ds.coords:
                        continue
                    coord_values = np.asarray(ds[coord].values).reshape(-1)
                    if coord_values.size > 0:
                        values.append(coord_values.astype(float))
            except Exception:
                continue

        if not values:
            return np.array([])
        return np.unique(np.concatenate(values))

    def _get_source_var(self, ds: xr.Dataset, var: str) -> Optional[xr.DataArray]:
        candidates = [var]
        if var == "e12":
            candidates = ["e120", "e12"]
        for name in candidates:
            if name in ds:
                return ds[name]
        return None

    def _spatial_horizontal_dims(self, dims: Iterable[str]) -> list[str]:
        return [dim for dim in dims if dim not in {"time", "zt", "zm"}]

    def _mean_spacing(self, coord_values: np.ndarray) -> Optional[float]:
        vals = np.asarray(coord_values, dtype=float).reshape(-1)
        if vals.size < 2:
            return None
        diffs = np.diff(np.sort(vals))
        diffs = diffs[np.abs(diffs) > 0.0]
        if diffs.size == 0:
            return None
        return float(np.mean(np.abs(diffs)))

    def _rolling_window(self, spacing: Optional[float]) -> int:
        if spacing is None or spacing <= 0.0:
            return 1
        cfg = self._periodic_cfg()
        if cfg is None:
            return 1
        scale = float(cfg.filter_scale_m)
        window = max(1, int(round(scale / spacing)))
        if window % 2 == 0:
            window += 1
        return window

    def _time_coord_to_seconds(self, time_coord: xr.DataArray) -> np.ndarray:
        values = np.asarray(time_coord.values)
        if np.issubdtype(values.dtype, np.datetime64):
            time0 = self.module.time0 or self.module.start
            if time0 is None:
                raise ValueError(
                    "Periodic perturbation interpolation needs do_openboundary.time0 or start when source times are datetime-like."
                )
            return (
                values.astype("datetime64[s]") - np.datetime64(time0, "s")
            ) / np.timedelta64(1, "s")
        return values.astype(float)

    def _resample_boundaries_to_time(
        self,
        boundaries: xr.Dataset,
        target_time_seconds: np.ndarray,
    ) -> xr.Dataset:
        if "time" not in boundaries.coords:
            return boundaries

        source_time_seconds = self._time_coord_to_seconds(boundaries["time"])
        boundaries_seconds = boundaries.assign_coords(
            time=("time", source_time_seconds)
        )

        if source_time_seconds.shape == target_time_seconds.shape and np.allclose(
            source_time_seconds, target_time_seconds
        ):
            # Same time axis: keep chunking and avoid building a full interp graph.
            return boundaries_seconds

        out = boundaries_seconds.interp(
            time=target_time_seconds,
            kwargs={"fill_value": "extrapolate"},
            assume_sorted=True,
        )
        return out.assign_coords(time=("time", target_time_seconds))

    def _add_noalign(self, base: xr.DataArray, delta: xr.DataArray) -> xr.DataArray:
        """Add two arrays with identical dims/shape without xarray index alignment."""
        if tuple(base.dims) != tuple(delta.dims):
            raise ValueError(
                f"Dimension mismatch for no-align add: {base.dims} vs {delta.dims}"
            )
        if base.shape != delta.shape:
            raise ValueError(
                f"Shape mismatch for no-align add: {base.shape} vs {delta.shape}"
            )
        return xr.DataArray(
            base.data + delta.data,
            dims=base.dims,
            coords={d: base[d] for d in base.dims if d in base.coords},
            attrs=base.attrs,
            name=base.name,
        )

    def _apply_relaxation_ramp(self, perturbation: xr.DataArray) -> xr.DataArray:
        cfg = self._periodic_cfg()
        tau = None if cfg is None else getattr(cfg, "tau", None)
        if tau is None:
            return perturbation

        tau = float(tau)
        if tau <= 0.0:
            return perturbation

        if "time" not in perturbation.dims:
            return perturbation

        t_seconds = self._time_coord_to_seconds(perturbation["time"])
        t_seconds = np.maximum(t_seconds, 0.0)
        ramp = 1.0 - np.exp(-t_seconds / tau)
        ramp_da = xr.DataArray(
            ramp,
            dims=("time",),
            coords={"time": perturbation["time"]},
        )
        return perturbation * ramp_da

    def _compute_perturbation(self, source_var: xr.DataArray) -> xr.DataArray:
        filtered = source_var
        for dim in self._spatial_horizontal_dims(source_var.dims):
            spacing = None
            if dim in source_var.coords:
                spacing = self._mean_spacing(source_var[dim].values)
            window = self._rolling_window(spacing)
            if window > 1:
                filtered = filtered.rolling(
                    {dim: window}, center=True, min_periods=1
                ).mean()

        return source_var - filtered

    def _align_to_template(
        self,
        source: xr.DataArray,
        template_meta: Dict[str, object],
    ) -> xr.DataArray:
        da = source.squeeze()
        template_dims = template_meta["dims"]
        template_coords = template_meta["coords"]
        stagger_map = {
            "xt": "xm",
            "xm": "xt",
            "yt": "ym",
            "ym": "yt",
            "zt": "zm",
            "zm": "zt",
        }

        for dim in template_dims:
            if dim == "time":
                continue

            target_coord = template_coords.get(dim)
            if target_coord is None:
                continue

            if dim in da.dims:
                da = da.interp(
                    {dim: target_coord},
                    kwargs={"fill_value": "extrapolate"},
                    assume_sorted=True,
                )
                continue

            alt = stagger_map.get(dim)
            if alt in da.dims:
                da = da.interp(
                    {alt: target_coord},
                    kwargs={"fill_value": "extrapolate"},
                    assume_sorted=True,
                )
                da = da.rename({alt: dim})
                continue

            if dim in da.coords:
                da = da.expand_dims({dim: np.asarray(da[dim].values).reshape(-1)})
                da = da.interp(
                    {dim: target_coord},
                    kwargs={"fill_value": "extrapolate"},
                    assume_sorted=True,
                )
                continue

            raise ValueError(
                f"Cannot map perturbation dims {da.dims} onto target dims {template_dims}."
            )

        if "time" in da.dims:
            pass
        else:
            if "time" in template_coords:
                da = da.expand_dims({"time": template_coords["time"].values})
            else:
                raise ValueError(
                    "Target template has no time coordinate while perturbation has no time dimension."
                )

        if tuple(da.dims) != tuple(template_dims):
            if set(da.dims) != set(template_dims):
                raise ValueError(
                    f"Aligned perturbation dims {da.dims} do not match target dims {template_dims}."
                )
            da = da.transpose(*template_dims)
        return da

    def _add_perturbations(
        self,
        boundaries: xr.Dataset,
        perturbations_by_name: Dict[str, xr.DataArray],
        target_time_seconds: Optional[np.ndarray],
    ) -> xr.Dataset:
        if not perturbations_by_name:
            return boundaries

        if target_time_seconds is not None:
            out = self._resample_boundaries_to_time(boundaries, target_time_seconds)
        else:
            out = boundaries.copy()

        for name, pert in perturbations_by_name.items():
            if name in out:
                pert_aligned = pert
                if target_time_seconds is not None and "time" in pert_aligned.coords:
                    pert_seconds = self._time_coord_to_seconds(pert_aligned["time"])
                    # Always normalize perturbation time coords to seconds so xarray
                    # alignment never mixes datetime64 and float indexes.
                    pert_aligned = pert_aligned.assign_coords(
                        time=("time", pert_seconds)
                    )
                    if not (
                        pert_seconds.shape == target_time_seconds.shape
                        and np.allclose(pert_seconds, target_time_seconds)
                    ):
                        pert_aligned = pert_aligned.interp(
                            time=target_time_seconds,
                            kwargs={"fill_value": "extrapolate"},
                            assume_sorted=True,
                        )

                base = out[name]
                if tuple(pert_aligned.dims) != tuple(base.dims):
                    if set(pert_aligned.dims) != set(base.dims):
                        raise ValueError(
                            f"Perturbation dims {pert_aligned.dims} do not match base dims {base.dims} for {name}."
                        )
                    pert_aligned = pert_aligned.transpose(*base.dims)

                pert_aligned = self._apply_relaxation_ramp(pert_aligned)

                out[name] = self._add_noalign(base, pert_aligned)
        return out
