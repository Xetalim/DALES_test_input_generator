"""Backrad profile helpers for DALES radiation schemes."""

from __future__ import annotations

from dataclasses import dataclass
import pathlib
from typing import Optional, Union

import numpy as np
from scipy.interpolate import interp1d
import xarray as xr


@dataclass
class BackradPressureProfile:
    """Pressure-coordinate background sounding for radiation input.

    Pressure is provided in Pascal and should be ordered from surface to top
    (high to low pressure).
    """

    pressure_pa: list[float]
    temperature_k: list[float]
    specific_humidity_kgkg: list[float]
    ozone_kgkg: Optional[list[float]] = None
    liquid_water_kgkg: Optional[list[float]] = None

    def normalized(self) -> "BackradPressureProfile":
        p = np.asarray(self.pressure_pa, dtype=float)
        t = np.asarray(self.temperature_k, dtype=float)
        q = np.asarray(self.specific_humidity_kgkg, dtype=float)

        if not (p.size == t.size == q.size):
            raise ValueError(
                "pressure_pa, temperature_k, and specific_humidity_kgkg must have equal length"
            )
        if p.size < 2:
            raise ValueError("Backrad profile must contain at least 2 pressure levels")
        if np.any(p <= 0):
            raise ValueError("pressure_pa must be positive")

        o3 = None
        if self.ozone_kgkg is not None:
            o3 = np.asarray(self.ozone_kgkg, dtype=float)
            if o3.size != p.size:
                raise ValueError("ozone_kgkg length must match pressure_pa")
        lwc = None
        if self.liquid_water_kgkg is not None:
            lwc = np.asarray(self.liquid_water_kgkg, dtype=float)
            if lwc.size != p.size:
                raise ValueError("liquid_water_kgkg length must match pressure_pa")

        order = np.argsort(p)[::-1]
        p = p[order]
        t = t[order]
        q = q[order]
        if o3 is not None:
            o3 = o3[order]
        if lwc is not None:
            lwc = lwc[order]

        return BackradPressureProfile(
            pressure_pa=p.tolist(),
            temperature_k=t.tolist(),
            specific_humidity_kgkg=q.tolist(),
            ozone_kgkg=None if o3 is None else o3.tolist(),
            liquid_water_kgkg=None if lwc is None else lwc.tolist(),
        )

    def to_netcdf_dataset(self) -> xr.Dataset:
        profile = self.normalized()
        ds = xr.Dataset(
            data_vars={
                "T": ("lev", np.asarray(profile.temperature_k, dtype=float)),
                "q": ("lev", np.asarray(profile.specific_humidity_kgkg, dtype=float)),
            },
            coords={"lev": ("lev", np.asarray(profile.pressure_pa, dtype=float))},
        )
        if profile.ozone_kgkg is not None:
            ds["o3"] = ("lev", np.asarray(profile.ozone_kgkg, dtype=float))
        return ds

    def to_ascii(self) -> str:
        profile = self.normalized()
        pressure = np.asarray(profile.pressure_pa, dtype=float)
        temp = np.asarray(profile.temperature_k, dtype=float)
        humidity = np.asarray(profile.specific_humidity_kgkg, dtype=float)
        ozone = (
            np.asarray(profile.ozone_kgkg, dtype=float)
            if profile.ozone_kgkg is not None
            else np.zeros_like(pressure)
        )
        liquid = (
            np.asarray(profile.liquid_water_kgkg, dtype=float)
            if profile.liquid_water_kgkg is not None
            else np.zeros_like(pressure)
        )

        lines = [f"{temp[0]:.6f} {pressure.size}"]
        for p, t, q, o3, ql in zip(pressure, temp, humidity, ozone, liquid):
            lines.append(f"{p:.6f} {t:.6f} {q:.8e} {o3:.8e} {ql:.8e}")
        lines.append("")
        return "\n".join(lines)

    @classmethod
    def from_netcdf(cls, path: pathlib.Path) -> "BackradPressureProfile":
        ds = xr.open_dataset(path)
        profile = cls(
            pressure_pa=np.asarray(ds["lev"].values, dtype=float).tolist(),
            temperature_k=np.asarray(ds["T"].values, dtype=float).tolist(),
            specific_humidity_kgkg=np.asarray(ds["q"].values, dtype=float).tolist(),
            ozone_kgkg=(
                np.asarray(ds["o3"].values, dtype=float).tolist()
                if "o3" in ds
                else None
            ),
        )
        return profile.normalized()

    @classmethod
    def from_ascii(cls, path: pathlib.Path) -> "BackradPressureProfile":
        with path.open("r", encoding="utf-8") as handle:
            rows = [line.strip() for line in handle.readlines() if line.strip()]

        if len(rows) < 2:
            raise ValueError(f"Invalid ASCII backrad file: {path}")

        first = rows[0].split()
        if len(first) < 2:
            raise ValueError(f"Invalid ASCII backrad header in {path}")
        nlev = int(float(first[1]))

        p: list[float] = []
        t: list[float] = []
        q: list[float] = []
        o3: list[float] = []
        ql: list[float] = []

        for row in rows[1 : 1 + nlev]:
            cols = row.split()
            if len(cols) < 3:
                raise ValueError(f"Invalid ASCII backrad row in {path}: {row}")
            p.append(float(cols[0]))
            t.append(float(cols[1]))
            q.append(float(cols[2]))
            o3.append(float(cols[3]) if len(cols) > 3 else 0.0)
            ql.append(float(cols[4]) if len(cols) > 4 else 0.0)

        return cls(
            pressure_pa=p,
            temperature_k=t,
            specific_humidity_kgkg=q,
            ozone_kgkg=o3,
            liquid_water_kgkg=ql,
        ).normalized()


@dataclass
class BackradInterpolatedProfile:
    """Interpolated-profile style builder for backrad pressure profiles.

    Provide anchor points on pressure levels and interpolate onto a target
    pressure grid. If ``target_pressure_pa`` is omitted, the anchor pressure
    levels are used as output levels.
    """

    pressure_pa: list[float]
    temperature_points: list[float]
    specific_humidity_points: list[float]
    ozone_points: Optional[list[float]] = None
    liquid_water_points: Optional[list[float]] = None
    target_pressure_pa: Optional[list[float]] = None
    fill_value: Union[float, str] = "extrapolate"

    def _interpolate(
        self, x: np.ndarray, y: np.ndarray, x_target: np.ndarray
    ) -> np.ndarray:
        order = np.argsort(x)
        x_sorted = x[order]
        y_sorted = y[order]
        return interp1d(
            x_sorted,
            y_sorted,
            fill_value=self.fill_value,
            bounds_error=False,
        )(x_target)

    def to_profile(
        self,
        template_profile: Optional[BackradPressureProfile] = None,
    ) -> BackradPressureProfile:
        p = np.asarray(self.pressure_pa, dtype=float)
        t = np.asarray(self.temperature_points, dtype=float)
        q = np.asarray(self.specific_humidity_points, dtype=float)

        if not (p.size == t.size == q.size):
            raise ValueError(
                "pressure_pa, temperature_points, and specific_humidity_points must have equal length"
            )
        if p.size < 2:
            raise ValueError(
                "BackradInterpolatedProfile requires at least two anchor points"
            )

        if self.target_pressure_pa is not None:
            target_p = np.asarray(self.target_pressure_pa, dtype=float)
        elif template_profile is not None:
            target_p = np.asarray(template_profile.pressure_pa, dtype=float)
        else:
            target_p = p.copy()

        t_out = self._interpolate(p, t, target_p)
        q_out = self._interpolate(p, q, target_p)

        o3_out = None
        if self.ozone_points is not None:
            o3 = np.asarray(self.ozone_points, dtype=float)
            if o3.size != p.size:
                raise ValueError("ozone_points length must match pressure_pa")
            o3_out = self._interpolate(p, o3, target_p)

        ql_out = None
        if self.liquid_water_points is not None:
            ql = np.asarray(self.liquid_water_points, dtype=float)
            if ql.size != p.size:
                raise ValueError("liquid_water_points length must match pressure_pa")
            ql_out = self._interpolate(p, ql, target_p)

        return BackradPressureProfile(
            pressure_pa=target_p.tolist(),
            temperature_k=t_out.tolist(),
            specific_humidity_kgkg=q_out.tolist(),
            ozone_kgkg=None if o3_out is None else o3_out.tolist(),
            liquid_water_kgkg=None if ql_out is None else ql_out.tolist(),
        ).normalized()


def profile_from_path(path: pathlib.Path) -> BackradPressureProfile:
    if path.suffix == ".nc":
        return BackradPressureProfile.from_netcdf(path)
    return BackradPressureProfile.from_ascii(path)


def default_profile() -> BackradPressureProfile:
    """Return a default profile, preferring pre-existing repository files."""
    candidates = (
        pathlib.Path.cwd() / "extra_data" / "backrad.inp.001.nc",
        pathlib.Path.cwd() / "extra_data" / "backrad.inp.001",
    )
    for candidate in candidates:
        if candidate.exists():
            return profile_from_path(candidate)

    pressure = [100000.0, 92500.0, 85000.0, 70000.0, 50000.0, 30000.0, 10000.0]
    temperature = [290.0, 285.0, 279.0, 265.0, 250.0, 230.0, 210.0]
    humidity = [1.2e-2, 8.0e-3, 4.0e-3, 1.2e-3, 3.0e-4, 8.0e-5, 1.0e-5]
    ozone = [3.0e-8, 3.0e-8, 4.0e-8, 1.2e-7, 3.5e-7, 7.0e-7, 1.0e-6]
    return BackradPressureProfile(
        pressure_pa=pressure,
        temperature_k=temperature,
        specific_humidity_kgkg=humidity,
        ozone_kgkg=ozone,
    )


def write_profile(profile: BackradPressureProfile, path: pathlib.Path) -> pathlib.Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".nc":
        profile.to_netcdf_dataset().to_netcdf(path)
    else:
        path.write_text(profile.to_ascii(), encoding="utf-8")
    return path
