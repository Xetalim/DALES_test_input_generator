import os
from typing import List, Any, TYPE_CHECKING, Mapping

import matplotlib.pyplot as plt
import netCDF4
import numpy as np

if TYPE_CHECKING:
    from modular_dales.vars import VariableDefinition


class AtmosphereProfileWriter:
    """Helper class for writing atmosphere profiles and time-dependent forcings."""

    def _ensure_init_axes(
        self,
        ncout,
        z: np.ndarray,
        times: np.ndarray | None = None,
    ) -> None:
        if "zh" not in ncout.dimensions:
            ncout.createDimension("zh", len(z))
        elif len(ncout.dimensions["zh"]) != len(z):
            raise ValueError(
                f"Existing init file has zh size {len(ncout.dimensions['zh'])}, expected {len(z)}"
            )

        if "zh" in ncout.variables:
            z_var = ncout.variables["zh"]
        else:
            z_var = ncout.createVariable("zh", "f8", ("zh",))
        z_var[:] = z

        if times is None:
            return

        if "time" not in ncout.dimensions:
            ncout.createDimension("time", len(times))
        elif len(ncout.dimensions["time"]) != len(times):
            raise ValueError(
                f"Existing init file has time size {len(ncout.dimensions['time'])}, expected {len(times)}"
            )

        if "time" in ncout.variables:
            t_var = ncout.variables["time"]
        else:
            t_var = ncout.createVariable("time", "f8", ("time",))
        t_var[:] = times
        t_var.long_name = "time"
        t_var.units = "s"

    def write_base_profiles(
        self,
        z: np.ndarray,
        var_dic: dict["VariableDefinition", Any],
        output_path,
        exp_id: int,
        times: List[float] | None = None,
        time_profile_series: (
            Mapping["VariableDefinition", Mapping[float, np.ndarray]] | None
        ) = None,
        extra_time_height_fields: Mapping[str, Mapping[str, Any]] | None = None,
    ) -> None:
        """Write base (time-independent) profiles to init.<exp_id>.nc.

        var_dic is expected to map variable names to containers that have:
          - .definition.long_name
          - .definition.unit
          - .values (1D array along z), or None to skip.
        """
        init_path = output_path / f"init.{exp_id:03d}.nc"
        with netCDF4.Dataset(init_path, "w", format="NETCDF3_CLASSIC") as ncout:
            times_arr = np.asarray(times, dtype=float) if times is not None else None
            if extra_time_height_fields and times_arr is None:
                raise ValueError(
                    "write_base_profiles requires 'times' when 'extra_time_height_fields' is provided"
                )
            self._ensure_init_axes(ncout, z, times_arr)

            time_profile_varnames = set()
            if time_profile_series:
                if times_arr is None:
                    raise ValueError(
                        "write_base_profiles requires 'times' when 'time_profile_series' is provided"
                    )
                for var_definition, series in time_profile_series.items():
                    series_times = sorted(float(t) for t in series.keys())
                    if len(series_times) != len(times_arr) or not np.allclose(
                        np.asarray(series_times, dtype=float),
                        times_arr,
                    ):
                        raise ValueError(
                            f"Timed init profile for '{var_definition.name}' has times {series_times}, expected {times_arr.tolist()}"
                        )

                    stacked = np.asarray(
                        [
                            np.asarray(series[float(t)], dtype=float)
                            for t in series_times
                        ],
                        dtype=float,
                    )
                    if stacked.ndim != 2 or stacked.shape[1] != len(z):
                        raise ValueError(
                            f"Timed init profile for '{var_definition.name}' must have shape (time, zh); got {stacked.shape}"
                        )

                    if var_definition.name in ncout.variables:
                        var = ncout.variables[var_definition.name]
                    else:
                        var = ncout.createVariable(
                            var_definition.name, "f8", ("time", "zh")
                        )
                    var[:, :] = stacked
                    var.long_name = var_definition.long_name + " (time-dependent)"
                    var.units = var_definition.unit
                    time_profile_varnames.add(var_definition.name)

            for var_definition, var_container in var_dic.items():
                if getattr(
                    var_container.definition, "must_only_be_time_dependent", False
                ):
                    raise ValueError(
                        f"Variable '{var_definition.name}' is marked as must_only_be_time_dependent but is being written to base profile file"
                    )
                if var_definition.name in time_profile_varnames:
                    continue
                if getattr(var_container, "values", None) is not None:
                    if var_definition.name in ncout.variables:
                        var = ncout.variables[var_definition.name]
                    else:
                        var = ncout.createVariable(var_definition.name, "f8", ("zh",))
                    var[:] = var_container.values
                    var.long_name = var_definition.long_name
                    var.units = var_definition.unit

            if extra_time_height_fields:
                for name, spec in extra_time_height_fields.items():
                    values = np.asarray(spec.get("values"), dtype=float)
                    if values.ndim != 2:
                        raise ValueError(
                            f"Init field '{name}' must be 2D; got shape {values.shape}"
                        )
                    if values.shape == (len(z), len(times_arr)):
                        values = values.T
                    if values.shape != (len(times_arr), len(z)):
                        raise ValueError(
                            f"Init field '{name}' must have shape {(len(times_arr), len(z))} or {(len(z), len(times_arr))}; got {values.shape}"
                        )

                    var = ncout.createVariable(name, "f8", ("time", "zh"))
                    var[:, :] = values
                    var.long_name = str(spec.get("long_name", name))
                    var.units = str(spec.get("units", "1"))

    def write_time_dependent_profiles(
        self,
        z: np.ndarray,
        times_with_zero: List[float],
        timedep_var_dic: dict["VariableDefinition", np.ndarray],
        output_path,
        exp_id: int,
    ) -> None:
        """Write time-dependent forcings to forcings.<exp_id>.nc.

        timedep_var_dic is expected to map variable-definition-like keys to
        arrays of shape (time,) or (z, time+1), including the t=0 profile
        which will be stripped before writing.
        The key objects must provide:
          - .time_dependent_name
          - .long_name
          - .unit
        """
        times = times_with_zero[1:]  # Exclude time=0 profile, which is in the base file
        with netCDF4.Dataset(
            output_path / f"forcings.{exp_id:03d}.nc",
            "w",
            format="NETCDF3_CLASSIC",
        ) as ncout:
            ncout.createDimension("zh", len(z))
            ncout.createDimension("time", len(times))

            nc_z_var = ncout.createVariable("zh", "f8", ("zh",))
            nc_z_var.long_name = "Height of vertical levels"
            nc_z_var.units = "m"
            nc_z_var[:] = z

            nc_t_var = ncout.createVariable("time", "f8", ("time",))
            nc_t_var.long_name = "Time validity of time-dependent values"
            nc_t_var.units = "s"
            nc_t_var[:] = np.asarray(times, dtype=float)

            for var_definition, var_container in timedep_var_dic.items():
                if var_container is None:
                    continue
                arr_np = np.asarray(var_container, dtype=float)
                if arr_np.ndim == 1:
                    # Exclude time=0 entry
                    arr_np = arr_np[1:]
                    if len(arr_np) != len(times):
                        raise ValueError(
                            f"Timed variable '{var_definition.time_dependent_name}' has wrong length {len(arr_np)}; expected {len(times)}"
                        )
                    nc_var = ncout.createVariable(
                        var_definition.time_dependent_name, "f8", ("time",)
                    )
                    nc_var[:] = arr_np
                    nc_var.long_name = var_definition.long_name
                    nc_var.units = var_definition.unit
                elif arr_np.ndim == 2:
                    # Exclude time=0 profile along time axis
                    arr_np = arr_np[:, 1:]
                    if arr_np.shape != (len(z), len(times)):
                        raise ValueError(
                            f"Timed variable '{var_definition.time_dependent_name}' has wrong shape {arr_np.shape}; expected {(len(z), len(times))}"
                        )
                    nc_var = ncout.createVariable(
                        var_definition.time_dependent_name, "f8", ("zh", "time")
                    )
                    nc_var[:] = arr_np
                    nc_var.long_name = var_definition.long_name + " (time-dependent)"
                    nc_var.units = var_definition.unit
                else:
                    raise ValueError(
                        f"Timed variable '{var_definition.time_dependent_name}' must be 1D or 2D; got shape {arr_np.shape}"
                    )

    def plot_profiles(
        self,
        z: np.ndarray,
        var_dic: dict["VariableDefinition", np.ndarray],
        output_path,
        exp_id: int,
    ) -> None:
        """Plot base profiles for quick inspection."""
        os.makedirs(output_path / ".." / "profiles", exist_ok=True)
        for var_definition, var_container in var_dic.items():
            if getattr(var_container, "values", None) is None:
                continue
            fig, ax = plt.subplots()
            ax.plot(var_container.values, z)
            label = var_definition.long_name
            ax.set_xlabel(label)
            ax.set_ylabel("z (m)")
            ax.set_title(f"exp {exp_id:03d} {label}")
            fig.savefig(
                output_path / ".." / "profiles" / f"profile_{var_definition.name}.png",
                dpi=300,
            )
            plt.close(fig)
