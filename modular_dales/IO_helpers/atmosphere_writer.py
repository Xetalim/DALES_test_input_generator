import os
from typing import List, Any, TYPE_CHECKING

import matplotlib.pyplot as plt
import netCDF4
import numpy as np

if TYPE_CHECKING:
    from modular_dales.vars import VariableDefinition


class AtmosphereProfileWriter:
    """Helper class for writing atmosphere profiles and time-dependent forcings."""

    def write_base_profiles(
        self,
        z: np.ndarray,
        var_dic: dict["VariableDefinition", Any],
        output_path,
        exp_id: int,
    ) -> None:
        """Write base (time-independent) profiles to init.<exp_id>.nc.

        var_dic is expected to map variable names to containers that have:
          - .definition.long_name
          - .definition.unit
          - .values (1D array along z), or None to skip.
        """
        with netCDF4.Dataset(
            output_path / f"init.{exp_id:03d}.nc", "w", format="NETCDF3_CLASSIC"
        ) as ncout:
            ncout.createDimension("zh", len(z))
            nc_z_var = ncout.createVariable("zh", "f8", ("zh",))
            nc_z_var[:] = z
            for var_definition, var_container in var_dic.items():
                if getattr(
                    var_container.definition, "must_only_be_time_dependent", False
                ):
                    raise ValueError(
                        f"Variable '{var_definition.name}' is marked as must_only_be_time_dependent but is being written to base profile file"
                    )
                if getattr(var_container, "values", None) is not None:
                    ncout.createVariable(var_definition.name, "f8", ("zh",))[
                        :
                    ] = var_container.values
                    ncout.variables[var_definition.name].long_name = (
                        var_definition.long_name
                    )
                    ncout.variables[var_definition.name].units = var_definition.unit

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
