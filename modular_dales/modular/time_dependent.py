from dataclasses import dataclass, field, fields, is_dataclass
from numbers import Number
from typing import TYPE_CHECKING, Any, Mapping, Optional

import numpy as np

from modular_dales.Atmosphere.ls2d_atmosphere import FromLS2D, LS2DAtmosphereModule
from modular_dales.Configuration import TimeModule
from modular_dales.IO_helpers.atmosphere_writer import AtmosphereProfileWriter
from modular_dales.MODULE_REGISTRY import register_module
from modular_dales.modular.time_dependent_scalars import TimeDependentScalar
from modular_dales.vars import get_var_by_name, VariableDefinition

from .simulation_module import simulation_module

if TYPE_CHECKING:
    from . import dales_simulation


@register_module
@dataclass
class TimedependentModule(simulation_module):
    """Collects and writes all time-dependent forcing series for atmosphere files."""

    sim: Optional["dales_simulation"] = field(default=None, repr=False)
    timesteps: list[float] = field(default_factory=list)
    timedep_var_dic: dict[VariableDefinition, np.ndarray] = field(
        default_factory=dict,
        init=False,
        repr=False,
        metadata={"serialize": False},
    )
    profile_writer: AtmosphereProfileWriter = field(
        default_factory=AtmosphereProfileWriter,
        init=False,
        repr=False,
        metadata={"serialize": False},
    )
    ntimedep: int = field(
        default=0,
        init=False,
        repr=False,
        metadata={
            "serialize": False,
            "nml": "PHYSICS",
            "key": "ntimedep",
        },
    )
    ltimedep: bool = field(
        default=True,
        init=True,
        repr=False,
        metadata={
            "serialize": False,
            "nml": "PHYSICS",
            "key": "ltimedep",
        },
    )
    usesLS2DforTime: FromLS2D = field(
        default=None,
        init=True,
        repr=False,
        metadata={"serialize": False},
    )
    LS2Dmodule: Optional["LS2DAtmosphereModule"] = field(
        default=None,
        init=False,
        repr=False,
        metadata={"serialize": False},
    )

    def __post_init__(self):
        super().__init__(self.sim)
        if isinstance(self.timesteps, np.ndarray):
            self.timesteps = self.timesteps.tolist()
        self.module_name = "TimedependentModule"

    def _merge_forcing_series(
        self,
        merged: dict[VariableDefinition, dict[float, Any]],
        var_definition: VariableDefinition,
        series: Mapping[float, Any],
    ):
        if var_definition not in merged:
            merged[var_definition] = {}
        for time_value, value in series.items():
            t = float(time_value)
            if t in merged[var_definition]:
                old_value = np.asarray(merged[var_definition][t], dtype=float)
                new_value = np.asarray(value, dtype=float)
                if not np.allclose(old_value, new_value):
                    raise ValueError(
                        f"Conflicting timed forcing for variable '{var_definition.name}' at t={t}"
                    )
            merged[var_definition][t] = value

    def _collect_timed_forcings(self) -> dict[VariableDefinition, dict[float, Any]]:
        merged = {}

        for module in self.sim.modules:
            if module is self:
                continue

            hook_forcings = module.get_timedep_atmosphere_forcings()
            for var_definition, series in hook_forcings.items():
                self._merge_forcing_series(merged, var_definition, series)

            if not is_dataclass(module):
                continue

            # We also check for any TimeDependentScalar fields in the module dataclass, which we interpret as providing a time series for the corresponding variable.
            for dataclass_field in fields(module):
                value = getattr(module, dataclass_field.name)
                if not isinstance(value, TimeDependentScalar):
                    continue
                if dataclass_field.name in get_var_by_name():
                    var_definition = get_var_by_name()[dataclass_field.name]
                else:
                    raise NotImplementedError(
                        f"TimedependentModule cannot determine variable definition for field {dataclass_field.name}"
                    )
                self._merge_forcing_series(merged, var_definition, value.to_time_map())

        return merged

    def _to_profile_column(self, value, z_len: int) -> np.ndarray:
        if isinstance(value, Number):
            return np.asarray(float(value), dtype=float)
        arr = np.asarray(value, dtype=float)
        if arr.ndim == 0:
            return arr
        if arr.ndim == 1 and len(arr) == z_len:
            return arr
        raise ValueError(
            f"Timed forcing must be scalar or 1D profile of length {z_len}; got {arr.shape}"
        )

    def _build_timedep_array(
        self,
        var_definition: VariableDefinition,
        series: Mapping[float, Any],
        time_values: list[float],
        z_len: int,
    ) -> np.ndarray:
        columns = [self._to_profile_column(series[t], z_len) for t in time_values]
        first = columns[0]

        if np.asarray(first).ndim == 0:
            for value in columns:
                if np.asarray(value).ndim != 0:
                    raise ValueError(
                        f"Timed forcing for '{var_definition.name}' mixes scalar and profile values"
                    )
            return np.asarray([float(v) for v in columns], dtype=float)

        for value in columns:
            if np.asarray(value).ndim != 1:
                raise ValueError(
                    f"Timed forcing for '{var_definition.name}' mixes profile and scalar values"
                )
        return np.column_stack(columns)

    def __add__(self, obj) -> "TimedependentModule":
        if not isinstance(obj, FromLS2D):
            return NotImplemented
        self.usesLS2DforTime = obj
        return self

    def __iadd__(self, other):
        return self.__add__(other)

    def check_settings(self):
        return None

    def prepare_calculation(self):
        if self.usesLS2DforTime:
            if self.LS2Dmodule is None:
                self.LS2Dmodule = self.sim.retrieve_module(LS2DAtmosphereModule)
                if not self.LS2Dmodule.prepare_calculation_done:
                    self.LS2Dmodule.prepare_calculation()
                    self.LS2Dmodule.prepare_calculation_done = True
            self.timesteps = list(getattr(self.LS2Dmodule, "_times_with_zero", []))
        timesteps = [float(t) for t in self.timesteps]
        if len(timesteps) != len(set(timesteps)):
            raise ValueError("TimedependentModule.timesteps contains duplicates")
        runtime = self.retrieve_module(TimeModule).runtime
        if timesteps and timesteps[-1] < runtime:
            raise ValueError(
                f"TimedependentModule.timesteps must include the final time {runtime}"
            )
        if timesteps != sorted(timesteps):
            raise ValueError("TimedependentModule.timesteps must be in ascending order")
        self.timedep_var_dic = {}
        if self.grid is None:
            return None

        all_forcings = self._collect_timed_forcings()
        if not all_forcings:
            return None

        if not self.timesteps:
            raise ValueError(
                "TimedependentModule.timesteps must be set when time-dependent forcings are present"
            )

        time_values = sorted(float(t) for t in self.timesteps)
        z_len = len(self.grid.zt)

        for var_definition, series in all_forcings.items():
            if not var_definition.time_dependent_name:
                # Some timed series are init-only (e.g. nudging profiles in init.nc).
                continue
            series_times = sorted(float(t) for t in series.keys())
            if len(series_times) != len(time_values) or not np.allclose(
                np.asarray(series_times, dtype=float),
                np.asarray(time_values, dtype=float),
            ):
                raise ValueError(
                    f"Timed forcing times for '{var_definition.name}' {series_times} do not match TimedependentModule.timesteps {time_values}"
                )

            self.timedep_var_dic[var_definition] = self._build_timedep_array(
                var_definition,
                series,
                time_values,
                z_len,
            )

        # we need to set ntimedep to the number of steps, but we don't include 0.
        self.ntimedep = len(time_values) - 1
        if self.ntimedep < 2:
            raise ValueError(
                f"TimedependentModule requires at least 3 timesteps (including time=0) for interpolation; got {len(time_values)}"
            )
        return None

    def write_files(self):
        if not self.timedep_var_dic:
            return None
        self.profile_writer.write_time_dependent_profiles(
            self.grid.zt,
            [float(t) for t in self.timesteps],
            self.timedep_var_dic,
            self.output_path / "input",
            self.exp_id,
        )
        return None
