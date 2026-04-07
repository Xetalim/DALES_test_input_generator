import logging

from dataclasses import dataclass, field
from typing import Union, Optional, List, Mapping, Any

import numpy as np
from scipy.interpolate import interp1d

from modular_dales.modular.simulation_module import simulation_module
from modular_dales.modular.time_dependent_scalars import TimeDependentScalar
from modular_dales.MODULE_REGISTRY import register_module
from modular_dales.vars import VariableDefinition, get_all_vars

from modular_dales.Atmosphere.shapes import SHAPE_FUNCTIONS
from modular_dales.IO_helpers import AtmosphereProfileWriter

logger = logging.getLogger(__name__)
logger.debug("Entered module: %s", __name__)

NUDGING_VAR_NAMES = {
    "ua_nudge",
    "va_nudge",
    "thl_nudge",
    "wa_nudge",
    "qt_nudge",
}


@dataclass
class AtmosphereVariable:
    """
    Represents an atmospheric variable with its definition and associated values.

    Attributes:
        definition: The variable definition containing metadata and specifications.
        is_time_dependent: Whether the variable varies over time. Defaults to False.
        values: Optional numpy array containing the variable's values.
    """

    definition: VariableDefinition
    is_time_dependent: bool = False
    values: Optional[np.ndarray] = field(
        default=None,
        repr=False,
        metadata={"serialize": False},
    )


def _build_default_variables() -> dict[VariableDefinition, AtmosphereVariable]:
    """Create default AtmosphereVariable instances for all known defs.

    The dict is keyed by the variable definition.
    """

    return build_default_variables(get_all_vars())


def build_default_variables(
    ALL_VARIABLES,
) -> dict[VariableDefinition, AtmosphereVariable]:
    return {
        var_def: AtmosphereVariable(definition=var_def)
        for var_def in ALL_VARIABLES
        if (var_def.is_profile and not var_def.must_only_be_time_dependent)
    }


def evaluate_profile_map(
    z: np.ndarray,
    profiles_by_var: Mapping[
        VariableDefinition, Union["AtmosphericProfile", "InterpolatedProfile"]
    ],
):
    var_dic = {
        var_def: None
        for var_def in get_all_vars()
        if var_def.is_profile and not var_def.must_only_be_time_dependent
    }
    for var_def, profile in profiles_by_var.items():
        var_dic[var_def] = profile.evaluate(z)
    return var_dic


@register_module
@dataclass
class AtmosphericProfile:
    """Single atmospheric profile configuration."""

    variable: VariableDefinition
    """Variable definition (see modular_dales.Atmosphere.vars)."""
    shape: str
    """Shape function (lin, exp, linmlsurf, expsinw, etc.)"""
    params: dict[str, Union[float, TimeDependentScalar]]
    """Shape-specific parameters"""

    def _resolve_params(self, time_value: Optional[float] = None) -> dict:

        resolved = {}
        for key, value in self.params.items():
            if isinstance(value, TimeDependentScalar):
                if time_value is None:
                    resolved[key] = value.value_at_start()
                else:
                    time_map = value.to_time_map()
                    t = float(time_value)
                    if t not in time_map:
                        raise ValueError(
                            f"TimeDependentScalar parameter '{key}' for '{self.variable.name}' has no value at t={t}"
                        )
                    resolved[key] = time_map[t]
            else:
                resolved[key] = value
        return resolved

    def _time_dependent_param_times(self) -> Optional[List[float]]:

        timed_params = [
            value
            for value in self.params.values()
            if isinstance(value, TimeDependentScalar)
        ]
        if not timed_params:
            return None

        reference = sorted(float(t) for t in timed_params[0].times)
        for timed_param in timed_params[1:]:
            current = sorted(float(t) for t in timed_param.times)
            if len(current) != len(reference) or not np.allclose(
                np.asarray(current, dtype=float),
                np.asarray(reference, dtype=float),
            ):
                raise ValueError(
                    f"All TimeDependentScalar params in '{self.variable.name}' must use identical time points"
                )
        return reference

    def evaluate(self, z: np.ndarray) -> np.ndarray:
        if self.shape not in SHAPE_FUNCTIONS:
            raise ValueError(
                f"Unknown shape function '{self.shape}' for variable '{self.variable.name}'"
            )
        return SHAPE_FUNCTIONS[self.shape](z, **self._resolve_params())

    def timed_forcing_series(self, z: np.ndarray) -> Optional[dict[float, np.ndarray]]:
        times = self._time_dependent_param_times()
        if not times:
            return None
        return {
            float(time_value): SHAPE_FUNCTIONS[self.shape](
                z,
                **self._resolve_params(time_value),
            )
            for time_value in times
        }


@register_module
@dataclass
class InterpolatedProfile:
    """Single atmospheric profile configuration."""

    variable: VariableDefinition
    """Variable definition (see modular_dales.Atmosphere.vars)."""
    z: list[float]
    """list of z coordinates for interpolation"""
    points: list[Union[float, TimeDependentScalar]]
    """Values at the specified z coordinates"""
    fill_value: Union[float, str] = "extrapolate"
    """Value to use for beyond the provided z range, or 'extrapolate' to extrapolate."""

    def _extract_points_at_time(self, time_value: Optional[float] = None) -> list:

        extracted = []
        for v in self.points:
            if isinstance(v, TimeDependentScalar):
                if time_value is None:
                    extracted.append(v.value_at_start())
                else:
                    time_map = v.to_time_map()
                    t = float(time_value)
                    if t not in time_map:
                        raise ValueError(
                            f"TimeDependentScalar point for '{self.variable.name}' has no value at t={t}"
                        )
                    extracted.append(time_map[t])
            else:
                extracted.append(v)
        return extracted

    def _time_dependent_point_times(self) -> Optional[List[float]]:

        timed_points = [
            point for point in self.points if isinstance(point, TimeDependentScalar)
        ]
        if not timed_points:
            return None

        reference = sorted(float(t) for t in timed_points[0].times)
        for timed_point in timed_points[1:]:
            current = sorted(float(t) for t in timed_point.times)
            if len(current) != len(reference) or not np.allclose(
                np.asarray(current, dtype=float),
                np.asarray(reference, dtype=float),
            ):
                raise ValueError(
                    f"All TimeDependentScalar points in '{self.variable.name}' must use identical time points"
                )
        return reference

    def evaluate(self, z: np.ndarray) -> np.ndarray:
        # Extract initial values from any TimeDependentScalar points
        points = self._extract_points_at_time()
        return interp1d(
            self.z,
            points,
            fill_value=self.fill_value,
        )(z)

    def timed_forcing_series(self, z: np.ndarray) -> Optional[dict[float, np.ndarray]]:
        times = self._time_dependent_point_times()
        if not times:
            return None
        return {
            float(time_value): interp1d(
                self.z,
                self._extract_points_at_time(time_value),
                fill_value=self.fill_value,
            )(z)
            for time_value in times
        }


@register_module
@dataclass
class TimedAtmosphereProfile:
    """Profile assignment at a specific simulation time (seconds)."""

    time: float
    profile: Union[AtmosphericProfile, InterpolatedProfile]


@register_module
@dataclass
class AtmosphereModule(simulation_module):
    """Atmosphere module for profile setup.

    Args:
        sim: Parent simulation instance
        shaped_profiles: Mapping of variable name to shaped profile configuration
        interpolated_profiles: Mapping of variable name to interpolated profile configuration
    """

    sim: Optional["simulation_module"] = field(default=None, repr=False)
    shaped_profiles: List[AtmosphericProfile] = field(default_factory=list)
    interpolated_profiles: List[InterpolatedProfile] = field(default_factory=list)
    timed_profiles: List[TimedAtmosphereProfile] = field(default_factory=list)
    variables: dict[VariableDefinition, AtmosphereVariable] = field(
        default_factory=_build_default_variables,
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

    collected_base_profiles: dict[
        VariableDefinition, Union[AtmosphericProfile, InterpolatedProfile]
    ] = field(
        default_factory=dict,
        init=False,
        repr=False,
        metadata={"serialize": False},
    )
    collected_timed_profiles_by_time: dict[
        float, dict[VariableDefinition, Union[AtmosphericProfile, InterpolatedProfile]]
    ] = field(
        default_factory=dict,
        init=False,
        repr=False,
        metadata={"serialize": False},
    )
    extra_init_time_height_fields: dict[str, dict[str, Any]] = field(
        default_factory=dict,
        init=False,
        repr=False,
        metadata={"serialize": False},
    )
    extra_init_times: Optional[List[float]] = field(
        default=None,
        init=False,
        repr=False,
        metadata={"serialize": False},
    )

    def __post_init__(self):
        super().__init__(self.sim)
        self.module_name = "AtmosphereModule"

    def _initialize_from_sim(self, sim):
        sim.has_atmosphere_module = True
        return super()._initialize_from_sim(sim)

    # Generic add method
    def __add__(
        self,
        other: Union[
            AtmosphericProfile,
            InterpolatedProfile,
            TimedAtmosphereProfile,
            list,
            tuple,
        ],
    ):
        if isinstance(other, AtmosphericProfile):
            self.shaped_profiles.append(other)
        elif isinstance(other, InterpolatedProfile):
            self.interpolated_profiles.append(other)
        elif isinstance(other, TimedAtmosphereProfile):
            self.timed_profiles.append(other)
        elif (
            isinstance(other, tuple)
            and len(other) == 2
            and isinstance(other[0], (int, float))
            and isinstance(other[1], (AtmosphericProfile, InterpolatedProfile))
        ):
            self.timed_profiles.append(
                TimedAtmosphereProfile(time=float(other[0]), profile=other[1])
            )
        elif isinstance(other, (list, tuple)):
            for item in other:
                self += item  # recursive
        else:
            raise TypeError(
                f"Cannot add object of type {type(other)} to AtmosphereModule"
            )
        return self

    def __iadd__(self, other):
        return self.__add__(other)

    def _validate_profile_variable(
        self, profile: Union[AtmosphericProfile, InterpolatedProfile]
    ):
        if profile.variable not in self.variables:
            raise ValueError(f"Unknown atmospheric variable '{profile.variable.name}'")

    def _collect_base_profiles(
        self,
    ) -> dict[VariableDefinition, Union[AtmosphericProfile, InterpolatedProfile]]:
        self.collected_base_profiles = {}
        for profile in self.shaped_profiles:
            self._validate_profile_variable(profile)
            self.collected_base_profiles[profile.variable] = profile
        for profile in self.interpolated_profiles:
            self._validate_profile_variable(profile)
            self.collected_base_profiles[profile.variable] = profile
        return self.collected_base_profiles

    def _collect_timed_profiles(
        self,
    ) -> dict[
        float, dict[VariableDefinition, Union[AtmosphericProfile, InterpolatedProfile]]
    ]:
        self.collected_timed_profiles_by_time = {}
        for timed_profile in self.timed_profiles:
            time = float(timed_profile.time)
            profile = timed_profile.profile
            self._validate_profile_variable(profile)
            var_definition = profile.variable
            if not var_definition.can_be_time_dependent:
                raise ValueError(
                    f"Variable '{var_definition.name}' cannot be time-dependent"
                )
            if time not in self.collected_timed_profiles_by_time:
                self.collected_timed_profiles_by_time[time] = {}
            self.collected_timed_profiles_by_time[time][var_definition] = profile
        return self.collected_timed_profiles_by_time

    def _build_timedep_profile_forcings(
        self,
        base_profiles: Mapping[
            VariableDefinition, Union[AtmosphericProfile, InterpolatedProfile]
        ],
        timed_profiles: Mapping[
            float,
            Mapping[VariableDefinition, Union[AtmosphericProfile, InterpolatedProfile]],
        ],
    ) -> dict[VariableDefinition, dict[float, np.ndarray]]:
        if not timed_profiles:
            return {}

        series_by_var = {}
        all_times = sorted(time for time in timed_profiles.keys())
        if all(time != 0 for time in all_times):
            raise ValueError("Timed profiles must include a profile at time 0")

        for time_value in all_times:
            merged_profiles = dict(base_profiles)
            merged_profiles.update(timed_profiles[time_value])
            profiles_at_time = evaluate_profile_map(self.grid.zt, merged_profiles)

            for var_definition in timed_profiles[time_value]:
                if var_definition not in series_by_var:
                    series_by_var[var_definition] = {}
                series_by_var[var_definition][time_value] = profiles_at_time[
                    var_definition
                ]

        return series_by_var

    def _build_timedep_param_forcings(
        self,
        base_profiles: Mapping[
            VariableDefinition, Union[AtmosphericProfile, InterpolatedProfile]
        ],
    ) -> dict[VariableDefinition, dict[float, np.ndarray]]:
        series_by_var: dict[VariableDefinition, dict[float, np.ndarray]] = {}
        for var_definition, profile in base_profiles.items():
            # Handle both AtmosphericProfile and InterpolatedProfile with time-dependent content
            timed_series = None
            if isinstance(profile, AtmosphericProfile):
                timed_series = profile.timed_forcing_series(self.grid.zt)
            elif isinstance(profile, InterpolatedProfile):
                timed_series = profile.timed_forcing_series(self.grid.zt)

            if not timed_series:
                continue

            if not var_definition.can_be_time_dependent:
                raise ValueError(
                    f"Variable '{var_definition.name}' cannot be time-dependent"
                )
            series_by_var[var_definition] = timed_series
        return series_by_var

    def prepare_calculation(self):
        """Prepare atmosphere profiles and update internal state."""
        if self.grid is None:
            raise ValueError("AtmosphereModule requires a simulation grid")

        base_profiles = self._collect_base_profiles()
        timed_profiles = self._collect_timed_profiles()
        timedep_param_forcings = self._build_timedep_param_forcings(base_profiles)

        built_profiles = evaluate_profile_map(self.grid.zt, base_profiles)
        timedep_varnames = set()
        for timedep_var in timed_profiles.values():
            timedep_varnames.update(timedep_var.keys())
        timedep_varnames.update(timedep_param_forcings.keys())

        for var_definition, variable in self.variables.items():
            variable.values = built_profiles[var_definition]
            variable.is_time_dependent = var_definition in timedep_varnames

        logger.info("AtmosphereModule: Profiles prepared")

    def check_settings(self):
        """Check atmosphere settings validity."""
        return None

    def write_files(self):
        """Write atmosphere profile files and generate plots."""
        if self.grid is None:
            return

        init_times = None
        timedep_forcings = self.get_timedep_atmosphere_forcings()

        nudging_series = {
            var_definition: series
            for var_definition, series in timedep_forcings.items()
            if var_definition.name in NUDGING_VAR_NAMES
        }
        configured_nudging_vars = {
            var_definition
            for var_definition in self.collected_base_profiles.keys()
            if var_definition.name in NUDGING_VAR_NAMES
        }
        if nudging_series:
            all_time_sets = {
                tuple(sorted(float(t) for t in series.keys()))
                for series in nudging_series.values()
            }
            if len(all_time_sets) != 1:
                raise ValueError(
                    "All nudging variables must use identical time points in AtmosphereModule"
                )
            init_times = list(next(iter(all_time_sets)))
            if 0.0 not in init_times:
                raise ValueError(
                    "Nudging variables must include time 0.0 in AtmosphereModule"
                )
        elif self.collected_timed_profiles_by_time:
            init_times = sorted(float(t) for t in self.collected_timed_profiles_by_time)

        # If a nudging variable is configured only as a base profile, treat it
        # as constant in time over the active init time axis.
        missing_timed_nudging = [
            var_definition
            for var_definition in configured_nudging_vars
            if var_definition not in nudging_series
        ]
        if missing_timed_nudging:
            if init_times is None:
                init_times = [0.0]
            for var_definition in missing_timed_nudging:
                base_values = self.variables[var_definition].values
                if base_values is None:
                    raise ValueError(
                        f"Nudging variable '{var_definition.name}' has no evaluated base profile"
                    )
                nudging_series[var_definition] = {
                    float(t): np.asarray(base_values, dtype=float) for t in init_times
                }

        if self.extra_init_time_height_fields:
            if self.extra_init_times is None:
                raise ValueError(
                    "AtmosphereModule has extra init fields but no extra_init_times"
                )
            extra_times = [float(t) for t in self.extra_init_times]
            if init_times is None:
                init_times = extra_times
            elif len(extra_times) != len(init_times) or not np.allclose(
                np.asarray(extra_times, dtype=float),
                np.asarray(init_times, dtype=float),
            ):
                raise ValueError(
                    "AtmosphereModule extra init field times do not match init_times"
                )

        output_input_path = self.output_path / "input"
        self.profile_writer.write_base_profiles(
            self.grid.zt,
            self.variables,
            output_input_path,
            self.exp_id,
            times=init_times,
            time_profile_series=nudging_series if nudging_series else None,
            extra_time_height_fields=self.extra_init_time_height_fields,
        )

        self.profile_writer.plot_profiles(
            self.grid.zt,
            self.variables,
            output_input_path,
            self.exp_id,
        )

    def get_timedep_atmosphere_forcings(
        self,
    ) -> dict[VariableDefinition, dict[float, np.ndarray]]:
        """
        :meta private:
        """
        base_profiles = self._collect_base_profiles()
        timed_profiles = self._collect_timed_profiles()
        timedep_forcings = self._build_timedep_profile_forcings(
            base_profiles, timed_profiles
        )
        timedep_param_forcings = self._build_timedep_param_forcings(base_profiles)

        for var_definition, series in timedep_param_forcings.items():
            if var_definition in timedep_forcings:
                for time_value, values in series.items():
                    if time_value in timedep_forcings[
                        var_definition
                    ] and not np.allclose(
                        timedep_forcings[var_definition][time_value], values
                    ):
                        raise ValueError(
                            f"Conflicting timed atmosphere forcing for '{var_definition.name}' at t={time_value}"
                        )
                    timedep_forcings[var_definition][time_value] = values
            else:
                timedep_forcings[var_definition] = series

        return timedep_forcings
