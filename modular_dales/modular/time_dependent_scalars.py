from dataclasses import dataclass, field

import numpy as np

from modular_dales.MODULE_REGISTRY import MODULE_REGISTRY, register_special_serializing


@register_special_serializing
@dataclass
class TimeDependentScalar:
    """Time-dependent scalar value defined on explicit timesteps."""

    times: list[float] = field(
        default_factory=list, metadata={"serialize": True}, init=True
    )
    values: list[float] = field(
        default_factory=list, metadata={"serialize": True}, init=True
    )

    def __post_init__(self):
        if isinstance(self.times, np.ndarray):
            self.times = self.times.tolist()
        if isinstance(self.values, np.ndarray):
            self.values = self.values.tolist()
        if len(self.times) != len(self.values):
            raise ValueError(
                f"TimeDependentScalar times(len={len(self.times)}) and values(len={len(self.values)}) must have equal length"
            )
        if len(self.times) == 0:
            raise ValueError("TimeDependentScalar requires at least one point")
        normalized = [float(t) for t in self.times]
        if len(set(normalized)) != len(normalized):
            raise ValueError("TimeDependentScalar times contain duplicates")
        if all(t != 0 for t in normalized):
            raise ValueError("TimeDependentScalar times must include 0")

    def to_time_map(self) -> dict[float, float]:
        return {float(t): float(v) for t, v in zip(self.times, self.values)}

    def value_at_start(self) -> float:
        time_map = self.to_time_map()
        if 0.0 in time_map:
            return time_map[0.0]
        first_time = min(time_map.keys())
        return time_map[first_time]


# Register TimeDependentScalar for YAML deserialization without making it a simulation module
MODULE_REGISTRY["TimeDependentScalar"] = TimeDependentScalar
