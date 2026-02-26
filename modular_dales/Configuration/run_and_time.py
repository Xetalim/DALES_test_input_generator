import logging
from dataclasses import dataclass, field
from typing import Optional

from modular_dales.MODULE_REGISTRY import register_module
from modular_dales.modular.dales_simulation import dales_simulation
from modular_dales.modular.simulation_module import simulation_module

logger = logging.getLogger(__name__)


@register_module
@dataclass
class RunModule(simulation_module):
    """Run configuration module for DALES simulation.

    Handles all run-related namelist parameters in the RUN section.

    Args:
        sim: Parent dales_simulation instance
        iexpnr: Experiment number
        runtime: Total simulation time in seconds
        ladaptive: Allow adaptive timestepping (highly recommended)
        dtmax: Maximum timestep
    """

    sim: Optional["dales_simulation"] = field(default=None, repr=False)
    iexpnr: Optional[int] = field(
        default=None, metadata={"nml": "RUN", "key": "iexpnr"}
    )
    runtime: Optional[float] = field(
        default=None, metadata={"nml": "RUN", "key": "runtime"}
    )
    ladaptive: Optional[bool] = field(
        default=None, metadata={"nml": "RUN", "key": "ladaptive"}
    )
    dtmax: Optional[float] = field(
        default=None, metadata={"nml": "RUN", "key": "dtmax"}
    )

    def __post_init__(self):
        super().__init__(self.sim)
        self.module_name = "RunModule"

    def prepare_calculation(self):
        """Set up run-related namelist parameters."""
        return None

    def check_settings(self):
        """Check run settings validity."""
        return None

    def write_files(self):
        """No files to write for run module."""
        return None


@register_module
@dataclass
class TimeModule(simulation_module):
    """Simulation time configuration module.

    Handles time-related parameters in the RUN and DOMAIN sections.

    Args:
        sim: Parent dales_simulation instance
        xtime: Time within the day in seconds (0-86400)
        xday: Day of year (1-365/366)
        xyear: Year
        runtime: Total runtime in seconds
        trestart: Restart file write interval in seconds
    """

    sim: Optional["dales_simulation"] = field(default=None, repr=False)
    xtime: float = field(
        default=0.0, metadata={"nml": "DOMAIN", "key": "xtime", "required": True}
    )
    xday: int = field(
        default=265, metadata={"nml": "DOMAIN", "key": "xday", "required": True}
    )
    xyear: int = field(
        default=2025, metadata={"nml": "DOMAIN", "key": "xyear", "required": True}
    )
    enable_datetime: bool = field(
        default=True,
        metadata={"nml": "namdatetime", "key": "l_datetime", "required": True},
    )
    timezone: int = field(default=0, metadata={"nml": "namdatetime", "key": "timezone"})
    startyear: int = field(
        default=2025,
        metadata={"nml": "namdatetime", "key": "startyear", "required": True},
    )
    startmonth: int = field(
        default=1,
        metadata={"nml": "namdatetime", "key": "startmonth", "required": True},
    )
    startday: int = field(
        default=1, metadata={"nml": "namdatetime", "key": "startday", "required": True}
    )

    runtime: Optional[float] = field(
        default=None, metadata={"nml": "RUN", "key": "runtime", "required": True}
    )
    trestart: Optional[float] = field(
        default=None, metadata={"nml": "RUN", "key": "trestart"}
    )

    def __post_init__(self):
        super().__init__(self.sim)
        self.module_name = "TimeModule"

    def do_config(self):
        """Configure time-related namelist parameters."""
        return None

    def check_settings(self):
        """Check time settings validity."""
        if self.xyear != self.startyear:
            logger.warning(
                "TimeModule: xyear (%d) does not match startyear (%d)",
                self.xyear,
                self.startyear,
            )
        return None

    def prepare_calculation(self):
        """No calculation work needed."""
        return None

    def write_files(self):
        """No files to write."""
        return None
