from typing import List, Optional, Union
from dataclasses import dataclass, field

from modular_dales.MODULE_REGISTRY import register_module
from modular_dales.modular.dales_simulation import dales_simulation
from modular_dales.modular.simulation_module import simulation_module


@register_module
@dataclass
class EasyOutputModule(simulation_module):
    """Output configuration module for DALES simulation.
    Outputs fielddumps, cape, cross-sections, and statistics at specified intervals.
    When added to a simulation, automatically updates the corresponding namelist parameters in the appropriate sections.

    Args:
        sim: Parent dales_simulation instance
        output_interval: Time interval(s) for output (can be a single value or list for different output types)
        enable_output: Whether to enable output (controls lfielddump, lcape, lcross, lstat in namelist)

    """

    sim: Optional["dales_simulation"] = field(default=None, repr=False)
    output_interval: Optional[Union[float, list[float]]] = field(
        default=60,
        metadata={
            "nml": [
                "namfielddump",
                "namcape",
                "namlsmcrosssection",
                "namgenstat",
                "namcrosssection",
                "namtimestat",
                "nambudget",
                "nambudget",
                "namgenstat",
                "namradfield",
                "namradfield",
            ],
            "key": [
                "dtav",
                "dtav",
                "dtav",
                "dtav",
                "dtav",
                "dtav",
                "dtav",
                "timeav",
                "timeav",
                "dtav",
                "timeav",
            ],
            "required": True,
        },
        init=True,
    )
    enable_output: bool = field(
        default=True,
        metadata={
            "nml": [
                "namfielddump",
                "namcape",
                "namlsmcrosssection",
                "namgenstat",
                "namcrosssection",
                "namtimestat",
                "nambudget",
                "namradfield",
            ],
            "key": [
                "lfielddump",
                "lcape",
                "lcross",
                "lstat",
                "lcross",
                "ltimestat",
                "lbudget",
                "lradfield",
            ],
            "required": True,
        },
    )

    def do_config(self):
        """Configure output-related namelist parameters."""
        return None

    def __post_init__(self):
        super().__init__(self.sim)
        self.module_name = "EasyOutputModule"

    def prepare_calculation(self):
        """No preparation needed."""
        return None

    def check_settings(self):
        """Check output settings validity."""
        return None

    def write_files(self):
        """Write output files if needed."""
        return None


@register_module
@dataclass
class IndependentOutputModule(simulation_module):
    """Output configuration module with independent settings per output type.

    This module is similar to EasyOutputModule, but every output type can be
    configured independently for both enable/disable flags and output timing.
    """

    sim: Optional["dales_simulation"] = field(default=None, repr=False)

    fielddump_enabled: bool = field(
        default=False,
        metadata={"nml": "namfielddump", "key": "lfielddump", "required": True},
    )
    fielddump_dtav: Optional[float] = field(
        default=60,
        metadata={"nml": "namfielddump", "key": "dtav", "required": True},
    )

    cape_enabled: bool = field(
        default=False,
        metadata={"nml": "namcape", "key": "lcape", "required": True},
    )
    cape_dtav: Optional[float] = field(
        default=60,
        metadata={"nml": "namcape", "key": "dtav", "required": True},
    )

    lsm_cross_enabled: bool = field(
        default=False,
        metadata={"nml": "namlsmcrosssection", "key": "lcross", "required": True},
    )
    lsm_cross_dtav: Optional[float] = field(
        default=60,
        metadata={"nml": "namlsmcrosssection", "key": "dtav", "required": True},
    )

    cross_enabled: bool = field(
        default=False,
        metadata={"nml": "namcrosssection", "key": "lcross", "required": True},
    )
    cross_dtav: Optional[float] = field(
        default=60,
        metadata={"nml": "namcrosssection", "key": "dtav", "required": True},
    )

    stats_enabled: bool = field(
        default=False,
        metadata={"nml": "namgenstat", "key": "lstat", "required": True},
    )
    stats_dtav: Optional[float] = field(
        default=60,
        metadata={"nml": "namgenstat", "key": "dtav", "required": True},
    )
    stats_timeav: Optional[float] = field(
        default=60,
        metadata={"nml": "namgenstat", "key": "timeav", "required": True},
    )

    timestat_enabled: bool = field(
        default=False,
        metadata={"nml": "namtimestat", "key": "ltimestat", "required": True},
    )
    timestat_dtav: Optional[float] = field(
        default=60,
        metadata={"nml": "namtimestat", "key": "dtav", "required": True},
    )

    budget_enabled: bool = field(
        default=False,
        metadata={"nml": "nambudget", "key": "lbudget", "required": True},
    )
    budget_dtav: Optional[float] = field(
        default=60,
        metadata={"nml": "nambudget", "key": "dtav", "required": True},
    )
    budget_timeav: Optional[float] = field(
        default=60,
        metadata={"nml": "nambudget", "key": "timeav", "required": True},
    )

    radfield_enabled: bool = field(
        default=False,
        metadata={"nml": "namradfield", "key": "lradfield", "required": True},
    )
    radfield_dtav: Optional[float] = field(
        default=60,
        metadata={"nml": "namradfield", "key": "dtav", "required": True},
    )
    radfield_timeav: Optional[float] = field(
        default=60,
        metadata={"nml": "namradfield", "key": "timeav", "required": True},
    )

    def do_config(self):
        """Configure output-related namelist parameters."""
        return None

    def __post_init__(self):
        super().__init__(self.sim)
        self.module_name = "IndependentOutputModule"

    def prepare_calculation(self):
        """No preparation needed."""
        return None

    def check_settings(self):
        """Check output settings validity."""
        return None

    def write_files(self):
        """Write output files if needed."""
        return None


@register_module
@dataclass
class CrossSectionOutputModule(simulation_module):
    """Output configuration module for cross-sections only.

    This module controls atmospheric cross-section output interval, selected
    planes, and their corresponding indices.
    """

    sim: Optional["dales_simulation"] = field(default=None, repr=False)

    cross_enabled: bool = field(
        default=True,
        metadata={"nml": "namcrosssection", "key": "lcross", "required": True},
    )
    cross_dtav: Optional[float] = field(
        default=60,
        metadata={"nml": "namcrosssection", "key": "dtav", "required": True},
    )
    xy: List[int] = field(
        default_factory=list,
        metadata={"nml": "namcrosssection", "key": "crossheight", "required": True},
    )
    xz: List[int] = field(
        default_factory=list,
        metadata={"nml": "namcrosssection", "key": "crossplane", "required": True},
    )
    yz: List[int] = field(
        default_factory=list,
        metadata={"nml": "namcrosssection", "key": "crossortho", "required": True},
    )
    xy_enabled: bool = field(
        default=False,
        metadata={"nml": "namcrosssection", "key": "lxy", "required": True},
    )
    xz_enabled: bool = field(
        default=False,
        metadata={"nml": "namcrosssection", "key": "lxz", "required": True},
    )
    yz_enabled: bool = field(
        default=False,
        metadata={"nml": "namcrosssection", "key": "lyz", "required": True},
    )

    def do_config(self):
        """Configure output-related namelist parameters."""
        return None

    def __post_init__(self):
        super().__init__(self.sim)
        self.module_name = "CrossSectionOutputModule"

    def prepare_calculation(self):
        if len(self.xy) > 0:
            self.xy_enabled = True
        else:
            self.xy_enabled = False
        if len(self.xz) > 0:
            self.xz_enabled = True
        else:
            self.xz_enabled = False
        if len(self.yz) > 0:
            self.yz_enabled = True
        else:
            self.yz_enabled = False
        """No preparation needed."""
        return None

    def check_settings(self):
        """Check output settings validity."""
        return None

    def write_files(self):
        """Write output files if needed."""
        return None


@register_module
@dataclass
class CheckSimulationModule(simulation_module):
    """
    Module for checking simulation validity and configuration.

    This module provides functionality to monitor and validate a DALES simulation
    during runtime. It allows output at an interval of `check_interval` seconds, at which DALES prints
    the currrent simulation time and ETA. It can also check for invalid values, verify tendencies, and
    optionally stop the simulation if invalid conditions are detected.

    Attributes:
        check_interval (Optional[int]): Time interval (in seconds) between output.
            Default is 60 seconds. Mapped to NAMCHECKSIM:tcheck in Fortran namelist.
            If 0, DALES prints every timestep.
        stop_on_invalid (bool): If True, stops the simulation when invalid values are detected.
            Default is False. Mapped to NAMCHECKSIM:lstop in Fortran namelist.
        check_tendencies (bool): If True, validates model tendencies during runtime checks.
            Default is False. Mapped to NAMCHECKSIM:lchecktend in Fortran namelist.
    """

    sim: Optional["dales_simulation"] = field(default=None, repr=False)
    check_interval: Optional[int] = field(
        default=60,
        metadata={
            "nml": "NAMCHECKSIM",
            "key": "tcheck",
            "serialize": True,
            "required": True,
        },
        init=True,
    )
    stop_on_invalid: bool = field(
        default=False,
        metadata={
            "nml": "NAMCHECKSIM",
            "key": "lstop",
            "serialize": True,
            "required": True,
        },
        init=True,
    )
    check_tendencies: bool = field(
        default=False,
        metadata={
            "nml": "NAMCHECKSIM",
            "key": "lchecktend",
            "serialize": True,
            "required": True,
        },
        init=True,
    )

    def __post_init__(self):
        super().__init__(self.sim)
        self.module_name = "CheckSimulationModule"

    def do_config(self):
        """No configuration needed."""
        return None

    def prepare_calculation(self):
        """Check simulation settings and configuration."""
        return None

    def check_settings(self):
        """No additional checks needed."""
        return None

    def write_files(self):
        """No files to write."""
        return None


@register_module
@dataclass
class SamplingModule(simulation_module):
    """Output configuration module for DALES simulation.
    Outputs samples at specified intervals.
    When added to a simulation, automatically updates the corresponding namelist parameters in the appropriate sections.

    Args:
        sim: Parent dales_simulation instance
        output_interval: Time interval(s) for output (can be a single value or list for different output types)
        enable_output: Whether to enable output (controls lsampcl, lsampco, lsampup, etc. in namelist)

    """

    sim: Optional["dales_simulation"] = field(default=None, repr=False)
    output_interval: Optional[Union[float, list[float]]] = field(
        default=60,
        metadata={
            "nml": ["namsampling", "namsampling"],
            "key": [
                "dtav",
                "timeav",
            ],
            "required": True,
        },
        init=True,
    )
    enable_output: bool = field(
        default=True,
        metadata={
            "nml": [
                "namsampling",
                "namsampling",
                "namsampling",
                "namsampling",
                "namsampling",
                "namsampling",
                "namsampling",
                "namsampling",
                "namsampling",
                "namsampling",
                "namsampling",
                "namsampling",
                "namsampling",
                "namsampling",
                "namsampling",
            ],
            "key": [
                "lsampcl",
                "lsampco",
                "lsampup",
                "lsampbuup",
                "lsampcldup",
                "lsamptend",
                "ltendleib",
                "lsamptendu",
                "lsamptendv",
                "lsamptendw",
                "lsamptendthl",
                "lsamptendqt",
                "lsamptendqr",
                "lsamptendnr",
                "lqlflux",
            ],
            "required": True,
        },
    )

    def do_config(self):
        """Configure output-related namelist parameters."""
        return None

    def __post_init__(self):
        super().__init__(self.sim)
        self.module_name = "SamplingModule"

    def prepare_calculation(self):
        """No preparation needed."""
        return None

    def check_settings(self):
        """Check output settings validity."""
        return None

    def write_files(self):
        """Write output files if needed."""
        return None
