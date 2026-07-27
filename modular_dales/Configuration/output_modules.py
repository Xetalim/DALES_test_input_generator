from typing import List, Optional, Union
from dataclasses import dataclass, field

from modular_dales.MODULE_REGISTRY import register_module
from modular_dales.modular.dales_simulation import dales_simulation
from modular_dales.modular.simulation_module import simulation_module


def _to_int_interval(value, field_name: str) -> int:
    """Normalize interval-like input to an integer number of seconds."""
    if value is None:
        raise ValueError(f"{field_name} cannot be None")
    try:
        ivalue = int(round(float(value)))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field_name} must be numeric, got {value!r}") from exc
    if ivalue <= 0:
        raise ValueError(f"{field_name} must be > 0, got {ivalue}")
    return ivalue


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
    output_interval: Optional[Union[int, list[int]]] = field(
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
        if isinstance(self.output_interval, list):
            self.output_interval = [
                _to_int_interval(v, "EasyOutputModule.output_interval")
                for v in self.output_interval
            ]
        else:
            self.output_interval = _to_int_interval(
                self.output_interval,
                "EasyOutputModule.output_interval",
            )

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
    fielddump_dtav: Optional[int] = field(
        default=60,
        metadata={"nml": "namfielddump", "key": "dtav", "required": True},
    )

    cape_enabled: bool = field(
        default=False,
        metadata={"nml": "namcape", "key": "lcape", "required": True},
    )
    cape_dtav: Optional[int] = field(
        default=60,
        metadata={"nml": "namcape", "key": "dtav", "required": True},
    )

    lsm_cross_enabled: bool = field(
        default=False,
        metadata={"nml": "namlsmcrosssection", "key": "lcross", "required": True},
    )
    lsm_cross_dtav: Optional[int] = field(
        default=60,
        metadata={"nml": "namlsmcrosssection", "key": "dtav", "required": True},
    )

    cross_enabled: bool = field(
        default=False,
        metadata={"nml": "namcrosssection", "key": "lcross", "required": True},
    )
    cross_dtav: Optional[int] = field(
        default=60,
        metadata={"nml": "namcrosssection", "key": "dtav", "required": True},
    )

    stats_enabled: bool = field(
        default=False,
        metadata={"nml": "namgenstat", "key": "lstat", "required": True},
    )
    stats_dtav: Optional[int] = field(
        default=60,
        metadata={"nml": "namgenstat", "key": "dtav", "required": True},
    )
    stats_timeav: Optional[int] = field(
        default=60,
        metadata={"nml": "namgenstat", "key": "timeav", "required": True},
    )

    timestat_enabled: bool = field(
        default=False,
        metadata={"nml": "namtimestat", "key": "ltimestat", "required": True},
    )
    timestat_dtav: Optional[int] = field(
        default=60,
        metadata={"nml": "namtimestat", "key": "dtav", "required": True},
    )

    budget_enabled: bool = field(
        default=False,
        metadata={"nml": "nambudget", "key": "lbudget", "required": True},
    )
    budget_dtav: Optional[int] = field(
        default=60,
        metadata={"nml": "nambudget", "key": "dtav", "required": True},
    )
    budget_timeav: Optional[int] = field(
        default=60,
        metadata={"nml": "nambudget", "key": "timeav", "required": True},
    )

    radfield_enabled: bool = field(
        default=False,
        metadata={"nml": "namradfield", "key": "lradfield", "required": True},
    )
    radfield_dtav: Optional[int] = field(
        default=60,
        metadata={"nml": "namradfield", "key": "dtav", "required": True},
    )
    radfield_timeav: Optional[int] = field(
        default=60,
        metadata={"nml": "namradfield", "key": "timeav", "required": True},
    )

    def do_config(self):
        """Configure output-related namelist parameters."""
        return None

    def __post_init__(self):
        super().__init__(self.sim)
        self.module_name = "IndependentOutputModule"
        for name in (
            "fielddump_dtav",
            "cape_dtav",
            "lsm_cross_dtav",
            "cross_dtav",
            "stats_dtav",
            "stats_timeav",
            "timestat_dtav",
            "budget_dtav",
            "budget_timeav",
            "radfield_dtav",
            "radfield_timeav",
        ):
            setattr(self, name, _to_int_interval(getattr(self, name), name))

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
    cross_dtav: Optional[int] = field(
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
        self.cross_dtav = _to_int_interval(
            self.cross_dtav,
            "CrossSectionOutputModule.cross_dtav",
        )

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
    output_interval: Optional[Union[int, list[int]]] = field(
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
        if isinstance(self.output_interval, list):
            self.output_interval = [
                _to_int_interval(v, "SamplingModule.output_interval")
                for v in self.output_interval
            ]
        else:
            self.output_interval = _to_int_interval(
                self.output_interval,
                "SamplingModule.output_interval",
            )

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
class FielddumpModule(simulation_module):
    """Dedicated output module for NAMFIELDDUMP with individual variable switches."""

    sim: Optional["dales_simulation"] = field(default=None, repr=False)

    dtav: int = field(default=60, metadata={"nml": "namfielddump", "key": "dtav", "required": True, "doc": "Sampling interval in seconds for field dumps."})
    lfielddump: bool = field(default=True, metadata={"nml": "namfielddump", "key": "lfielddump", "required": True, "doc": "Enable field dump output."})
    ldiracc: bool = field(default=False, metadata={"nml": "namfielddump", "key": "ldiracc", "doc": "Use direct-access output file mode."})
    lbinary: bool = field(default=False, metadata={"nml": "namfielddump", "key": "lbinary", "doc": "Write binary fielddump files instead of text."})
    klow: int = field(default=1, metadata={"nml": "namfielddump", "key": "klow", "doc": "Lowest model level index included in dump."})
    khigh: int = field(default=0, metadata={"nml": "namfielddump", "key": "khigh", "doc": "Highest model level index included in dump; 0 lets DALES decide."})
    ncoarse: int = field(default=1, metadata={"nml": "namfielddump", "key": "ncoarse", "doc": "Horizontal coarsening factor for dumped fields."})
    tmin: float = field(default=0.0, metadata={"nml": "namfielddump", "key": "tmin", "doc": "Start time in seconds for writing dumps."})
    tmax: float = field(default=1.0e30, metadata={"nml": "namfielddump", "key": "tmax", "doc": "End time in seconds for writing dumps."})

    lu: bool = field(default=False, metadata={"nml": "namfielddump", "key": "lu", "doc": "Dump u velocity field."})
    lv: bool = field(default=False, metadata={"nml": "namfielddump", "key": "lv", "doc": "Dump v velocity field."})
    lw: bool = field(default=False, metadata={"nml": "namfielddump", "key": "lw", "doc": "Dump w velocity field."})
    lqt: bool = field(default=False, metadata={"nml": "namfielddump", "key": "lqt", "doc": "Dump total water mixing ratio field."})
    lql: bool = field(default=False, metadata={"nml": "namfielddump", "key": "lql", "doc": "Dump liquid water mixing ratio field."})
    lthl: bool = field(default=False, metadata={"nml": "namfielddump", "key": "lthl", "doc": "Dump liquid water potential temperature field."})
    lbuoy: bool = field(default=False, metadata={"nml": "namfielddump", "key": "lbuoy", "doc": "Dump buoyancy field."})
    lcli: bool = field(default=False, metadata={"nml": "namfielddump", "key": "lcli", "doc": "Dump cloud ice field."})
    lclw: bool = field(default=False, metadata={"nml": "namfielddump", "key": "lclw", "doc": "Dump cloud liquid water field."})
    lta: bool = field(default=False, metadata={"nml": "namfielddump", "key": "lta", "doc": "Dump absolute temperature field."})
    lplw: bool = field(default=False, metadata={"nml": "namfielddump", "key": "lplw", "doc": "Dump longwave pressure-related field."})
    lpli: bool = field(default=False, metadata={"nml": "namfielddump", "key": "lpli", "doc": "Dump ice-phase pressure-related field."})
    lhus: bool = field(default=False, metadata={"nml": "namfielddump", "key": "lhus", "doc": "Dump specific humidity field."})
    lhur: bool = field(default=False, metadata={"nml": "namfielddump", "key": "lhur", "doc": "Dump relative humidity field."})
    ltntr: bool = field(default=False, metadata={"nml": "namfielddump", "key": "ltntr", "doc": "Dump total tendency field (resolved)."})
    ltntrs: bool = field(default=False, metadata={"nml": "namfielddump", "key": "ltntrs", "doc": "Dump subgrid tendency field."})
    ltntrl: bool = field(default=False, metadata={"nml": "namfielddump", "key": "ltntrl", "doc": "Dump large-scale tendency field."})
    le12: bool = field(default=False, metadata={"nml": "namfielddump", "key": "le12", "doc": "Dump subgrid TKE field."})
    lekh: bool = field(default=False, metadata={"nml": "namfielddump", "key": "lekh", "doc": "Dump heat diffusivity field."})
    lekm: bool = field(default=False, metadata={"nml": "namfielddump", "key": "lekm", "doc": "Dump momentum diffusivity field."})
    lsv: bool = field(default=False, metadata={"nml": "namfielddump", "key": "lsv", "doc": "Dump scalar tracer fields."})

    def __post_init__(self):
        super().__init__(self.sim)
        self.module_name = "FielddumpModule"
        self.dtav = _to_int_interval(self.dtav, "FielddumpModule.dtav")

    def do_config(self):
        return None

    def prepare_calculation(self):
        return None

    def check_settings(self):
        if self.khigh != 0 and self.khigh < self.klow:
            raise ValueError("FielddumpModule: khigh must be >= klow, or 0 for DALES default")
        return None

    def write_files(self):
        return None
