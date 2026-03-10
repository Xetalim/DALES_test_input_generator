"""Open Boundary Conditions module for lateral boundary forcing."""

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np
import xarray as xr
import dask

from modular_dales.Geometry import GridDalesOpenBC
from modular_dales.modular import simulation_module
from modular_dales.MODULE_REGISTRY import register_module
from modular_dales.Atmosphere import AtmosphereModule

from modular_dales.LBC.nest_dales_in_dales import (
    boundary_fields_fine,
    initial_fields_fine,
)
from modular_dales.LBC.nest_dales_in_HARMONIE import (
    boundary as harmonie_boundary,
    initfields,
    prep_harmonie,
)

from modular_dales.vars import get_all_vars, get_var_by_name

logger = logging.getLogger(__name__)
logger.debug("Entered module: %s", __name__)
### FOR DEBUGGING, SET THE DASK SCHEDULER TO SINGLE-THREADED
### ELSE LEAVE COMMENTED OUT
# dask.config.set(
#     scheduler="single-threaded"
# )  # overwrite default with single-threaded scheduler

# if you're on a large cluster, you can for example set up a dask SLURM cluster which has a lot of memory so you can do any open BC calculations better.
# for example, here is how you can setup a SLURM cluster at ATOS. (set your own account), with this you can set tchunk=50 in OPENBC.
# from dask_jobqueue import SLURMCluster
# from dask.distributed import Client
# cluster = SLURMCluster(
#     cores=16,
#     processes=8,
#     memory="128GB",
#     account=ACCOUNT,
#     walltime="00:20:00",
#     job_script_prologue=[
#         'export LANG="en_US.utf8"',
#         'export LANGUAGE="en_US.utf8"',
#         'export LC_ALL="en_US.utf8"',
#     ],
#     job_extra_directives=['--qos="nf"'],
# )
# client = Client(cluster)


@dataclass
class Nest_in_Harmonie:
    """Nest a DALES simulation inside a HARMONIE model environment.

    This class configures and manages the coupling between a DALES
    simulation and an enclosing HARMONIE model, handling the exchange
    of necessary meteorological and surface fields.

    Args:
        ml_glob : Optional[str]
            Glob pattern or path template for HARMONIE model-level fields
            (e.g., 3D atmospheric state) used to drive or constrain DALES.
        sfc_glob : Optional[str]
            Glob pattern or path template for HARMONIE surface fields
            (e.g., surface fluxes, skin temperature) required for the nesting.
    """

    ml_glob: Optional[str] = field(
        default=None, repr=True, metadata={"serialize": True}, init=True
    )
    sfc_glob: Optional[str] = field(
        default=None, repr=True, metadata={"serialize": True}, init=True
    )


@dataclass
class Nest_in_Dales:
    """Nest DALES inside another DALES simulation

    Args:
        outpath_coarse : Optional[str]
            Path to the output directory of the coarse parent DALES simulation
            used to source boundary fields for all time steps except t=0.
        outpath_coarse_old : Optional[str]
            Path to the output directory of the previous DALES simulation,
            used to source initial boundary fields at t=0
            from the last time step in the output of the previous simulation
        inpath_coarse : Optional[str]
            Path to the input directory of the previous DALES simulation,
            used to source initial boundary fields at t=0
            from the initfields.nc input of the previous simulation
        inpath : Optional[str]
            Path to the input directory of the current DALES simulation,
            used to source initial fields for the current simulation
    """

    outpath_coarse: Optional[str] = field(
        default=None, repr=True, metadata={"serialize": True}, init=True
    )
    outpath_coarse_old: Optional[str] = field(
        default=None, repr=True, metadata={"serialize": True}, init=True
    )
    inpath_coarse: Optional[str] = field(
        default=None, repr=True, metadata={"serialize": True}, init=True
    )
    inpath: Optional[str] = field(
        default=None, repr=True, metadata={"serialize": True}, init=True
    )


@dataclass
class Nest_in_AtmosphereProfiles:
    """Nest DALES in a horizontally homogeneous atmosphere from profiles.

    This mode does not require a separate coarse DALES or HARMONIE run.
    Instead, it reuses the vertical profiles configured in the
    :class:`AtmosphereModule` of the *same* simulation to construct
    open boundary fields that are uniform in x and y.

    By default, the following mapping is used between open-boundary
    variables and atmospheric profile variables (see ``vars.ATMO_VARS_BY_NAME``):

        - u   <- ug
        - v   <- vg
        - w   <- w
        - thl <- thetal
        - qt  <- qt
        - e12 <- tke

    You can override this mapping via ``variable_mapping`` if needed.

        This mode does *not* use any existing ``AtmosphereModule`` that is
        part of the :class:`dales_simulation` modules. Instead, you must
        provide an explicit :class:`AtmosphereModule` instance via
        ``atmosphere_module``. This instance does not need to be added to
        the simulation's module list, but it must have a valid ``sim`` with
        a grid attached.
    """

    variable_mapping: Dict[str, str] = field(
        default_factory=lambda: {
            "u": "ug",
            "v": "vg",
            "w": "w",
            "thl": "thetal",
            "qt": "qt",
            "e12": "tke",
        },
        repr=True,
        metadata={"serialize": True},
        init=True,
    )
    add_to_top_thl: Optional[float] = field(
        default=None,
        repr=True,
        metadata={"serialize": True},
        init=True,
    )
    noise_std: Optional[float] = field(
        default=None,
        repr=True,
        metadata={"serialize": True},
        init=True,
    )
    noise_seed: Optional[int] = field(
        default=None,
        repr=True,
        metadata={"serialize": True},
        init=True,
    )
    noise_boundaries: Optional[List[str]] = field(
        default=None,
        repr=True,
        metadata={"serialize": True},
        init=True,
    )
    noise_variables: Optional[List[str]] = field(
        default=None,
        repr=True,
        metadata={"serialize": True},
        init=True,
    )
    noise_minzt: Optional[float] = field(
        default=None,
        repr=True,
        metadata={"serialize": True},
        init=True,
    )
    noise_maxzt: Optional[float] = field(
        default=None,
        repr=True,
        metadata={"serialize": True},
        init=True,
    )
    atmosphere_module_name: Optional[str] = field(
        default=None,
        repr=True,
        metadata={"serialize": True},
        init=True,
    )
    atmosphere_module: Optional[AtmosphereModule] = field(
        default=None,
        repr=False,
        metadata={"serialize": True},
        init=True,
    )


@register_module
@dataclass
class do_openboundary(simulation_module):
    sim: Optional["simulation_module"] = field(default=None, repr=False)
    openBCgrid: GridDalesOpenBC = field(
        default=None, repr=False, metadata={"serialize": False}, init=False
    )
    indices: Optional[dict] = field(
        default=None, repr=False, metadata={"serialize": False}, init=False
    )
    harmonieprepper: Optional[prep_harmonie.harmoniePrepper] = field(
        default=None, repr=False, metadata={"serialize": False}, init=False
    )
    boundaries: Optional[Any] = field(
        default=None, repr=False, metadata={"serialize": False}, init=False
    )
    initfields: Optional[Any] = field(
        default=None, repr=False, metadata={"serialize": False}, init=False
    )
    nest_in_harmonie: Optional[Nest_in_Harmonie] = field(
        default=None, repr=False, metadata={"serialize": True}, init=True
    )
    nest_in_dales: Optional[Nest_in_Dales] = field(
        default=None, repr=False, metadata={"serialize": True}, init=True
    )
    nest_in_atmosphere: Optional[Nest_in_AtmosphereProfiles] = field(
        default=None, repr=False, metadata={"serialize": True}, init=True
    )

    e12: Optional[float] = field(
        default=None, repr=True, metadata={"serialize": True}, init=True
    )
    tracernames: Optional[List[str]] = field(
        default_factory=list, repr=True, metadata={"serialize": True}, init=True
    )
    tchunk: Optional[int] = field(
        default=None, repr=True, metadata={"serialize": True}, init=True
    )
    start: Optional[str] = field(
        default=None, repr=True, metadata={"serialize": True}, init=True
    )
    time0: Optional[str] = field(
        default=None, repr=True, metadata={"serialize": True}, init=True
    )
    end: Optional[str] = field(
        default=None, repr=True, metadata={"serialize": True}, init=True
    )
    lopenbc: bool = field(
        default=True,
        repr=True,
        metadata={"nml": "OPENBC", "key": "lopenbc", "serialize": True},
        init=True,
    )
    linithetero: bool = field(
        default=True,
        repr=True,
        metadata={"nml": "OPENBC", "key": "linithetero", "serialize": True},
        init=True,
    )
    lper: bool = field(
        default=False,
        repr=True,
        metadata={"nml": "OPENBC", "key": "lper", "serialize": True},
        init=True,
    )
    lbuoytop: bool = field(
        default=True,
        repr=True,
        metadata={"nml": "OPENBC", "key": "lbuoytop", "serialize": True},
        init=True,
    )
    dxint: int = field(
        default=100,
        repr=True,
        metadata={"nml": "OPENBC", "key": "dxint", "serialize": True},
        init=True,
    )
    dyint: int = field(
        default=100,
        repr=True,
        metadata={"nml": "OPENBC", "key": "dyint", "serialize": True},
        init=True,
    )
    dzint: int = field(
        default=-1,
        repr=True,
        metadata={"nml": "OPENBC", "key": "dzint", "serialize": True},
        init=True,
    )
    dxturb: int = field(
        default=-1,
        repr=True,
        metadata={"nml": "OPENBC", "key": "dxturb", "serialize": True},
        init=True,
    )
    dyturb: int = field(
        default=-1,
        repr=True,
        metadata={"nml": "OPENBC", "key": "dyturb", "serialize": True},
        init=True,
    )
    taum: int = field(
        default=0,
        repr=True,
        metadata={"nml": "OPENBC", "key": "taum", "serialize": True},
        init=True,
    )
    tauh: int = field(
        default=20,
        repr=True,
        metadata={"nml": "OPENBC", "key": "tauh", "serialize": True},
        init=True,
    )
    pbc: int = field(
        default=3,
        repr=True,
        metadata={"nml": "OPENBC", "key": "pbc", "serialize": True},
        init=True,
    )
    lsynturb: bool = field(
        default=False,
        repr=True,
        metadata={"nml": "OPENBC", "key": "lsynturb", "serialize": True},
        init=True,
    )
    iturb: int = field(
        default=0,
        repr=True,
        metadata={"nml": "OPENBC", "key": "iturb", "serialize": True},
        init=True,
    )
    tau: int = field(
        default=180,
        repr=True,
        metadata={"nml": "OPENBC", "key": "tau", "serialize": True},
        init=True,
    )
    lambda_: int = field(
        default=1,
        repr=True,
        metadata={"nml": "OPENBC", "key": "lambda", "serialize": True},
        init=True,
    )
    nmodes: int = field(
        default=100,
        repr=True,
        metadata={"nml": "OPENBC", "key": "nmodes", "serialize": True},
        init=True,
    )
    lambdas: int = field(
        default=1,
        repr=True,
        metadata={"nml": "OPENBC", "key": "lambdas", "serialize": True},
        init=True,
    )
    lambdas_x: int = field(
        default=1,
        repr=True,
        metadata={"nml": "OPENBC", "key": "lambdas_x", "serialize": True},
        init=True,
    )
    lambdas_y: int = field(
        default=1,
        repr=True,
        metadata={"nml": "OPENBC", "key": "lambdas_y", "serialize": True},
        init=True,
    )
    lambdas_z: int = field(
        default=1,
        repr=True,
        metadata={"nml": "OPENBC", "key": "lambdas_z", "serialize": True},
        init=True,
    )
    igrw_damp: int = field(
        default=0,
        repr=True,
        metadata={"nml": "PHYSICS", "key": "igrw_damp", "serialize": True},
        init=True,
    )
    lconstexner: bool = field(
        default=True,
        repr=True,
        metadata={
            "nml": "thermodynamics",
            "key": "lconstexner",
            "serialize": True,
            "raise_conflict": True,
        },
        init=True,
    )
    lbaseexner: bool = field(
        default=True,
        repr=True,
        metadata={
            "nml": "thermodynamics",
            "key": "lbaseexner",
            "serialize": True,
            "raise_conflict": True,
        },
        init=True,
    )

    def __post_init__(self):
        super().__init__(self.sim)
        self.module_name = "do_openboundary"

    def __add__(self, obj) -> "do_openboundary":
        """Add configurations to open boundary module.

        Args:
            obj: Nest_in_Harmonie, Nest_in_Dales, Nest_in_AtmosphereProfiles

        Returns:
            self for chaining
        """
        if isinstance(obj, Nest_in_Harmonie):
            self.nest_in_harmonie = obj
        elif isinstance(obj, Nest_in_Dales):
            self.nest_in_dales = obj
        elif isinstance(obj, Nest_in_AtmosphereProfiles):
            self.nest_in_atmosphere = obj

        else:
            raise TypeError(
                f"Expected Nest_in_Harmonie/Nest_in_Dales/Nest_in_AtmosphereProfiles, got {type(obj)}"
            )
        return self

    def __iadd__(self, modification) -> "do_openboundary":
        """In-place addition of modification."""
        return self.__add__(modification)

    def prepare_calculation(self):
        if not type(self.grid) is GridDalesOpenBC:
            logger.warning(
                "Passing wrong grid type to openboundary, converting to GridDalesOpenBC with as_openBC."
            )
            self.openBCgrid = self.grid.as_openbc()
        else:
            self.openBCgrid = self.grid
        if self.nest_in_harmonie is not None:
            self._prepare_from_harmonie()
        elif self.nest_in_dales is not None:
            self._prepare_from_dales()
        elif self.nest_in_atmosphere is not None:
            self._prepare_from_atmosphere()
        else:
            raise ValueError("Unknown source for open boundary conditions")

    def _prepare_from_harmonie(self) -> None:
        config = {
            "openboundary": {
                "e12": self.e12,
                "tracernames": self.tracernames,
                "tchunk": self.tchunk,
                # "iexpnr": self.exp_id,
                "start": self.start,
                "author": "author",
                "time0": self.time0,
                "end": self.end,
                "HARMONIE_ml_glob": self.nest_in_harmonie.ml_glob,
                "HARMONIE_sfc_glob": self.nest_in_harmonie.sfc_glob,
            }
        }
        self.harmonieprepper = prep_harmonie.harmoniePrepper(
            config["openboundary"], self.openBCgrid
        )
        self.harmonieprepper.load_data()
        data, transform = self.harmonieprepper.prep_harmonie()
        # we need to use the right surface pressure as calculated from the input data
        logger.info("Setting namelist NAMSURFACE:ps to %f", self.harmonieprepper.ps)

        self.set_nml_section(
            "NAMSURFACE",
            "ps",
            self.harmonieprepper.ps,
        )
        # we need to use the right skin liq. pot. temperature from the input data
        logger.info("Setting namelist NAMSURFACE:thls to %f", self.harmonieprepper.thls)
        self.set_nml_section(
            "NAMSURFACE",
            "thls",
            self.harmonieprepper.thls,
        )
        (data,) = dask.optimize(data)

        logger.debug("Setting up boundary fields")
        self.boundaries = harmonie_boundary.boundary_fields(
            config["openboundary"],
            self.openBCgrid,
            data,
            output_path=self.output_path,
        )
        logger.debug("Setting up initial fields")
        self.initfields = initfields.initial_fields(
            config["openboundary"],
            self.openBCgrid,
            data,
            transform,
            output_path=self.output_path,
        )
        logger.debug("Setup all openBC fields, optimizing fields now..")
        self.boundaries, self.initfields = dask.optimize(
            self.boundaries, self.initfields
        )
        logger.debug("Optimized fields")

    def _prepare_from_dales(self) -> None:
        config = {
            "openboundary": {
                "e12": self.e12,
                "tracernames": self.tracernames,
                "tchunk": self.tchunk,
                "start": self.start,
                "time0": self.time0,
                "author": "author",
                # "iexpnr": self.exp_id,
                "end": self.end,
                # source of boundary fields for everything but t=0
                "outpath_coarse": self.nest_in_dales.outpath_coarse,
                # Source of initial boundary fields from a previous simulation, specifically,
                # the last time step in the output of the previous simulation
                "outpath_coarse_old": self.nest_in_dales.outpath_coarse_old,
                # source of initial boundary fields from a previous simulation,
                # specifically, t=0 of the previous simulation
                "inpath_coarse": self.nest_in_dales.inpath_coarse,
                # Source of initial fields for the current simulation, specifically
                # the initfields.nc input of the previous simulation
                "inpath": self.nest_in_dales.inpath,
            }
        }
        if self.nest_in_dales.inpath is not None:
            self.initfields = initial_fields_fine.initial_fields_fine(
                config["openboundary"],
                grid=self.openBCgrid,
                output_path=self.output_path,
            )
        else:

            self.initfields = xr.Dataset(coords={"time": [0]})

        self.boundaries = boundary_fields_fine.boundary_fields_fine(
            config["openboundary"],
            grid=self.openBCgrid,
            output_path=self.output_path,
            grid_indices=self.indices,
        )

    def _prepare_from_atmosphere(self) -> None:
        # Build open boundaries directly from vertical profiles configured
        # in the AtmosphereModule of this simulation. This provides
        # horizontally homogeneous (uniform) lateral boundary conditions
        # without requiring a separate coarse DALES run.

        atmo_module: Optional[AtmosphereModule] = (
            self.nest_in_atmosphere.atmosphere_module
        )
        if atmo_module is None:
            raise ValueError(
                "Nest_in_AtmosphereProfiles requires 'atmosphere_module' to be set; "
                "using AtmosphereModule instances from dales_simulation is not supported."
            )
        if atmo_module.sim is None or atmo_module.sim.grid is None:
            atmo_module._initialize_from_sim(self.sim)

        # Ensure profiles are evaluated; in the normal pipeline this
        # should already be true, but we guard against custom usage.
        if not atmo_module.prepare_calculation_done:
            atmo_module.prepare_calculation()
            atmo_module.prepare_calculation_done = True

        mapping = self.nest_in_atmosphere.variable_mapping or {}
        default_mapping = {
            "u": "ug",
            "v": "vg",
            "w": "w",
            "thl": "thetal",
            "qt": "qt",
            "e12": "tke",
        }
        # Fill in any missing keys from defaults without overwriting
        for key, val in default_mapping.items():
            mapping.setdefault(key, val)

        profiles_1d: Dict[str, np.ndarray] = {}
        for obc_var, atmo_name in mapping.items():
            atmo_definition = get_var_by_name()[atmo_name]
            if atmo_definition not in atmo_module.variables:
                raise ValueError(
                    f"AtmosphereModule is missing required profile '{atmo_name}' for open boundary variable '{obc_var}'"
                )
            var_container = atmo_module.variables[atmo_definition]
            if var_container.values is None:
                raise ValueError(
                    f"AtmosphereModule variable '{atmo_name}' has no evaluated values; "
                    "ensure AtmosphereModule.prepare_calculation() has run."
                )
            profiles_1d[obc_var] = np.asarray(var_container.values, dtype=float)

        # Sanity check vertical length
        nz = len(self.openBCgrid.zt)
        for obc_var, prof in profiles_1d.items():
            if prof.shape[0] != nz:
                raise ValueError(
                    f"Profile for '{obc_var}' has length {prof.shape[0]}, expected {nz} (len(grid.zt))"
                )

        # Helper to build a single boundary field for a given variable
        def _build_uniform_boundary(
            prof_z: np.ndarray, var: str, bnd: str
        ) -> xr.DataArray:
            # Reproduce the dimensional layout used by load_any_boundary_var
            var_dims_dic = {
                "u": {"x": "xm", "y": "yt", "z": "zt"},
                "v": {"x": "xt", "y": "ym", "z": "zt"},
                "w": {"x": "xt", "y": "yt", "z": "zm"},
                "default": {"x": "xt", "y": "yt", "z": "zt"},
            }
            boundary_var_assign_dic = {
                "west": {"default": ["z", "y"]},
                "east": {"default": ["z", "y"]},
                "south": {"default": ["z", "x"]},
                "north": {"default": ["z", "x"]},
                "top": {"default": ["y", "x"]},
            }
            assign_grid_dic = {
                "xt": self.openBCgrid.xt,
                "xm": self.openBCgrid.xm,
                "yt": self.openBCgrid.yt,
                "ym": self.openBCgrid.ym,
                "zt": self.openBCgrid.zt,
                "zm": self.openBCgrid.zm,
            }

            if var in var_dims_dic:
                var_dims = var_dims_dic[var]
            else:
                var_dims = var_dims_dic["default"]

            if var in boundary_var_assign_dic[bnd]:
                base_dims = boundary_var_assign_dic[bnd][var]
            else:
                base_dims = boundary_var_assign_dic[bnd]["default"]

            dims = [var_dims[d] for d in base_dims]

            # Determine vertical dimension, if any, and map profile to it
            vertical_dim = None
            for d in ("zt", "zm"):
                if d in dims:
                    vertical_dim = d
                    break

            if vertical_dim is None:
                # Top boundary: use value at model top
                fill_value = float(prof_z[-1])
                shape = [len(assign_grid_dic[d]) for d in dims]
                data = np.full(shape, fill_value, dtype=float)
            else:
                if vertical_dim == "zt":
                    prof_on_vert = prof_z
                    vert_coords = assign_grid_dic["zt"]
                else:  # zm
                    # Interpolate from zt -> zm
                    prof_on_vert = np.interp(
                        assign_grid_dic["zm"], self.openBCgrid.zt, prof_z
                    )
                    vert_coords = assign_grid_dic["zm"]

                shape = [len(assign_grid_dic[d]) for d in dims]
                data = np.zeros(shape, dtype=float)
                vert_index = dims.index(vertical_dim)
                for k in range(len(vert_coords)):
                    idx = [slice(None)] * len(shape)
                    idx[vert_index] = k
                    data[tuple(idx)] = prof_on_vert[k]

            coords = {d: assign_grid_dic[d] for d in dims}
            return xr.DataArray(data, coords=coords, dims=dims)

        # Determine time sampling from AtmosphereModule timedep forcings (if any).
        # AtmosphereModule currently returns a dict keyed by VariableDefinition,
        # so we normalise keys to plain variable names here.
        timed_forcings = atmo_module.get_timedep_atmosphere_forcings() or {}
        timed_forcings_by_name: Dict[str, Dict[float, np.ndarray]] = {}
        for key, series in timed_forcings.items():
            if isinstance(key, str):
                name = key
            else:
                name = getattr(key, "name", None)
            if name is None:
                continue
            # Merge series for the same variable name, later entries override.
            target = timed_forcings_by_name.setdefault(name, {})
            for t, values in series.items():
                target[float(t)] = np.asarray(values, dtype=float)

        times_set = {0.0}
        for obc_var, atmo_name in mapping.items():
            series = timed_forcings_by_name.get(atmo_name)
            if series:
                times_set.update(float(t) for t in series.keys())

        all_times = sorted(times_set)

        base_time_str = self.time0 or self.start
        if base_time_str is None:
            raise ValueError(
                "Nest_in_AtmosphereProfiles requires 'time0' or 'start' to be set on do_openboundary"
            )

        time_points = [
            np.datetime64(base_time_str) + np.timedelta64(int(round(t)), "s")
            for t in all_times
        ]

        ds = xr.Dataset(coords={"time": ("time", time_points)})
        ds = ds.assign_coords(
            {
                "xt": ("xt", self.openBCgrid.xt),
                "xm": ("xm", self.openBCgrid.xm),
                "yt": ("yt", self.openBCgrid.yt),
                "ym": ("ym", self.openBCgrid.ym),
                "zt": ("zt", self.openBCgrid.zt),
                "zm": ("zm", self.openBCgrid.zm),
            }
        )

        boundaries = ["west", "east", "south", "north", "top"]
        base_vars = ["u", "v", "w", "thl", "qt", "e12"]
        for var in mapping.keys():
            if var not in base_vars:
                logger.info(f"Adding variable '{var}' to base_vars")
                base_vars.append(var)

        add_to_top_thl = getattr(self.nest_in_atmosphere, "add_to_top_thl", None)
        for var in base_vars:
            atmo_name = mapping[var]
            series = timed_forcings_by_name.get(atmo_name, {})
            var_slices = []
            for t in all_times:
                if t in series:
                    prof_z = np.asarray(series[t], dtype=float)
                else:
                    prof_z = profiles_1d[var]
                for bnd in boundaries:
                    if var == "thl" and bnd == "top" and add_to_top_thl is not None:
                        prof_z = prof_z.copy()
                        prof_z[-1] += add_to_top_thl
                    da2d = _build_uniform_boundary(prof_z, var, bnd)
                    da3d = da2d.expand_dims(
                        {
                            "time": [
                                np.datetime64(base_time_str)
                                + np.timedelta64(int(round(t)), "s")
                            ]
                        }
                    )
                    var_slices.append((bnd, da3d))

            # Combine along time for each boundary
            by_boundary: Dict[str, List[xr.DataArray]] = {}
            for bnd, da in var_slices:
                by_boundary.setdefault(bnd, []).append(da)
            for bnd, da_list in by_boundary.items():
                ds[f"{var}{bnd}"] = xr.concat(da_list, dim="time")

        # Optionally add random noise to the Atmosphere-based open boundaries.
        # This is useful to break perfect homogeneity and seed variability
        # in the inflow. Noise is applied independently at each (time, x, y, z)
        # point, with a Gaussian distribution.
        noise_std = getattr(self.nest_in_atmosphere, "noise_std", None)
        noise_minzt = getattr(self.nest_in_atmosphere, "noise_minzt", None)
        noise_maxzt = getattr(self.nest_in_atmosphere, "noise_maxzt", None)
        if noise_std is not None and (
            noise_minzt is not None or noise_maxzt is not None
        ):
            zt = self.openBCgrid.zt
            zm = self.openBCgrid.zm
            if noise_minzt is not None:
                mask = zt >= noise_minzt
                mask_zm = zm >= noise_minzt
            else:
                mask = np.ones_like(zt, dtype=bool)
                mask_zm = np.ones_like(zm, dtype=bool)
            if noise_maxzt is not None:
                mask &= zt <= noise_maxzt
                mask_zm &= zm <= noise_maxzt
            if (not np.any(mask)) or (not np.any(mask_zm)):
                logger.warning(
                    "Noise std is set but no vertical levels are within the specified min/max zt bounds; skipping noise addition."
                )
                noise_std = None
            else:
                logger.info(
                    "Applying noise with std=%.3f to levels where zt is between %.2f and %.2f",
                    noise_std,
                    zt[mask][0],
                    zt[mask][-1],
                )
        if noise_std is not None and noise_std > 0.0:
            rng = np.random.default_rng(
                getattr(self.nest_in_atmosphere, "noise_seed", None)
            )
            requested_bounds = getattr(
                self.nest_in_atmosphere, "noise_boundaries", None
            )
            requested_vars = getattr(self.nest_in_atmosphere, "noise_variables", None)
            if requested_bounds is None:
                active_bounds = set(boundaries)
            else:
                active_bounds = {b for b in requested_bounds if b in boundaries}
            if requested_vars is None:
                active_vars = set(base_vars)
            else:
                active_vars = {v for v in requested_vars if v in base_vars}
            for var in active_vars:
                for bnd in active_bounds:
                    name = f"{var}{bnd}"
                    if name not in ds:
                        continue
                    arr = ds[name]
                    mask_3d = xr.ones_like(arr, dtype=bool)
                    # Apply vertical mask to avoid adding noise outside the specified zt bounds
                    if "zt" in arr.dims:
                        mask_3d = mask_3d.where(ds.zt >= noise_minzt, other=False)
                        mask_3d = mask_3d.where(ds.zt <= noise_maxzt, other=False)
                    if "zm" in arr.dims:
                        mask_3d = mask_3d.where(ds.zm >= noise_minzt, other=False)
                        mask_3d = mask_3d.where(ds.zm <= noise_maxzt, other=False)
                    noise = rng.normal(loc=0.0, scale=noise_std, size=arr.shape)
                    # Only add noise to the active region defined by the vertical mask
                    noise *= mask_3d
                    ds[name] = arr + noise
        # Tracers: currently default to zero fields with scalar-like layout
        # if self.tracernames:
        #     zero_prof = np.zeros_like(self.openBCgrid.zt, dtype=float)
        #     for tracer in self.tracernames:
        #         tracer_slices = []
        #         for t in all_times:
        #             for bnd in boundaries:
        #                 da2d = _build_uniform_boundary(zero_prof, "qt", bnd)
        #                 da3d = da2d.expand_dims(
        #                     {
        #                         "time": [
        #                             np.datetime64(base_time_str)
        #                             + np.timedelta64(int(round(t)), "s")
        #                         ]
        #                     }
        #                 )
        #                 tracer_slices.append((bnd, da3d))
        #         by_boundary: Dict[str, List[xr.DataArray]] = {}
        #         for bnd, da in tracer_slices:
        #             by_boundary.setdefault(bnd, []).append(da)
        #         for bnd, da_list in by_boundary.items():
        #             ds[f"{tracer}{bnd}"] = xr.concat(da_list, dim="time")

        config = {
            "openboundary": {
                "e12": self.e12,
                "tracernames": self.tracernames,
                "tchunk": self.tchunk,
                "start": self.start,
                "time0": self.time0,
                "author": "author",
                "end": self.end,
            }
        }

        # Reuse attribute/encoding setup from existing helper
        self.initfields = xr.Dataset(coords={"time": [0]})
        self.boundaries = boundary_fields_fine.set_openboundary_attrs(
            config["openboundary"], ds
        )

    def write_files(self):
        # Save data
        logger.debug("Optimizing initial fields")
        self.initfields = self.initfields.assign_attrs(
            {
                "title": f"initfields.inp.{self.exp_id:03d}.nc",
            }
        )
        initfields_writer = self.initfields.to_netcdf(
            path=self.output_path / "input" / self.initfields.attrs["title"],
            mode="w",
            format="netcdf4",
            compute=False,
        )
        (initfields_writer,) = dask.optimize(initfields_writer)
        logger.debug("Writing initial fields")

        initfields_writer.compute()

        self.boundaries = self.boundaries.assign_attrs(
            {
                "title": f"openboundaries.inp.{self.exp_id:03d}.nc",
            }
        )
        logger.debug("Optimizing boundaries")
        openboundaries_writer = self.boundaries.to_netcdf(
            path=self.output_path / "input" / self.boundaries.attrs["title"],
            mode="w",
            format="netcdf4",
            compute=False,
        )
        (openboundaries_writer,) = dask.optimize(openboundaries_writer)
        logger.debug("Writing boundaries")

        openboundaries_writer.compute()
