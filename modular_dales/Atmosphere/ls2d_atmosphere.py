import logging
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, List, Optional

import netCDF4 as nc4
import numpy as np
import time
import logging
import random
from typing import Callable, Type, Tuple
from modular_dales.MODULE_REGISTRY import register_module, register_singleton
from modular_dales.modular.simulation_module import simulation_module
from modular_dales.modular.time_dependent_scalars import TimeDependentScalar
from modular_dales.vars import (
    VariableDefinition,
    thls,
    ua,
    va,
    w,
    thetal,
    qt,
    tke,
    ug,
    vg,
    wa,
    dqtdxls,
    dqtdyls,
    tnqt_adv,
    tnthetal_rad,
    dudt_ls,
    dvdt_ls,
    wtsurf,
    wqsurf,
    qtsurf,
    psurf,
    qnetav,
)


logger = logging.getLogger(__name__)


def retry_exponential(
    func: Callable,
    *,
    max_attempts: int = 5,
    max_total_time: float = 60.0,
    base_delay: float = 1.0,
    max_delay: float = 30.0,
):
    """
    Retry a function with exponential backoff.

    Args:
        func: Function to execute
        max_attempts: Maximum number of attempts
        max_total_time: Maximum total retry time (seconds)
        base_delay: Initial delay in seconds
        max_delay: Maximum delay between retries
    """
    start_time = time.monotonic()
    attempt = 0

    while True:
        attempt += 1
        if func():
            return
        else:
            elapsed = time.monotonic() - start_time

            if attempt >= max_attempts:
                logger.error(
                    "Max attempts reached (%s). Giving up download.",
                    max_attempts,
                )
                raise RuntimeError("Max download attempts reached. Giving up download.")

            if elapsed >= max_total_time:
                logger.error(
                    "Max total retry time exceeded (%.2fs). Giving up download.",
                    max_total_time,
                )
                raise RuntimeError("Max total retry time exceeded. Giving up download.")

            delay = min(base_delay * (2 ** (attempt - 1)), max_delay)

            logger.warning(
                "Download attempt %s failed. (data not yet ready) Retrying in %.2fs...",
                attempt,
                delay,
            )

            time.sleep(delay)
@register_singleton
@register_module
@dataclass
class FromLS2D:
    """Marker class to enable LS2D-driven soil/roughness in LSM.
    Also enables injection of LS2D time series into atmosphere.

    When added to :class:`LSMModule` or :class:`TimeDependentModule` (``lsm += FromLS2D()``), the
    module will, if an :class:`LS2DAtmosphereModule` is present,
    override soil temperature/moisture and soil type index from LS2D
    and also set the bulk roughness lengths ``z0mav`` and ``z0hav``
    from LS2D time-mean values.
    """


@register_module
@dataclass
class LS2DAtmosphereModule(simulation_module):
    """Atmosphere forcing module using LS2D ERA5 processing.

      * Uses LS2D to download and process ERA5 data.
      * Calls ``era.calculate_forcings`` in LS2D.
      * Interpolates the resulting forcings onto the DALES
        vertical grid given by :class:`GridDales` (``grid.zt``).
            * Provides large-scale forcings and base profiles that can be
                injected into :class:`AtmosphereModule` and surface modules.
    """

    sim: Optional["simulation_module"] = field(default=None, repr=False)

    # LS2D configuration (similar to jupyter_tests/jupyter_test.py)
    central_lat: Optional[float] = None
    central_lon: Optional[float] = None
    area_size: float = 1.0
    case_name: str = "ls2d_case"
    era5_path: Optional[str] = None
    start_date: Optional[datetime] = None
    end_date: Optional[datetime] = None
    write_log: bool = True
    data_source: str = "CDS"

    # LS2D forcing calculation options
    n_av: int = 0
    method: str = "2nd"

    # Initial turbulence level (cf. jupyter_tests/jupyter_test.py)
    init_tke: float = 0.1

    # Runtime-only helpers (not serialized)
    les_input: Any = field(
        default=None,
        init=False,
        repr=False,
        metadata={"serialize": False},
    )
    _times_with_zero: List[float] = field(
        default_factory=list,
        init=False,
        repr=False,
        metadata={"serialize": False},
    )
    _timedep_var_dic: Dict[VariableDefinition, np.ndarray] = field(
        default_factory=dict,
        init=False,
        repr=False,
        metadata={"serialize": False},
    )

    def __post_init__(self):
        super().__init__(self.sim)
        self.module_name = "LS2DAtmosphereModule"

    def check_settings(self):
        """Basic validation of configuration before running LS2D."""

        if self.grid is None:
            raise ValueError(
                "LS2DAtmosphereModule requires a GridDales grid (sim.grid) to be set"
            )

        # Derive central_lat/central_lon from the grid if not set
        if self.central_lat is None and getattr(self.grid, "xlat", None) is not None:
            self.central_lat = float(self.grid.xlat)
        if self.central_lon is None and getattr(self.grid, "xlon", None) is not None:
            self.central_lon = float(self.grid.xlon)

        # Derive case_name from the parent simulation if not set
        if (not self.case_name) and getattr(self.sim, "case_name", None):
            self.case_name = str(self.sim.case_name)

        missing: List[str] = []
        for name in ("central_lat", "central_lon", "era5_path", "start_date", "end_date"):
            if getattr(self, name) is None:
                missing.append(name)
        if missing:
            raise ValueError(
                "LS2DAtmosphereModule missing required settings: " + ", ".join(missing)
            )
        return None

    def _import_ls2d(self):
        """Import ls2d lazily so tests that don't use it still run.

        Raises a clear error when ls2d is not available.
        """

        try:
            import ls2d  # type: ignore[import]
        except Exception as exc:  # pragma: no cover - environment dependent
            raise RuntimeError(
                "LS2DAtmosphereModule requires the 'ls2d' package to be installed and importable"
            ) from exc
        return ls2d

    def prepare_calculation(self):
        """Run LS2D, build ``les_input`` on ``GridDales.zt`` and cache forcings.

        The resulting time-dependent fields are stored internally in a
        form that can be used both for writing ``forcings.<exp_id>.nc``
        and for injecting time series into other modules.

        * Surface scalars are stored as 1D arrays with length ``nt``.
        * Profile variables are stored as 2D arrays of shape
            ``(nz, nt)`` where first axis is height and second axis is time.

        The time axis is taken from ``les_input.time_sec`` and cached in
        :attr:`_times_with_zero`.
        """

        if self.grid is None:
            raise ValueError(
                "LS2DAtmosphereModule requires grid (GridDales) to be initialized before prepare_calculation"
            )

        ls2d = self._import_ls2d()

        settings = {
            "central_lat": self.central_lat,
            "central_lon": self.central_lon,
            "area_size": self.area_size,
            "case_name": self.case_name,
            "era5_path": self.era5_path,
            "start_date": self.start_date,
            "end_date": self.end_date,
            "write_log": self.write_log,
            "data_source": self.data_source,
        }

        logger.info("LS2DAtmosphereModule: downloading ERA5 via LS2D")
        
        retry_exponential(lambda: ls2d.download_era5(settings, exit_when_waiting=False), max_attempts=30, max_total_time=3600, base_delay=60, max_delay=180)

        logger.info("LS2DAtmosphereModule: reading ERA5 via LS2D.Read_era5")
        era = ls2d.Read_era5(settings)
        era.calculate_forcings(n_av=self.n_av, method=self.method)

        # Interpolate ERA5 variables and forcings onto the DALES grid
        z_dales = np.asarray(self.grid.zt, dtype=float)
        logger.info(
            "LS2DAtmosphereModule: building les_input on DALES grid (k=%d)",
            z_dales.size,
        )
        les_input = era.get_les_input(z_dales)
        self.les_input = les_input

        times_sec = np.asarray(les_input.time_sec.values, dtype=float)
        if times_sec.ndim != 1 or times_sec.size < 2:
            raise ValueError(
                "LS2DAtmosphereModule expects les_input.time_sec to be 1D with at least 2 entries"
            )
        self._times_with_zero = times_sec.tolist()

        nt = times_sec.size
        nz = z_dales.size

        # ------------------------------------------------------------------
        # Base (time-independent) profiles from LS2D at time index 0
        # ------------------------------------------------------------------
        base_profiles: Dict[VariableDefinition, np.ndarray] = {}

        # Helper to extract the first time slice of a (time, z) field
        def _first_time_profile(name: str) -> Optional[np.ndarray]:
            if not hasattr(les_input, name):
                return None
            values_ = np.asarray(getattr(les_input, name).values, dtype=float)
            if values_.ndim == 2:
                # Expecting (time, z)
                if values_.shape[0] == nt and values_.shape[1] == nz:
                    return values_[0, :]
                if values_.shape[1] == nt and values_.shape[0] == nz:
                    # Transposed convention (z, time)
                    return values_[:, 0]
                raise ValueError(
                    f"LS2DAtmosphereModule: unexpected shape for initial les_input.{name}: "
                    f"{values_.shape}, expected ({nt}, {nz}) or ({nz}, {nt})"
                )
            raise ValueError(
                f"LS2DAtmosphereModule: initial profile les_input.{name} must be 2D, got {values_.shape}"
            )

        prof_thetal = _first_time_profile("thl")
        prof_qt = _first_time_profile("qt")
        prof_u = _first_time_profile("u")
        prof_v = _first_time_profile("v")
        prof_w = _first_time_profile("wls")

        if prof_thetal is not None:
            base_profiles[thetal] = prof_thetal
        if prof_qt is not None:
            base_profiles[qt] = prof_qt
        if prof_u is not None:
            base_profiles[ua] = prof_u
        if prof_v is not None:
            base_profiles[va] = prof_v
        if prof_w is not None:
            base_profiles[w] = prof_w

        # Vertical velocity is not provided by LS2D as an initial
        # profile; follow jupyter_test convention and set it to zero.
        base_profiles[w] = np.zeros(nz, dtype=float)

        # Turbulent kinetic energy: jupyter_test uses a constant
        # value over height.
        base_profiles[tke] = np.ones(nz, dtype=float) * float(self.init_tke)

        def _reshape_profile(name: str) -> Optional[np.ndarray]:
            """Return profile as (nz, nt) or scalar series as (nt,)."""

            if not hasattr(les_input, name):
                return None

            values = np.asarray(getattr(les_input, name).values, dtype=float)

            if values.ndim == 2:
                # LS2D convention: (time, z) -> transpose to (z, time)
                if values.shape == (nt, nz):
                    return values.T
                if values.shape == (nz, nt):
                    return values
                raise ValueError(
                    f"LS2DAtmosphereModule: unexpected shape for les_input.{name}: {values.shape}, "
                    f"expected ({nt}, {nz}) or ({nz}, {nt})"
                )

            if values.ndim == 1:
                if values.size != nt:
                    raise ValueError(
                        f"LS2DAtmosphereModule: unexpected length for les_input.{name}: "
                        f"{values.size}, expected {nt}"
                    )
                return values

            raise ValueError(
                f"LS2DAtmosphereModule: unsupported dimensions for les_input.{name}: {values.shape}"
            )

        timedep: Dict[VariableDefinition, np.ndarray] = {}

        # Large-scale wind and tendencies (profile variables)
        arr_ug = _reshape_profile("ug")
        if isinstance(arr_ug, np.ndarray) and arr_ug.ndim == 2:
            timedep[ug] = arr_ug

        arr_vg = _reshape_profile("vg")
        if isinstance(arr_vg, np.ndarray) and arr_vg.ndim == 2:
            timedep[vg] = arr_vg

        arr_wls = _reshape_profile("wls")
        if isinstance(arr_wls, np.ndarray) and arr_wls.ndim == 2:
            timedep[wa] = arr_wls

        # Advection tendencies for moisture and temperature
        arr_dtqt = _reshape_profile("dtqt_advec")
        if isinstance(arr_dtqt, np.ndarray) and arr_dtqt.ndim == 2:
            timedep[tnqt_adv] = arr_dtqt

        arr_dtthl = _reshape_profile("dtthl_advec")
        if isinstance(arr_dtthl, np.ndarray) and arr_dtthl.ndim == 2:
            timedep[tnthetal_rad] = arr_dtthl

        # Momentum tendencies
        arr_dtu = _reshape_profile("dtu_advec")
        if isinstance(arr_dtu, np.ndarray) and arr_dtu.ndim == 2:
            timedep[dudt_ls] = arr_dtu

        arr_dtv = _reshape_profile("dtv_advec")
        if isinstance(arr_dtv, np.ndarray) and arr_dtv.ndim == 2:
            timedep[dvdt_ls] = arr_dtv

        # Surface pressure and fluxes/temperatures from LS2D where available.
        # These correspond conceptually to ps, Ts, wth and qt time series.
        arr_ps = _reshape_profile("ps")
        if isinstance(arr_ps, np.ndarray) and arr_ps.ndim == 1:
            timedep[psurf] = arr_ps

        # Surface (skin) temperature: map to ``thls`` timedep variable and,
        # when possible, into a SurfaceModule via TimeDependentScalar.
        arr_ts = _reshape_profile("ts")
        if isinstance(arr_ts, np.ndarray) and arr_ts.ndim == 1:
            timedep[thls] = arr_ts
            self._inject_surface_series("thls", self._times_with_zero, arr_ts)

        # Sensible and latent heat fluxes: map to kinematic flux variables
        # and inject into ConstantFlux-style surface modules when present.
        arr_wth = _reshape_profile("wth")
        if isinstance(arr_wth, np.ndarray) and arr_wth.ndim == 1:
            timedep[wtsurf] = arr_wth
            self._inject_surface_series("wtsurf", self._times_with_zero, arr_wth)

        arr_wq = _reshape_profile("wq")
        if isinstance(arr_wq, np.ndarray) and arr_wq.ndim == 1:
            timedep[wqsurf] = arr_wq
            self._inject_surface_series("wqsurf", self._times_with_zero, arr_wq)

        # Surface total water mixing ratio: take lowest model level as proxy
        # when a dedicated surface field is not available.
        arr_qt_prof = _reshape_profile("qt")
        if isinstance(arr_qt_prof, np.ndarray) and arr_qt_prof.ndim == 2:
            timedep[qtsurf] = arr_qt_prof[0, :]

        # Net radiative flux at the surface is not provided directly by LS2D;
        # keep it as zeros for now so the variable exists.
        timedep.setdefault(qnetav, np.zeros(nt, dtype=float))

        # Placeholders for horizontal moisture gradients if they ever
        # become available in les_input (currently unused in jupyter_test).
        _ = dqtdxls, dqtdyls  # referenced to avoid "unused" warnings

        # Cache for later use (writing forcings and/or exposing to
        # TimedependentModule via hooks).
        self._timedep_var_dic = timedep

        # Inject LS2D-derived base profiles into an AtmosphereModule
        # if one is present in the simulation. That module will then
        # write ``init.<exp_id>.nc`` using its normal machinery.
        try:
            self._inject_base_profiles_into_atmosphere(base_profiles, z_dales)
        except Exception as exc:  # pragma: no cover - robustness
            logger.warning(
                "LS2DAtmosphereModule: failed to inject base profiles into AtmosphereModule: %s",
                exc,
            )
        try:
            self._inject_timed_profiles(timedep, z_dales)
        except Exception as exc:  # pragma: no cover - robustness
            logger.warning(
                "LS2DAtmosphereModule: failed to inject timed profiles into AtmosphereModule: %s",
                exc,
            )
        logger.info(
            "LS2DAtmosphereModule: prepared LS2D base profiles and forcings for %d times and %d levels",
            nt,
            nz,
        )
        return None
    def _inject_timed_profiles(
        self,
        timedep: Dict[VariableDefinition, np.ndarray],
        z_dales: np.ndarray,
    ) -> None:
        """Inject LS2D time-dependent profiles into an AtmosphereModule if present.

        For each variable in ``timedep`` that does not already have a
        configured profile in the AtmosphereModule, this helper appends a
        :class:`InterpolatedProfile` using the LS2D values at each time.
        """

        if not timedep or self.sim is None:
            return

        # Local import to avoid circular imports at module load time
        from modular_dales.Atmosphere.atmosphere import (
            AtmosphereModule,
            InterpolatedProfile,
        )

        if not self.module_exists(AtmosphereModule):
            return

        atmo: AtmosphereModule = self.retrieve_module(AtmosphereModule)

        z_list = [float(z) for z in z_dales]
        for time_idx, time in enumerate(self._times_with_zero):
            if not float(time) in atmo.collected_timed_profiles_by_time:
                atmo.collected_timed_profiles_by_time[float(time)] = {}
            for var_def, values in timedep.items():
                if var_def in atmo.collected_timed_profiles_by_time[float(time)]:
                    # Respect user-configured profiles; they can override LS2D.
                    continue
                if values.ndim != 2:
                    continue
                if values.ndim != 2 or values.shape[0] != len(z_list):
                    logger.warning(
                        "LS2DAtmosphereModule: skipping injection of timed profile for %s due to unexpected shape %s",
                        var_def.name,
                        values.shape,
                    )
                    continue

                profile_at_time = values[:,time_idx]
                atmo.collected_timed_profiles_by_time[float(time)][var_def] = InterpolatedProfile(
                    variable=var_def,
                    z=z_list,
                    points=[float(v) for v in profile_at_time],  # transpose to (time, z)
                )
    def write_files(self):
        """Write LS2D-derived files.

        This module no longer writes the ``init.<exp_id>.nc`` file
        itself; base profiles are injected into an
        :class:`AtmosphereModule` (when present), which then writes the
        init file. The radiation background is
        written to ``backrad.inp.<exp_id>.nc``.
        """

        if self.grid is None:
            return None

        output_input_path = self.output_path / "input"
        # Optional LS2D radiation background for RRTMG-style schemes
        # (backrad.inp.<exp_id>.nc). Only written when LS2D has been
        # run and provided the required variables.
        try:
            self._write_backrad_file(output_input_path)
        except Exception as exc:  # pragma: no cover - robustness
            logger.warning("LS2DAtmosphereModule: skipping backrad output: %s", exc)

        return None

    # ------------------------------------------------------------------
    # Extra helpers: backrad
    # ------------------------------------------------------------------

    def _write_backrad_file(self, output_input_path):
        """Write LS2D-based radiation background ``backrad.inp.<id>.nc``.

        This mirrors :func:`dales_ls2d_tools.create_backrad` used in the
        LS2D examples: time-mean profiles of ``p_lay``, ``t_lay`` and
        ``h2o_lay`` are written on a single vertical coordinate ``lev``.
        """

        if self.les_input is None:
            return None

        les_input = self.les_input
        if not all(hasattr(les_input, name) for name in ("p_lay", "t_lay", "h2o_lay")):
            return None

        p = np.asarray(les_input["p_lay"].values, dtype=float)
        T = np.asarray(les_input["t_lay"].values, dtype=float)
        q = np.asarray(les_input["h2o_lay"].values, dtype=float)

        # Expect (time, lev); fall back to averaging over the first axis
        # that matches the time dimension length.
        nt = len(self._times_with_zero) if self._times_with_zero else p.shape[0]

        def _time_mean(arr: np.ndarray) -> np.ndarray:
            if arr.ndim == 2:
                if arr.shape[0] == nt:
                    return arr.mean(axis=0)
                if arr.shape[1] == nt:
                    return arr.mean(axis=1)
            # Fallback: mean over first axis
            return arr.mean(axis=0)

        p_prof = _time_mean(p)
        T_prof = _time_mean(T)
        q_prof = _time_mean(q)

        filename = output_input_path / f"backrad.inp.{self.exp_id:03d}.nc"
        with nc4.Dataset(filename.as_posix(), "w") as nc_file:
            nc_file.createDimension("lev", p_prof.size)

            p_var = nc_file.createVariable("lev", "f4", ("lev",))
            T_var = nc_file.createVariable("T", "f4", ("lev",))
            q_var = nc_file.createVariable("q", "f4", ("lev",))

            p_var[:] = p_prof
            T_var[:] = T_prof
            q_var[:] = q_prof

        logger.info(
            "LS2DAtmosphereModule: wrote backrad.inp.%03d.nc to %s",
            self.exp_id,
            output_input_path,
        )
        return None

    # ------------------------------------------------------------------
    # Internal helpers for injection into other modules
    # ------------------------------------------------------------------

    def _inject_base_profiles_into_atmosphere(
        self,
        base_profiles: Dict[VariableDefinition, np.ndarray],
        z_dales: np.ndarray,
    ) -> None:
        """Inject LS2D base profiles into an AtmosphereModule if present.

        For each variable in ``base_profiles`` that does not already
        have a configured profile in the AtmosphereModule, this helper
        appends an :class:`InterpolatedProfile` using the LS2D values.
        """

        if not base_profiles or self.sim is None:
            return

        # Local import to avoid circular imports at module load time
        from modular_dales.Atmosphere.atmosphere import (
            AtmosphereModule,
            InterpolatedProfile,
        )

        if not self.module_exists(AtmosphereModule):
            return

        atmo = self.retrieve_module(AtmosphereModule)

        z_list = [float(z) for z in z_dales]

        for var_def, values in base_profiles.items():
            if var_def in atmo.collected_base_profiles:
                # Respect user-configured profiles; they can override LS2D.
                continue
            atmo.collected_base_profiles[var_def] = InterpolatedProfile(
                variable=var_def,
                z=z_list,
                points=[float(v) for v in values],
            )

    def _inject_surface_series(
        self,
        field_name: str,
        times: List[float],
        values: np.ndarray,
    ) -> None:
        """Inject a time series into a surface module, if compatible.

        When a surface module with a matching dataclass field name is
        present (e.g. ``thls`` for ConstantSurfaceTemperatureModule or
        ``wtsurf``/``wqsurf`` for flux modules), this helper assigns a
        :class:`TimeDependentScalar` constructed from ``times`` and
        ``values``. Existing non-None values are left untouched to
        allow manual overrides.
        """

        if self.sim is None:
            return

        # Construct TimeDependentScalar once
        try:
            series = TimeDependentScalar(times=list(times), values=list(values))
        except Exception:
            return

        for module in self.sim.modules:
            if module is self:
                continue
            if not hasattr(module, field_name):
                continue
            current = getattr(module, field_name)
            if current is not None:
                # Respect explicit user configuration
                continue
            setattr(module, field_name, series)
            # Only set on the first suitable module
            break
