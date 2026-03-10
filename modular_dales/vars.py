from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple, Dict


@dataclass(frozen=True)
class VariableDefinition:
    """Definition of an atmospheric variable.

    This is the user-facing object you should import and pass to
    AtmosphericProfile / InterpolatedProfile instead of raw strings.
    """

    name: str
    long_name: str
    unit: str
    can_be_time_dependent: bool = False
    must_only_be_time_dependent: bool = False
    is_profile: bool = True
    time_dependent_name: str | None = None


thls = VariableDefinition(
    "thls",
    "Surface liquid water potential temperature",
    "K",
    is_profile=False,
    can_be_time_dependent=True,
    time_dependent_name="thlsurf_timedep",
)
wtsurf = VariableDefinition(
    "wtsurf",
    "Surface kinematic heat flux",
    "W/m^2",
    is_profile=False,
    can_be_time_dependent=True,
    time_dependent_name="wtsurf_timedep",
)
wqsurf = VariableDefinition(
    "wqsurf",
    "Surface kinematic moisture flux",
    "kg/m^2/s",
    is_profile=False,
    can_be_time_dependent=True,
    time_dependent_name="wqsurf_timedep",
)
qtsurf = VariableDefinition(
    "qtsurf",
    "Surface total water mixing ratio",
    "kg/kg",
    is_profile=False,
    can_be_time_dependent=True,
    time_dependent_name="qtsurf_timedep",
)
psurf = VariableDefinition(
    "psurf",
    "Surface pressure",
    "Pa",
    is_profile=False,
    can_be_time_dependent=True,
    time_dependent_name="psurf_timedep",
)
qnetav = VariableDefinition(
    "qnetav",
    "Net radiative flux at the surface",
    "W/m^2",
    is_profile=False,
    can_be_time_dependent=True,
    time_dependent_name="qnetavsurf_timedep",
)
ua = VariableDefinition("ua", "Initial eastward velocity profile", "m/s")
va = VariableDefinition("va", "Initial northward velocity profile", "m/s")
w = VariableDefinition("w", "Initial vertical velocity profile", "m/s")
thetal = VariableDefinition(
    "thetal",
    "Initial liquid water potential temperature profile",
    "K",
)
qt = VariableDefinition(
    "qt",
    "Initial total water mixing ratio profile.",
    "kg/kg",
)
tke = VariableDefinition(
    "tke",
    "Initial profile of the square root of the turbulence kinetic energy (TKE)",
    "m/s",
)
ug = VariableDefinition(
    "ug",
    "Geostrophic eastward wind",
    "m/s",
    can_be_time_dependent=True,
    time_dependent_name="ug_timedep",
)
vg = VariableDefinition(
    "vg",
    "Geostrophic northward wind",
    "m/s",
    can_be_time_dependent=True,
    time_dependent_name="vg_timedep",
)
dpdx = VariableDefinition(
    "dpdx",
    "Eastward pressure gradient",
    "Pa/m",
    can_be_time_dependent=True,
    time_dependent_name="dpdx_ls_timedep",
)
dpdy = VariableDefinition(
    "dpdy",
    "Northward pressure gradient",
    "Pa/m",
    can_be_time_dependent=True,
    time_dependent_name="dpdy_ls_timedep",
)
wa = VariableDefinition(
    "wa",
    "Large-scale subsidence.",
    "m/s",
    can_be_time_dependent=True,
    time_dependent_name="wf_ls_timedep",
)
dqtdxls = VariableDefinition(
    "dqtdxls",
    "Eastward gradient of the total water mixing ratio due to advection",
    "kg/kg/m",
    can_be_time_dependent=True,
    time_dependent_name="dqtdx_ls_timedep",
)
dqtdyls = VariableDefinition(
    "dqtdyls",
    "Northward gradient of the total water mixing ratio due to advection",
    "kg/kg/m",
    can_be_time_dependent=True,
    time_dependent_name="dqtdy_ls_timedep",
)
tnqt_adv = VariableDefinition(
    "tnqt_adv",
    "Tendency of the total water mixing ratio",
    "kg/kg/s",
    can_be_time_dependent=True,
    time_dependent_name="dqtdt_ls_timedep",
)
tnthetal_rad = VariableDefinition(
    "tnthetal_rad",
    "Tendency of the liquid water potential temperature due to radiative heating",
    "K/s",
    can_be_time_dependent=True,
    time_dependent_name="dthl_rad_timedep",
)
dudt_ls = VariableDefinition(
    "dudt_ls",
    "Tendency of the eastward velocity due to large-scale forcing",
    "m/s^2",
    can_be_time_dependent=True,
    must_only_be_time_dependent=True,
    time_dependent_name="dudt_ls_timedep",
)
dvdt_ls = VariableDefinition(
    "dvdt_ls",
    "Tendency of the northward velocity due to large-scale forcing",
    "m/s^2",
    can_be_time_dependent=True,
    must_only_be_time_dependent=True,
    time_dependent_name="dvdt_ls_timedep",
)

ALL_VARIABLES: List[VariableDefinition] = [
    thls,
    ua,
    va,
    w,
    thetal,
    qt,
    tke,
    ug,
    vg,
    dpdx,
    dpdy,
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
]

ATMO_VARS_BY_NAME: Dict[str, VariableDefinition] = {v.name: v for v in ALL_VARIABLES}


def register_var(var):
    global ALL_VARIABLES, ATMO_VARS_BY_NAME
    ALL_VARIABLES.append(var)
    ATMO_VARS_BY_NAME = {v.name: v for v in ALL_VARIABLES}


def get_var_by_name() -> Dict[str, VariableDefinition]:
    return ATMO_VARS_BY_NAME


def get_all_vars() -> List[VariableDefinition]:
    return ALL_VARIABLES


__all__ = [
    "VariableDefinition",
    "thls",
    "ua",
    "va",
    "w",
    "thetal",
    "qt",
    "tke",
    "ug",
    "vg",
    "dpdx",
    "dpdy",
    "wa",
    "dqtdxls",
    "dqtdyls",
    "tnqt_adv",
    "tnthetal_rad",
    "dudt_ls",
    "dvdt_ls",
    "wtsurf",
    "wqsurf",
    "qtsurf",
    "psurf",
    "qnetav",
    "ALL_VARIABLES",
    "ATMO_VARS_BY_NAME",
]
