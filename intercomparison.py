"""Stiperski IBM intercomparison case"""

import logging
import os
import subprocess
import pytest
import yaml
import numpy as np
import xarray as xr
from modular_dales.Atmosphere import (
    AtmosphereModule,
    AtmosphericProfile,
    InterpolatedProfile,
)
import modular_dales.Atmosphere as Atmosphere
import modular_dales.Configuration as Configuration
import modular_dales.Geometry as Geometry
import modular_dales.Radiation as Radiation
import modular_dales.modular as modular
import modular_dales.logging_wrapper as logging_wrapper
import modular_dales.Surface as Surface
from modular_dales.modular.simulation_module import set_nml_section
import modular_dales.vars as vars

with open("machine_conf.yaml", "r") as f:
    machine_conf = yaml.safe_load(f)
sim = modular.dales_simulation("intercomparison", machine_conf=machine_conf)
dz0 = 10


@np.vectorize
def dz(n):
    if n < 200:
        return dz0
    else:
        return dz0 * 1.023275 ** (n - 200)


grid = Geometry.GridDales(
    itot=512,
    jtot=1024,
    kmax=290,
    xsize=512 * 20,  # 20 m grid spacing
    ysize=1024 * 20,
    x0=0,
    y0=0,
    alpha=1,
    dz0=10,
)
grid.dz = np.zeros(grid.kmax)
grid.zt = np.zeros(grid.kmax)
grid.zm = np.zeros(grid.kmax + 1)
grid.dz[:] = dz(np.arange(grid.kmax) + 1)
grid.zm[1:] = np.cumsum(grid.dz)
grid.zt[:] = 0.5 * (grid.zm[1:] + grid.zm[:-1])
grid.zsize = grid.zm[-1]
print(grid.zt)
print(grid.zm)
print(grid.zsize)

# 10 meter resolution below 2000 m, then stretched grid above that with a stretching factor of 1.023275, which for 290 points gives a grid top at 5 km.
sim += grid

sim += Surface.ConstantFluxesModule(wtsurf=0.12, z0hav=0.1, z0mav=0.1, ps=100000)
set_nml_section(
    sim.nml, sim.nml_docs, "moisture_off", "thermodynamics", "lmoist", False
)
set_nml_section(sim.nml, sim.nml_docs, "theta_randomization", "run", "randthl", 0.25)
set_nml_section(sim.nml, sim.nml_docs, "theta_randomization", "run", "dtmax", 0.2)
