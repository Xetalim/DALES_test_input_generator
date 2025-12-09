# Creates DALES grid
import numpy as np
from dataclasses import dataclass


class GridDales:
    def __init__(self, input_dic):
        self.input_dic = input_dic
        self.xsize = input_dic["xsize"]
        self.ysize = input_dic["ysize"]

        self.itot = input_dic["itot"]
        self.jtot = input_dic["jtot"]
        self.kmax = input_dic["kmax"]

        self.imax = self.itot
        self.jmax = self.jtot
        self.i1 = self.imax + 1
        self.j1 = self.jmax + 1
        self.k1 = self.kmax + 1
        self.k2 = self.kmax + 2
        self.i2 = self.imax + 2
        self.j2 = self.jmax + 2

        self.x0 = input_dic["x0"]
        self.y0 = input_dic["y0"]

        self.dx = input_dic["xsize"] / self.itot
        self.dy = input_dic["ysize"] / self.jtot

        self.dz0 = input_dic["dz0"]
        self.alpha = input_dic["alpha"]

        i_arr = np.arange(self.itot)
        j_arr = np.arange(self.jtot)

        self.xt = self.x0 + (i_arr + 0.5) * self.dx
        self.yt = self.y0 + (j_arr + 0.5) * self.dy

        self.xm = self.x0 + i_arr * self.dx
        self.ym = self.y0 + j_arr * self.dy

        # self.xt = self.x0 + np.arange(0.5 * self.dx, self.xsize, self.dx)
        # self.yt = self.y0 + np.arange(0.5 * self.dy, self.ysize, self.dy)

        self.xt_openbc = self.xt
        self.yt_openbc = self.yt

        self.xm_openbc = self.x0 + np.arange(0, self.xsize + self.dx, self.dx)
        self.ym_openbc = self.y0 + np.arange(0, self.ysize + self.dy, self.dy)

        # self.xt_ghost = self.x0 + (i_arr + 0.5) * self.dx
        # self.yt_ghost = self.y0 + (j_arr + 0.5) * self.dy

        # self.xm_ghost = self.x0 + i_arr * self.dx
        # self.ym_ghost = self.y0 + j_arr * self.dy
        if self.alpha:
            # Stretched height grid
            self.dz = np.zeros(self.kmax)
            self.zt = np.zeros(self.kmax)
            self.zm = np.zeros(self.kmax + 1)
            self.dz[:] = self.dz0 * (self.alpha) ** np.arange(self.kmax)
            self.zm[1:] = np.cumsum(self.dz)
            self.zt[:] = 0.5 * (self.zm[1:] + self.zm[:-1])
            self.zsize = self.zm[-1]
        else:
            # Equidistant height grid
            self.dz = np.ones(self.kmax) * self.dz0
            self.zsize = self.kmax * self.dz0
            self.zt = np.arange(0.5 * self.dz0, self.zsize, self.dz0)
            self.zm = np.arange(0, self.zsize + self.dz0, self.dz0)
        kmax = self.kmax
        k1 = kmax + 1

        import helper_scripts.LBC.boundary_info as bi

        self.res = bi.openbc_counts_idx(self, fortran_indexing=False)

        print(self.res)
        # zh = np.ones(k1)
        # zf = np.ones(k1)

        # zf[:kmax] = self.zt[:]
        # zh[0] = 0.0
        # for k in range(0, kmax):
        #     zh[k + 1] = zh[k] + 2.0 * (zf[k] - zh[k])
        # zf[k1 - 1] = zf[kmax - 1] + 2.0 * (zh[k1 - 1] - zf[kmax - 1])

        # print(zf)
        # print(zh)

    def as_openbc(self):
        return GridDalesOpenBC(self.input_dic)


class GridDalesOpenBC(GridDales):
    def __init__(self, input_dic):
        super().__init__(input_dic)

        self.xt = self.x0 + np.arange(0.5 * self.dx, self.xsize, self.dx)
        self.yt = self.y0 + np.arange(0.5 * self.dy, self.ysize, self.dy)

        self.xm = self.x0 + np.arange(0, self.xsize + self.dx, self.dx)
        self.ym = self.y0 + np.arange(0, self.ysize + self.dy, self.dy)


def get_domain_info_nml(nml, config):
    return {
        "x0": config["profile"]["grid"]["x0"],
        "y0": config["profile"]["grid"]["y0"],
        "itot": nml["DOMAIN"]["itot"],
        "jtot": nml["DOMAIN"]["jtot"],
        "dx": nml["DOMAIN"]["xsize"] / nml["DOMAIN"]["itot"],
        "dy": nml["DOMAIN"]["ysize"] / nml["DOMAIN"]["jtot"],
        "xsize": nml["DOMAIN"]["xsize"],
        "ysize": nml["DOMAIN"]["ysize"],
        "kmax": nml["DOMAIN"]["kmax"],
        "alpha": config["profile"]["grid"]["alpha_stretch"],
        "dz0": config["profile"]["grid"]["dz"],
        #         stretch: false
        # alpha_stretch: 1.01
        # dz: 5
    }


def generate_dales_domain(domain):
    """
    Initialise a land surface grid with the dimensions of the DALES grid

    Parameters
    ----------
    domain : dict
        disctionary with domain size and resolution
    lutypes : dict
        disctionary with land use types
    parnames : list
        List of parameters to process

    Returns
    -------
    lsm_input : LSM_input_DALES
        Class containing Dales input parameters for all LU types.
    nn_dominant : int
        Number of grid points (+/-) used in "dominant" interpolation method.
    nblockx : int
        Number of blocks in x-direction.
    nblocky : int
        Number of blocks in y-direction.

    """
    # x0, y0 are in RD coordinates
    x0 = domain["x0"]
    y0 = domain["y0"]
    dx = domain["dx"]
    dy = domain["dy"]
    itot = domain["itot"]
    jtot = domain["jtot"]
    xsize = itot * dx
    ysize = jtot * dy

    # LES grid in RD coordinates
    x_rd = np.arange(x0 + dx / 2, x0 + xsize, dx)
    y_rd = np.arange(y0 + dy / 2, y0 + ysize, dy)
    return x_rd, y_rd, itot, jtot

    # ix_west = nesting_idx.ix_west
    # ix_east = nesting_idx.ix_east
    # iy_south = nesting_idx.iy_south
    # iy_north = nesting_idx.iy_north


@dataclass
class nesting_idx:
    ix_west: int
    ix_east: int
    iy_south: int
    iy_north: int
