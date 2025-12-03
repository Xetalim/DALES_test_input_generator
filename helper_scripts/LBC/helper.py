import math
import numpy as np
import xarray as xr


def calcBaseprof(zt, thls, ps, pref0=1e5):
    # constants
    lapserate = np.array([-6.5 / 1000.0, 0.0, 1.0 / 1000, 2.8 / 1000])
    zmat = np.array([11000.0, 20000.0, 32000.0, 47000.0])
    grav = 9.81
    rd = 287.04
    cp = 1004.0
    zsurf = 0
    k1 = len(zt)
    # Preallocate
    pmat = np.zeros(4)
    tmat = np.zeros(4)
    rhobf = np.zeros(k1)
    pb = np.zeros(k1)
    tb = np.zeros(k1)
    # DALES code
    tsurf = thls * (ps / pref0) ** (rd / cp)
    pmat[0] = np.exp(
        (
            np.log(ps) * lapserate[0] * rd
            + np.log(tsurf + zsurf * lapserate[0]) * grav
            - np.log(tsurf + zmat[0] * lapserate[0]) * grav
        )
        / (lapserate[0] * rd)
    )
    tmat[0] = tsurf + lapserate[0] * (zmat[0] - zsurf)
    for j in np.arange(1, 4):
        if abs(lapserate[j]) < 1e-10:
            pmat[j] = np.exp(
                (
                    np.log(pmat[j - 1]) * tmat[j - 1] * rd
                    + zmat[j - 1] * grav
                    - zmat[j] * grav
                )
                / (tmat[j - 1] * rd)
            )
        else:
            pmat[j] = np.exp(
                (
                    np.log(pmat[j - 1]) * lapserate[j] * rd
                    + np.log(tmat[j - 1] + zmat[j - 1] * lapserate[j]) * grav
                    - np.log(tmat[j - 1] + zmat[j] * lapserate[j]) * grav
                )
                / (lapserate[j] * rd)
            )
        tmat[j] = tmat[j - 1] + lapserate[j] * (zmat[j] - zmat[j - 1])

    for k in range(k1):
        if zt[k] < zmat[0]:
            pb[k] = np.exp(
                (
                    np.log(ps) * lapserate[0] * rd
                    + np.log(tsurf + zsurf * lapserate[0]) * grav
                    - np.log(tsurf + zt[k] * lapserate[0]) * grav
                )
                / (lapserate[0] * rd)
            )
            tb[k] = tsurf + lapserate[0] * (zt[k] - zsurf)
        else:
            j = 0
            while zt[k] >= zmat[j]:
                j = j + 1
            tb[k] = tmat[j - 1] + lapserate[j] * (zt[k] - zmat[j - 1])
            if abs(lapserate[j]) < 1e-99:
                pb[k] = np.exp(
                    (
                        np.log(pmat[j - 1]) * tmat[j - 1] * rd
                        + zmat[j - 1] * grav
                        - zt[k] * grav
                    )
                    / (tmat[j - 1] * rd)
                )
            else:
                pb[k] = np.exp(
                    (
                        np.log(pmat[j - 1]) * lapserate[j] * rd
                        + np.log(tmat[j - 1] + zmat[j - 1] * lapserate[j]) * grav
                        - np.log(tmat[j - 1] + zt[k] * lapserate[j]) * grav
                    )
                    / (lapserate[j] * rd)
                )
        rhobf[k] = pb[k] / (rd * tb[k])  # dry estimate
    return rhobf


def differentiate(data, coord, order, acc=6):
    ncoef = int(2 * np.floor((order + 1) / 2) - 1 + acc)
    out = xr.zeros_like(data)
    out = out.where(out != 0).load()
    x = data.coords[coord].values
    ipoints = np.arange((ncoef - 1) / 2, len(x) - (ncoef - 1) / 2, dtype=int)
    b = np.zeros((ncoef))
    b[order] = 1.0
    for ip in ipoints:
        A = np.zeros((ncoef, ncoef))
        for j in range(ncoef):
            ip2 = ip - int((ncoef - 1) / 2) + j
            dx = x[ip2] - x[ip]
            for i in range(ncoef):
                A[i, j] = 1 / math.factorial(i) * dx**i
        coef = np.linalg.solve(A, b)
        out[{coord: ip}] = 0.0
        for i in range(ncoef):
            ip2 = ip - int((ncoef - 1) / 2) + i
            out[{coord: ip}] = out[{coord: ip}] + coef[i] * data[{coord: ip2}]
    return out


# @jit(nopython=True, nogil=True)
def interp_z(z, data, z_int):
    data_int = np.zeros(
        (np.shape(data)[0], len(z_int), np.shape(data)[2], np.shape(data)[3])
    )
    # Reverse data if height is descending
    if z[0, 1, 0, 0] < z[0, 0, 0, 0]:
        data = data[:, ::-1, :, :]
        z = z[:, ::-1, :, :]
    for it in range(np.shape(data)[0]):
        for iy in range(np.shape(data)[2]):
            for ix in range(np.shape(data)[3]):
                data_int[it, :, iy, ix] = np.interp(
                    z_int, z[it, :, iy, ix], data[it, :, iy, ix]
                )
    return data_int


def load_data(var, index, drop=False):
    return var.isel(index, drop=drop).values
