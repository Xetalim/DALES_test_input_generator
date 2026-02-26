import numpy as np


def lin(z, surf_val, ddz):
    return surf_val + ddz * z


def exp(z, surf_val, lambda_val, z_ml=500):
    u = surf_val * np.exp(-(z - z_ml) / lambda_val)
    u[z <= z_ml] = surf_val
    return u


def linmlsurf(z, lapse_rate, surf_val, offset_val=1.25, z_ml=500):
    u = np.zeros(z.shape)
    u[z <= z_ml] = surf_val - offset_val
    u[z > z_ml] = surf_val - offset_val + (z[z > z_ml] - z_ml) * lapse_rate
    return u


def expsinw(z, surf_val, H, amp, Hp):
    wbase = -surf_val * (1 - np.exp(-z / H))
    wonion = amp * np.sin(2.0 * np.pi / Hp * z)
    wonion[z > Hp] = 0.0
    return wbase + wonion


SHAPE_FUNCTIONS = {
    "lin": lin,
    "expsinw": expsinw,
    "linmlsurf": linmlsurf,
    "exp": exp,
}
