#!/usr/bin/env python


import numpy as np
from scipy.interpolate import interp1d
import os
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
# linear profile => Use for u
def lin(z, surf_val, ddz):
    return surf_val + ddz*z

# With mixed layer, where u0 is set at zml => Use for qt
def exp(z, surf_val, lambda_val, z_ml=500):
    u = surf_val * np.exp(-(z-z_ml) / lambda_val)
    u[z<=z_ml] = surf_val
    return u

# With mixed layer with fixed offset (du0) from surface value (u0),
# and sloping FT (dudz) from fixed mixed layer height (zml)
# Use for thl
def linmlsurf(z, lapse_rate, surf_val, offset_val=1.25, z_ml=500):
    u = np.zeros(z.shape)
    u[z<=z_ml] = surf_val - offset_val # Positive offsets are reductions w.r.t surface
    u[z>z_ml] = surf_val - offset_val + (z[z>z_ml] - z_ml)*lapse_rate
    return u

# Sum of exponential dropoff and sine => Use for subsidence
def expsinw(z, surf_val, H, amp, Hp):
    wbase = -surf_val*(1-np.exp(-z/H))
    wonion = amp*np.sin(2.*np.pi/Hp*z)
    wonion[z>Hp] = 0.
    return wbase + wonion

def output_profiles(profile_config, nml, out_dir, plot=False):
    exp_id = nml["RUN"]["iexpnr"]
    case = "case"


    Nz = nml['DOMAIN']['kmax']
    stretch = profile_config["grid"]["stretch"]
    alpha = profile_config["grid"]["alpha_stretch"]
    dz = profile_config["grid"]["dz"]


    if not stretch:
        # equally spaced levels, first level is 1/2 dz over ground
        z = (np.arange(Nz) +.5) * dz
    else:
        # stretched vertical levels. First level is of size dz and located at dz/2, dz then stretched by factor alpha
        z = np.zeros(Nz)
        z[0] = dz/2
        dz_ = dz
        for i in range(1,Nz):
            z[i] = z[i-1] + dz_
            dz_ *= alpha

    #total height:
    if stretch and alpha > 1:
        Ztot = dz * (1-alpha**Nz) / (1-alpha)
    else:
        Ztot = dz*Nz
    var_dic = {"u":None,"v":None,"w":None,"thl":None,"qt":None,"tke":None,"thl_tend":None,"qt_tend":None}
    shape_function_dic = {"lin":lin,"expsinw":expsinw,"linmlsurf":linmlsurf,"exp":exp}
    for var in var_dic.keys():
        if var in profile_config["profiles"]:
            var_config = profile_config["profiles"][var]
            shape_function = shape_function_dic[var_config["shape"]]
            shape_config = var_config["params"]
            var_dic[var] = shape_function(z, **shape_config)
        elif var in profile_config["interpolated"]:
            var_config = profile_config["interpolated"][var]
            z_points = var_config["z"]
            var_points = var_config["points"]
            fill_value = var_config["fill_value"]
            var_dic[var] =  interp1d(z_points,  var_points, fill_value=fill_value)(z)


    with netCDF4.Dataset(out_dir/ f"init.{exp_id:03d}.nc", "w", format="NETCDF3_CLASSIC") as ncout:
        # create time dimension
        nc_z_dim = ncout.createDimension('zh', Nz)

        # normal profiles
        nc_z_var = ncout.createVariable('zh', 'f8', ('zh',))
        nc_thetal_var = ncout.createVariable('thetal', 'f8', ('zh',))
        nc_qt_var = ncout.createVariable('qt', 'f8', ('zh',))
        nc_ua_var = ncout.createVariable('ua', 'f8', ('zh',))
        nc_va_var = ncout.createVariable('va', 'f8', ('zh',))
        nc_tke_var = ncout.createVariable('tke', 'f8', ('zh',))

        # geostrophic
        nc_ug_var = ncout.createVariable('ug', 'f8', ('zh',))
        nc_vg_var = ncout.createVariable('vg', 'f8', ('zh',))
        nc_w_subs_var = ncout.createVariable('wa', 'f8', ('zh',))
        nc_qt_tend_var = ncout.createVariable('tnqt_adv', 'f8', ('zh',))
        nc_dthlrad_var = ncout.createVariable('tnthetal_rad', 'f8', ('zh',))

        nc_z_var[:] = z
        nc_thetal_var[:] = var_dic["thl"]
        nc_qt_var[:] = var_dic["qt"]
        nc_ua_var[:] = var_dic["u"]
        nc_va_var[:] = var_dic["v"]
        nc_tke_var[:] = var_dic["tke"]

        nc_ug_var[:] = var_dic["u"]
        nc_vg_var[:] = var_dic["v"]
        nc_w_subs_var[:] = var_dic["w"]
        nc_qt_tend_var[:] = var_dic["qt_tend"]
        nc_dthlrad_var[:] = var_dic["thl_tend"]
        
        
    if plot:
        os.makedirs(out_dir / ".." / "profiles", exist_ok=True)
        for var, varname in zip(["thl", "qt", "u", "v", "tke", "w", "qt_tend", "thl_tend"], ["Theta_l","q_t","u","v","tke","w_subsidence","Tendency_q_t", "Tendency_theta_l_radiation"]):
            fig, ax = plt.subplots()
            ax.plot(var_dic[var], z)
            ax.set_xlabel(varname)
            ax.set_ylabel("z (m)")
            ax.set_title(f"exp {exp_id:03d} {varname}")
            fig.savefig(out_dir / ".." / "profiles" / f"profile_{varname}.png", dpi=300)
            plt.close()