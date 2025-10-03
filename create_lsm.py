import matplotlib.pyplot as plt
import numpy as np
import sys, os
# Custom Python scripts/tools/...
from vegetation_properties import ifs_vegetation, top10_to_ifs
from lsm_input_dales import LSM_input_DALES
from landuse_types import lu_types_depac
import do_profiles
import os
import sys
import f90nml
import shutil




# Correction factor for aspect ratio of plots
ASPECT_CORR = 2


def init_dales_grid(domain, lutypes, parnames):
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
    x0 = domain['x0']
    y0 = domain['y0']
    dx = domain['dx']
    dy = domain['dy']
    itot = domain['itot']
    jtot = domain['jtot']
    xsize = itot*dx
    ysize = jtot*dy

    # LES grid in RD coordinates
    x_rd = np.arange(x0+dx/2, x0+xsize, dx)
    y_rd = np.arange(y0+dy/2, y0+ysize, dy)
    # Instance of `LSM_input` class, which defines/writes the DALES LSM input:
    lsm_input = LSM_input_DALES(itot, jtot, 4, lutypes, parnames, debug=False)

    lsm_input.x[:] = x_rd
    lsm_input.y[:] = y_rd

    return lsm_input

def fill_lu_types(lu_types, lsm_input):
    """
    Fills lu_types

    Parameters
    ----------
    lu_types : dict
        properties of each land use type.
    lsm_input : LSM_input_DALES
        Class containing Dales input parameters for all LU types.

    Returns
    -------
    lsm_input : LSM_input_DALES
        Class containing Dales input parameters for all LU types.
    lu_types : dict
        properties of each land use type.
    """

    luname  = [x[1]['lu_long'] for x in lu_types.items()]
    lushort = [x[1]['lu_short'] for x in lu_types.items()]
    lveg    = [x[1]['lveg'] for x in lu_types.items()]
    laqu    = [x[1]['laqu'] for x in lu_types.items()]
    ilu     = np.arange(len(lu_types)) + 1

    setattr(lsm_input, 'luname', luname)
    setattr(lsm_input, 'lushort', lushort)
    setattr(lsm_input, 'lveg', lveg)
    setattr(lsm_input, 'laqu', laqu)
    setattr(lsm_input, 'ilu', ilu)

    # set LU cover for each grid cell
    for lu in lu_types:
        lu_types[lu]['lu_domid'] = np.ones((lsm_input.jtot, lsm_input.itot))*lu_types[lu]['lu_ids'][0]
        lu_types[lu]['lu_frac'] = np.ones((lsm_input.jtot, lsm_input.itot)) * 0
        setattr(lsm_input, 'c_'+lu, lu_types[lu]['lu_frac'])

    return lsm_input, lu_types


def init_lutypes_ifs(lsm_input, lu_dict, lu_types, parnames_lsm ):
    """Assign surface properties to DALES land use types based on ECMWF
       lookup table.

    Parameters
    ----------
    lsm_input : LSM_input_DALES
        Class containing Dales input parameters for all LU types.
    lu_dict : dict
        LU type properties.
    parnames_lsm : list
        List of land use parameters to process

    Returns
    -------
    lsm_input : LSM_input_DALES
        Class containing Dales input parameters for all LU types.

    """
    #
    # Init land use
    #
    shape = (lsm_input.jtot, lsm_input.itot)
    for lu in lu_dict.keys():
        # print('\n processing', lu_dict[lu]['lu_long'])
        for parname in parnames_lsm:
            if parname == 'cover' or parname == 'c_veg':
                parfield = lu_dict[lu]['lu_frac'].copy()
            else:
                # parfield = np.full(shape, np.nan)
                parfield = np.full(shape, 0.0)

            for vt in lu_dict[lu]['lu_ids']:
                iv = top10_to_ifs[vt]     # Index in ECMWF lookup table

                if parname == 'lutype':
                    parfield[:] = iv  # LG: Only apply mask to cover and c_veg (DALES crashes when zeros or nans are in
                                      # the array)
                elif parname == 'tskin':
                    parfield[:] = 273.15  # LG: Only apply mask to cover and c_veg (DALES crashes when zeros or nans
                                          # are in the array)
                elif parname in ["cover","c_veg"]:
                    continue
                else:
                    if parname =='ar':
                        parname_ifs = 'a_r'
                    elif parname =='br':
                        parname_ifs = 'b_r'
                    else:
                        parname_ifs = parname
                    # parfield[mask] = getattr(ifs_vegetation, parname_ifs) [iv]
                    parfield[:] = getattr(ifs_vegetation, parname_ifs) [iv]  # LG: Only apply mask to cover and c_veg
            setattr(lsm_input, '_'.join([parname, lu]), parfield)

    totcover = calc_totcover(lsm_input, lu_types, 'cover')
    setattr(lsm_input, 'cover_tot', totcover)

    totcveg = calc_totcover(lsm_input, lu_types, 'c_veg')
    setattr(lsm_input, 'c_veg_tot', totcveg)

    # TODO: more consistent way to check for LU type with bare soil
    bs_name = [k for k in lu_types.keys() if 'bar' in lu_types[k]['lu_long'].lower()][0]

    cover     = getattr(lsm_input, 'cover_tot')
    cover_bs0 = getattr(lsm_input, 'cover_'+bs_name)
    cover_bs = np.round(1 - cover + cover_bs0, 6)

    setattr(lsm_input, 'cover_' + bs_name, cover_bs)

    #recalculate
    totcover = calc_totcover(lsm_input, lu_types, 'cover')
    setattr(lsm_input, 'cover_tot', totcover)

    totcveg = calc_totcover(lsm_input, lu_types, 'c_veg')
    setattr(lsm_input, 'c_veg_tot', totcveg)

    return lsm_input


def calc_totcover(lsm_input, lu_types, ctype):
    """
    Calculate sum over cover of individual LU types to check if it sums up to 1

    Parameters
    ----------
    lsm_input : LSM_input_DALES
        Class containing Dales input parameters for all LU types.
    lu_types : dict
        LU type properties.
    ctype : str
        LU cover type to be summed.

    Returns
    -------
    totcover : np.array
        Total LU cover.

    """
    covers = [ctype + '_' + s for s in lu_types.keys()]
    totcover = np.zeros([lsm_input.jtot,lsm_input.itot])
    for c in covers:
        totcover+=getattr(lsm_input, c)

    return totcover

def some_plots(lsm_input, plotvars):
    """
    Generate some standard plots of Land Surface Model input data

    Parameters
    ----------
    lsm_input : LSM_input_DALES
        Class containing Dales input parameters for all LU types.
    plotvars : list
        List of variables to plot.

    Returns
    -------
    None.

    """
    import xarray as xr
    data_vars = dict()

    for plotvar in plotvars:
        data = getattr(lsm_input, plotvar)
        data_vars[plotvar] = (['y','x'], data)

    coords = dict(x=(['x'], lsm_input.x),
                  y=(['y'], lsm_input.y))
    ds_lsm = xr.Dataset(data_vars=data_vars,
                        coords=coords
                       )

    for plotvar in list(ds_lsm.variables):
        if plotvar == 'x' or plotvar == 'y': continue
        fig, ax = plt.subplots(1)
        ax.set_aspect(abs((lsm_input.y[-1] - lsm_input.y[0])/(lsm_input.x[-1] - lsm_input.x[0])) * ASPECT_CORR)
        ds_lsm[plotvar].plot(ax=ax, cmap='Reds', vmin=0, vmax=None)
        plt.tight_layout()

    plt.show()

    return


def process_input(lu_types, domain, output_path, exp_id, lplot, modify_func=None):
    """Function that connects all processing steps:
    Init DALES grid
    Write output files in netCDF and binary (to be phased out) format
    Make some standard plots (optional)

    Parameters
    ----------
    lu_types : dict
        LU type properties.
    parnames : list
        List of parameter names to process.
    domain : dict
        Dales domain settings.
    output_path : str
        Dir to write output to.
    start_date : datetime.datetime
        Time stamp of Dales run start.
    exp_id : int
        Experiment ID.
    lplot : bool
        Flag to plot the output.

    Returns
    -------
    lsm_input : LSM_input_DALES
        Class containing Dales input parameters for all LU types.

    """

    # land use parameters
    parnames = ['cover','c_veg','z0m','z0h','lai','ar','br',
                    'lambda_s','lambda_us','rs_min','gD','tskin','lutype']
    lsm_input = init_dales_grid(domain, lu_types, parnames)

    lsm_input.index_soil[:,:] = 2

    # represents average over index_soil=2
    for i in range(4):
        lsm_input.theta_soil[i,:,:] = [0.36867549, 0.25300502, 0.14997292, 0.16459982][i]
        lsm_input.t_soil[i,:,:] = [283.26038083, 286.79894009, 290.88998902, 288.09079126][i]

    # Python -> Fortran indexing
    lsm_input.index_soil += 1

    lsm_input, lu_dict  = fill_lu_types(lu_types, lsm_input)

    if modify_func:
        lu_types, lu_dict, lsm_input = modify_func(lu_types, lu_dict, lsm_input)

    lsm_input = init_lutypes_ifs(lsm_input, lu_dict, lu_types, parnames)

    lsm_input.save_netcdf(f'{output_path}/lsm.inp_{exp_id:03d}.nc')

    if lplot:
        plotvars = ['cover_'+ s for s in lu_types.keys()]
        plotvars.append('cover_tot')
        some_plots(lsm_input, plotvars)

    return lsm_input