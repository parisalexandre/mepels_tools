#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Concatenate all the xbeach outputs contained in the different directories
into a single one xboutput_global.nc

Author: Alexandre Paris
"""

import os
import glob
import xarray
import numpy as np
from dask.diagnostics import ProgressBar

all_output_files = glob.glob('*/xboutput.nc')
os.system('rm */xboutput_bis.nc') # WARNING! THIS COMMAND IS DANGEROUS !!!!

time = np.zeros(len(all_output_files))
for i, output_file in enumerate(all_output_files):
    cdir = output_file[:-len('xboutput.nc')]
    ds = xarray.open_mfdataset(output_file)
    if i != 0:
        ds = ds.assign_coords(globaltime=(ds.globaltime + time[i-1])) # time adjustment to be the real one
        ds = ds.assign_coords(meantime=(ds.meantime + time[i-1]))
    time[i] = time[i] + float(ds.globaltime[-1]) # last time of previous simulation becomes the first time of the next
    #ds = ds.drop_vars(['k_mean', 'k_var', 'k_min', 'k_max', 'hh_mean', 'hh_var', 'hh_min', 'hh_max',
    #              'Sutot_mean', 'Sutot_var', 'Sutot_min', 'Sutot_max', 'Svtot_mean', 'Svtot_var',
    #              'Svtot_min', 'Svtot_max', 'cctot_mean', 'cctot_var', 'cctot_min', 'cctot_max',
    #              'cx_mean', 'cx_var', 'cx_min', 'cx_max', 'cgx_mean', 'cgx_var', 'cgx_min', 'cgx_max',
    #              'zs_mean', 'zs_var', 'zs_min', 'zs_max', 'H_mean', 'H_var', 'H_min', 'H_max',
    #              'urms_mean', 'urms_var', 'urms_min', 'urms_max', 'As_mean', 'As_var', 'As_min', 'As_max',
    #              'Sk_mean', 'Sk_var', 'Sk_min', 'Sk_max', 'u_mean', 'u_var', 'u_min', 'u_max',
    #              'v_mean', 'v_var', 'v_min', 'v_max', 'ue_mean', 'ue_var', 'ue_min', 'ue_max',
    #              've_mean', 've_var', 've_min', 've_max'])

    # We drop the meantime dimension but there must be a better way
    ds = ds.drop_dims('meantime')
    ds.to_netcdf(path='{}xboutput_bis.nc'.format(cdir), mode='w', format='NETCDF3_64BIT')

all_output_files_bis = glob.glob('*/xboutput_bis.nc')
ds = xarray.open_mfdataset(all_output_files_bis)

ds.to_netcdf(path='./xboutput_global.nc', mode='w', format='NETCDF3_64BIT')

