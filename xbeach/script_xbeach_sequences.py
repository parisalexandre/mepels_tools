#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Launch a series of xbeach simulations following the dynamic sequences
contained in file_dynamic. This file should contain three columns:
- date in days
- date in yyyy-mm-dd
- dynamic: erosion, accretion, or equilibrium

Author: Alexandre Paris
"""

import os
import sys
import subprocess
import time
import pandas
from scipy.io import netcdf
import numpy as np


file_dynamic = './dynamic_omegas_duck_2019_daily.csv'

###############################################################################
# Change the following values with your parameters
facua_accretion = 0.2
waveform_accretion = 'vanthiel'
facua_erosion = 0.15
waveform_erosion = 'ruessink_vanrijn'
facua_equi = 0.15
waveform_equi = 'vanthiel'

def facua_waveform(dyn):
    if dyn == 'accretion':
        facua = facua_accretion
        waveform = waveform_accretion
    elif dyn == 'erosion':
        facua = facua_erosion
        waveform = waveform_erosion
    elif dyn == 'equilibrium':
        facua = facua_equi
        waveform = waveform_equi
    return facua, waveform


###############################################################################
# Read the dynamic file
serie = pandas.read_csv(file_dynamic, sep=' ', names=['date', 'year', 'dynamic']).values

tide = np.loadtxt('tide_2019_entire.txt')
wave = np.loadtxt('wave_2019_entire.txt')

date_initial = serie[0, 0]
date_initial_yformat = serie[0, 1]

for i in range(1, len(serie)):
    if serie[i, 2] != serie[i-1, 2]:
        date_final_yformat = serie[i-1, 1]
        date_final = serie[i-1, 0]
        print('initial date: ', date_initial_yformat, 'final date: ', date_final_yformat)
        facua, waveform = facua_waveform(serie[i-1, 2])
        tstart = (date_initial - 737426)*24*3600 # 737426 corresponds to 2019-01-01 and to the zero of tide and wave files
        tstop = (date_final - 737426)*24*3600
        diff = tstop - tstart
        if diff == 0:
            tstop = tstop + 24*3600

        os.system('mkdir {}_{}'.format(date_initial_yformat, date_final_yformat))
        os.system('cp params.txt launch_xbeach.oar ./{}_{}'.format(date_initial_yformat, date_final_yformat))
        os.chdir('./{}_{}'.format(date_initial_yformat, date_final_yformat))
        
        ##########
        # Replacement of values in params.txt and launch_xbeach.oar

        with open('params.txt', 'r') as f:
            params = f.read()
            params = params.replace('facua        = 0.05', 'facua        = {}'.format(facua))
            params = params.replace('waveform = vanthiel', 'waveform = {}'.format(waveform))
            if date_initial_yformat != '2019-02-01':
                params = params.replace('depfile   = ../bathy_2019-02-01.dep',
                                        'depfile   = ../bathy_{}.dep'.format(date_initial_yformat))
                #params = params.replace('tstart       = 2678400', 'tstart       = {}'.format(tstart))
            params = params.replace('tstop     = 0', 'tstop     = {}'.format(tstop - tstart))
        with open('params.txt', 'w') as f:
            f.write(params)
        
        with open('launch_xbeach.oar', 'r') as f:
            launch = f.read()
            launch = launch.replace('#OAR -n  sequence_xbeach',
                                    '#OAR -n  sequence_{}_{}'.format(date_initial_yformat, date_final_yformat))
        with open('launch_xbeach.oar', 'w') as f:
            f.write(launch)

        ##########
        # Modification of tide and wave files
        new_tide = tide[tide[:, 0] >= tstart]
        new_tide = new_tide[new_tide[:, 0] <= tstop]
        new_tide[:, 0] = new_tide[:, 0] - tstart
        np.savetxt('tide_{}_{}.txt'.format(date_initial_yformat, date_final_yformat), new_tide, fmt='%1.16g')

        new_wave = wave[int(tstart/3600 - 1):int(tstop/3600)]
        np.savetxt('wave_{}_{}.txt'.format(date_initial_yformat, date_final_yformat), new_wave, fmt='%1.16g')

        with open('params.txt', 'r') as f:
            params = f.read()
            params = params.replace('zs0file   = ../tide_2019_entire.txt',
                                    'zs0file   = ./tide_{}_{}.txt'.format(date_initial_yformat, date_final_yformat))
            params = params.replace('bcfile     = ../wave_2019_entire.txt',
                                    'bcfile     = ./wave_{}_{}.txt'.format(date_initial_yformat, date_final_yformat))

        with open('params.txt', 'w') as f:
            f.write(params)

        ##########
        # XBeach launch
        fail = 1
        while fail == 1:
            os.system('oarsub -S ./launch_xbeach.oar')
        
            while not os.path.exists('./XBlog.txt'):
                time.sleep(2)

            os.system('oarstat -u > oarstatfile')
            while os.stat('oarstatfile').st_size != 0:
                time.sleep(60)
                os.system('oarstat -u > oarstatfile')
        
            ##########
            # Check if the job is done
            with open('XBlog.txt', 'r') as f:
                last_line = f.readlines()[-1]

            if last_line != '  End of program xbeach\n':
                #print('ERROR: Check the log files and investigate')
                #sys.exit()
                for filename in os.listdir(os.getcwd()):
                    if filename[-6:] == 'stderr':
                        with open(filename, 'r') as f:
                            if 'failed' in f.read():
                                os.system('rm OAR* XB* E* e* q* xboutput.nc')
                                fail = 1
                                print('Something went wrong. The job is submitted again')
            else:
                fail = 0

        ##########
        # The last date becomes the initial date for the next simulation
        if diff == 0:
            date_initial = date_final + 1
            date_initial_yformat = date_final_yformat[:9] + '{}'.format(int(date_final_yformat[9]) + 1)
        else:
            date_initial = date_final
            date_initial_yformat = date_final_yformat
            
        ##########
        # Extraction of final bathymetry
        file2read = netcdf.NetCDFFile('./xboutput.nc', 'r')
        zb_s = file2read.variables['zb'][:].copy()
        zb = np.squeeze(zb_s)
        file2read.close()
        with open('../bathy_{}.dep'.format(date_initial_yformat), 'w') as f: # new first date
            f.write('-8 ') # We rewrite the first value removed by Xbeach
            for z in list(zb[-1, 1:-1]):
                f.write(str(z) + ' ')
            f.write('2.41840626') # We rewrite the last too

        ##########
        os.chdir('../')


