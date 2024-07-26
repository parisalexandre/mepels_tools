# /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Read NetCDF wave files and create a .dat file with wave data
time Hs Tp Dir Level
"""

import os
import subprocess
import datetime
import netCDF4 as nc
import csv
import pandas
import xarray as xr


def _from_ordinal(serial_dn):
    """
    This function converts a serial date number into a datetime object

    The reference for datetime.fromordinal is 0001/01/01
    """
    i_x = int(serial_dn)
    classic_date = datetime.datetime.fromordinal(i_x - 366)
    remainder = float(serial_dn) - i_x
    hour, remainder = divmod(24 * remainder, 1)
    minute, remainder = divmod(60 * remainder, 1)
    second, remainder = divmod(60 * remainder, 1)
    microsecond = int(1e6 * remainder)
    if microsecond < 10:
        microsecond = 0  # compensate for rounding errors
    classic_date = datetime.datetime(
        classic_date.year,
        classic_date.month,
        classic_date.day,
        int(hour),
        int(minute),
        int(second),
        microsecond,
    )

    if microsecond > 999990:  # compensate for rounding errors
        classic_date += datetime.timedelta(microseconds=1e6 - microsecond)

    return classic_date


def datenum(date):
    """
    Equivalent of datenum Matlab function
    """
    date_patterns = ["%Y-%m-%d", "%Y%m%d%H%M", '%Y-%m-%d %H:%M:%S']
    for pattern in date_patterns:
        try:
            d = datetime.datetime.strptime(date, pattern)
        except:
            pass

    return (
        366
        + d.toordinal()
        + (d - datetime.datetime.fromordinal(d.toordinal())).total_seconds()
        / (24 * 60 * 60)
    )


def read_netcdf(nc_file, year, month):
    """
    Reads each FRF file and creates a temp file with wave data
    """
    nc = xr.open_dataset(nc_file, drop_variables=['station_name', 'nominalDepth', 'depth', 'gaugeDepth',
                                                  'waveTm1', 'waveTm2', 'wavePeakDirectionPeakFrequency',
                                                  'waveMeanDirectionPeakFrequency', 'wavePrincipleDirection',
                                                  'waveFrequency', 'waveDirectionBins', 'waveEnergyDensity',
                                                  'directionalWaveEnergyDensity', 'qcFlagE', 'qcFlagD',
                                                  'directionalPeakSpread', 'spectralWidthParameter',
                                                  'waveDirectionEstimator', 'waveA1Value', 'waveB1Value',
                                                  'waveA2Value', 'waveB2Value', 'airPressure', 'sourceZ'])
    data = nc.to_dataframe()
    # Simple Green law
    if '17m' in nc_file:
        print(nc_file)
        data['waveHs'] = data['waveHs'] * (17/8)**(1/4)
    if '26m' in nc_file:
        print(nc_file)
        data['waveHs'] = data['waveHs'] * (26/8)**(1/4)

    data = data.dropna()
    #data_daily = data.resample('H').mean()
    data_daily = data
    data_daily['date'] = data_daily.index.round('S')
    data_daily['date'] = data_daily['date'].apply(str)
    data_daily['date'] = data_daily['date'].apply(datenum)
    data_daily = data_daily.dropna()
    
    if month < 10:
        if 'waterLevel' in data_daily:
            data_daily.to_csv('temp_{}_0{}.dat'.format(year, month), index=False, header=False,
                              sep=' ', columns=['date', 'waveHs', 'waveTp', 'waveMeanDirection', 'waterLevel'])#,
       #                       na_rep=float('NaN'))
        else:
            data_daily['waterLevel'] = float('NaN')
            data_daily.to_csv('temp_{}_0{}.dat'.format(year, month), index=False, header=False,
                              sep=' ', columns=['date', 'waveHs', 'waveTp', 'waveMeanDirection', 'waterLevel'])#,
       #                       na_rep=float('NaN'))
    else:
        if 'waterLevel' in data_daily:
            data_daily.to_csv('temp_{}_{}.dat'.format(year, month), index=False, header=False,
                              sep=' ', columns=['date', 'waveHs', 'waveTp', 'waveMeanDirection', 'waterLevel'])#,
       #                       na_rep=float('NaN'))
        else:
            data_daily['waterLevel'] = float('NaN')
            data_daily.to_csv('temp_{}_{}.dat'.format(year, month), index=False, header=False,
                              sep=' ', columns=['date', 'waveHs', 'waveTp', 'waveMeanDirection', 'waterLevel'])#,
       #                       na_rep=float('NaN'))


# Year by year and month by month
for year in range(1990, 2023):  # Change dates if needed
    for month in range(1, 13):
        if month < 10:
            try:
                nc_file = 'FRF-ocean_waves_8m-array_{}0{}.nc'.format(year, month)
                read_netcdf(nc_file, year, month)
            except FileNotFoundError:
                pass
            #if os.path.isfile('FRF-ocean_waves_8m-array_{}0{}.nc'.format(year, month)) is False:
            #    if os.path.isfile('../waves_17m/FRF-ocean_waves_waverider-17m_{}0{}.nc'.format(year, month)) is False:
            #        nc_file = '../waves_26m/FRF-ocean_waves_waverider-26m_{}0{}.nc'.format(year, month)
            #    else:
            #        nc_file = '../waves_17m/FRF-ocean_waves_waverider-17m_{}0{}.nc'.format(year, month)
            #else:
            #    nc_file = 'FRF-ocean_waves_8m-array_{}0{}.nc'.format(year, month)
        else:
            #if os.path.isfile('FRF-ocean_waves_8m-array_{}{}.nc'.format(year, month)) is False:
            #    if os.path.isfile('../waves_17m/FRF-ocean_waves_waverider-17m_{}{}.nc'.format(year, month)) is False:
            #        nc_file = '../waves_26m/FRF-ocean_waves_waverider-26m_{}{}.nc'.format(year, month)
            #    else:
            #        nc_file = '../waves_17m/FRF-ocean_waves_waverider-17m_{}{}.nc'.format(year, month)
            #else:
            #    nc_file = 'FRF-ocean_waves_8m-array_{}{}.nc'.format(year, month)
            try:
                nc_file = 'FRF-ocean_waves_8m-array_{}{}.nc'.format(year, month)
                read_netcdf(nc_file, year, month)
            except FileNotFoundError:
                pass
        #read_netcdf(nc_file, year, month)

# Cleaning output file
os.system('cat temp_* > temp.dat')
out_file = open('2_temp.dat', "w")
subprocess.call(['sed', '-e', 's/\"1/1/g', '-e', 's/\"2/2/g', '-e', 's/\" / /g','temp.dat'], stdout=out_file)
out_file.close()
os.system('rm temp*')

data = pandas.read_csv('2_temp.dat', sep=' ', names=['date', 'waveHs', 'waveTp', 'waveMeanDirection', 'waterLevel'])#, index_col='date')
os.system('rm 2_temp.dat')

data['time'] = data['date'].apply(_from_ordinal)
data = data.set_index('time')
data = data.drop(columns='date')
#data = data.resample('H').mean()
#data = data.resample('1D').mean()
#data = data.interpolate()
data = data.reset_index()
#data['date'] = data['time'].apply(str)
#data['date'] = data['date'].apply(datenum)
#data = data.set_index('date')
#data = data.set_index('time')
#data = data.dropna()
data.to_csv('forcing_duck_daily_8m_1990-2022.dat', header=False, sep=' ', na_rep=float('NaN'), index=False,
            columns=['time', 'waveHs', 'waveTp', 'waveMeanDirection', 'waterLevel'])

