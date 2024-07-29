# /usr/bin/python
# -*- coding: utf-8 -*-

"""
Some statistics on wave time series
time Hs Tp Dir Level

Author: Alexandre Paris
"""

import sys
import time
import pathlib
import datetime
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas
import yaml
import xarray as xr
import windrose
from matplotlib import dates


start_time = time.time()


###################################################################################################
def _from_ordinal(serial_dn):
    """
    This function converts a serial date number into a datetime object

    The reference for datetime.fromordinal is 0001/01/01

    Example:
    _from_ordinal(738303)
    >>> 2021-05-07 00:00:00
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

    Example:
    datenum('2009-11-11')
    >>> 734088.0
    """
    date_patterns = ["%Y-%m-%d", "%Y%m%d%H%M", "%Y-%m-%d %H:%M:%S"]
    for pattern in date_patterns:
        try:
            date_d = datetime.datetime.strptime(date, pattern)
            return (366 + date_d.toordinal()
                    + (date_d -
                       datetime.datetime.fromordinal(date_d.toordinal())).total_seconds()
                    / (24 * 60 * 60)
                    )
        except:
            pass

    print('Error: the date pattern is not {}, {} or {}'.format(date_patterns[0],
                                                               date_patterns[1],
                                                               date_patterns[2]))
    sys.exit()


def omega_eq(serie):
    deltat = round(serie['delta_t'].mean())
    num_per_day = round(86400 / deltat)
    phi = int(input('Enter the value of phi calculated by ShoreFor: '))
    ddd = 2*phi*num_per_day
    qqq = 10**(-1/(phi*num_per_day))
    deno = (1-qqq**(ddd+1))/(1-qqq)
    omega_eq = np.zeros(len(serie))
    omega = serie['omega'].values
    vqq = np.ones(((ddd+1), 1))
    for i in range(1, ddd+1):
        vqq[i] = qqq**i
    vqq = np.flipud(vqq)
    for k in range(ddd+1, len(serie)):
        omega_eq[k] = np.matmul(omega[(k-ddd-1):k], vqq)
    serie['omega_eq'] = omega_eq
    serie['omega_eq'] = serie['omega_eq'] / deno
    serie['year'] = serie['date'].apply(_from_ordinal)
    hs98 = serie['Hs'].quantile(q=0.98)
    serie['mask_hs98'] = serie['Hs'].gt(hs98).map({True: r'$Hs > Hs_{98}$',
                                                   False: False})
    serie['Dynamic'] = (serie['omega'] < serie['omega_eq']).map({False: 'erosion',
                                                                 True: 'accretion'})
    serie_eq = serie.copy(deep=True)
    serie_eq.drop(serie_eq[serie_eq['omega_eq'] == 0].index, inplace=True)
    
    return serie_eq, serie


###################################################################################################
# Configuration choice
# Choose study and load data from the configuration file

#beach_letter = input(
#    "Beach you can used :"
#    "\n biscarrosse as 1,"
#    "\n trucvert as 2,"
#    "\n gravelines as 3,"
#    "\n duck as 4,"
#    "\n delilah as 5,"
#    "\n 6: other,"
#    "\n Enter the corresponding number : "
#)
#other = False
#if beach_letter == "1":
#    BEACH_CONFIG = "biscarrosse"
#elif beach_letter == "2":
#    BEACH_CONFIG = "trucvert"
#elif beach_letter == "3":
#    BEACH_CONFIG = "gravelines"
#elif beach_letter == "4":
#    BEACH_CONFIG = "duck"
#elif beach_letter == "5":
#    BEACH_CONFIG = "delilah"
#elif beach_letter == "6":
#    other = True
#    BEACH_CONFIG = input("Enter the file name you want to use:")
#
#if other == True:
#    serie = pandas.read_csv(BEACH_CONFIG)
#    print('File chosen: ', BEACH_CONFIG)
#else:
#    serie = pandas.read_csv("./serie_" + BEACH_CONFIG + ".csv")
#    print("File chosen: ", 'serie_'+BEACH_CONFIG+'.csv')

serie = pandas.read_csv('./stats_waves_duck_daily_8m_1990-2022.csv')
serie.drop(serie[np.isnan(serie['omega'])].index, inplace=True)
sns.set_theme()
sns.set_context("talk")

serie_eq, serie = omega_eq(serie)

## Plot
#fig1, ax1 = plt.subplots(nrows=3, ncols=1)
## Omega timeseries and equilibrium omega
#sns.lineplot(data=serie, x=serie['year'], y=serie['omega'],
#             linewidth=0.5, ax=ax1[1])
#ax1[1].set(xlabel='Year', ylabel='Dean number')
#ax1[1].set_xlim(left=serie_eq['year'].iloc[0], right=serie_eq['year'].iloc[-1])
#ax1[1].tick_params(axis='x', rotation=70)
#sns.lineplot(data=serie, x=serie['year'], y=serie['omega_eq'],
#             linewidth=2, color='red',
#             label='Equilibrium Dean number', ax=ax1[1])
#ax1[1].legend(loc='upper right')
#
## Omega timeseries and dynamic around OmegaEq
#sns.lineplot(data=serie_eq, x=serie['year'], y=serie['omega'],
#             linewidth=0.5, ax=ax1[2],
#             hue='Dynamic')
#ax1[2].set(xlabel='Year', ylabel='Dean number')
#ax1[2].set_xlim(left=serie_eq['year'].iloc[0], right=serie_eq['year'].iloc[-1])
#ax1[2].tick_params(axis='x', rotation=70)
#ax1[2].legend(loc='upper right')
#
## Hs timeseries and Hs > Hs_98
#hue_order = [False, r'$Hs > Hs_{98}$']
#sns.scatterplot(data=serie, x=serie['year'], y=serie['Hs'],
#                linewidth=0.5, ax=ax1[0],
#                hue='mask_hs98', s=3, hue_order=hue_order,
#                edgecolor='none')
#handles, labels = ax1[0].get_legend_handles_labels()
#ax1[0].legend(handles=handles[1:], labels=labels[1:], loc='upper right')
#ax1[0].set(xlabel='Year', ylabel='Hs (m)')
#ax1[0].set_xlim(left=serie_eq['year'].iloc[0], right=serie_eq['year'].iloc[-1])
#ax1[0].tick_params(axis='x', rotation=70)

####################################################################################################
#events = serie_eq[['date', 'year', 'omega', 'state', 'omega_eq', 'Dynamic']].copy()
#events = events.reset_index()
#events['annee'] = events['year'].dt.year
##events.drop(events[events['annee'] != 2019].index, inplace=True)
#events['year'] = events['date'].apply(_from_ordinal)
#events = events.reset_index()
#events['mask'] = events['omega'] > events['omega_eq']
#events['total_duration'] = pandas.to_timedelta(events['year']
#                                               - events['year'][0]).astype('timedelta64[s]')
#events['event_duration'] = events['total_duration']
#events = events.assign(k=[1]*len(events.index))
#for i in range(1, len(events.index)):
#    if events.loc[i, 'mask'] == events.loc[i-1, 'mask']:
#        events.loc[i, 'k'] = events.loc[i-1, 'k'] + 1
#    events.loc[i-1,
#               'event_duration'] = events.loc[i-1,
#                                              'total_duration'] - events.loc[i-events.loc[i-1, 'k'],
#                                                                             'total_duration']
#
#events.iloc[-1, events.columns.get_loc('event_duration')] = events.iloc[-1, events.columns.get_loc('total_duration')] -\
#                                                          events.iloc[-events.iloc[-1, events.columns.get_loc('k')],
#                                                                     events.columns.get_loc('total_duration')]
#events['event_duration'] = events['event_duration'] / 3600 / 24  # Converts into days
#events['real_duration'] = events['event_duration']
#
#for i in range(len(events.index)-1, 1, -1):  # question du zero non traitée
#    if events.loc[i, 'event_duration'] > events.loc[i-1, 'event_duration']:
#        events.loc[i-1, 'real_duration'] = events.loc[i, 'real_duration']
#
## Keep the state duration of the last element for each cluster
#events['real_duration'] = events['real_duration'].mask(events['real_duration'].shift(-1) == events['real_duration'])
#
#events['duration_levels'] = pandas.cut(events['real_duration'], [0, 4, 8, 50],
#                                       labels=['less than 4 days', 'between 4 and 8 days', 'more than 8 days'])
#
##print(events)

##########
# Erosion, accretion or equilibrium
lim = 0.5
conditions = [
        (serie['omega'] <= serie['omega_eq'] - lim),
        (serie['omega'] > serie['omega_eq'] - lim) & (serie['omega'] <= serie['omega_eq'] + lim),
        (serie['omega'] > serie['omega_eq'] + lim)
        ]

morphodyn = ['accretion', 'equilibrium', 'erosion']
serie['morphodyn'] = np.select(conditions, morphodyn)

print(serie.columns)
#serie.drop(serie[serie['year'].dt.year != 2008].index, inplace=True)
serie.drop(serie[serie['year'].dt.year != 2019].index, inplace=True)
#serie_first = serie[serie['year'] <= '2008-08-31']#.index, inplace=True
#serie_second = serie[serie['year'] >= '2008-08-31']

##########
# Plot
fig9, ax9 = plt.subplots(2)
sns.scatterplot(data=serie, x='year', y='omega', hue='state', style='state', ax=ax9[0], s=50, #s=200,
                hue_order=['reflective', 'intermediate', 'dissipative'],
                style_order=['reflective', 'intermediate', 'dissipative'])
legend_labels, _ = ax9[0].get_legend_handles_labels()
ax9[0].legend(legend_labels, ['reflective',
                              'intermediate',
                              'dissipative'], bbox_to_anchor=(1, 1))

sns.scatterplot(data=serie, x='year', y='omega', hue='morphodyn', style='morphodyn', ax=ax9[1], s=50, # s=200,
                hue_order=['accretion', 'equilibrium', 'erosion'],
                style_order=['accretion', 'equilibrium', 'erosion'])
legend_labels, _ = ax9[1].get_legend_handles_labels()
ax9[1].legend(legend_labels, ['accretion',
                              'equilibrium',
                              'erosion'], loc='best', fontsize='xx-small')#, bbox_to_anchor=(1, 1))

#ax9[0].set(xlabel='Date', ylabel='Dean number')
ax9[0].set(ylabel='Dean number')
ax9[1].set(xlabel='Date', ylabel='Dean number')
#ax9[1].set(xlabel='Month', ylabel='Dean number')
ax9[0].plot(serie['year'], serie['omega_eq'], linewidth=4, linestyle='--', color='red')
ax9[1].plot(serie['year'], serie['omega_eq'], linewidth=4, linestyle='--', color='red')
#first = ['2019-02-01', '2019-02-15', '2019-03-25', '2019-04-11', '2019-04-17', '2019-05-20', '2019-06-25', '2019-07-16', '2019-08-15']
#second = ['2019-09-03', '2019-09-09', '2019-09-17', '2019-09-23', '2019-10-15', '2019-10-25', '2019-11-14', '2019-11-22', '2019-12-06']
#for date in first:
#    ax9[0].axvline(serie.loc[serie['year'] == date, 'year'].values, color='magenta')
#    ax9[1].axvline(serie.loc[serie['year'] == date, 'year'].values, color='magenta')
#for date in second:
#    ax9[0].axvline(serie.loc[serie['year'] == date, 'year'].values, color='magenta')
#    ax9[1].axvline(serie.loc[serie['year'] == date, 'year'].values, color='magenta')

#ax9[0].xaxis.set_major_locator(dates.DayLocator(interval=10))
#ax9[0].xaxis.set_major_formatter(dates.DateFormatter("%Y-%m-%d"))
ax9[0].tick_params(axis='x', rotation=90)#, labelsize=20)
ax9[1].xaxis.set_major_locator(dates.DayLocator(interval=30))
ax9[1].xaxis.set_major_formatter(dates.DateFormatter("%m-%d"))
ax9[1].tick_params(axis='x', rotation=90)#, labelsize=20)
#ax9[0].set_xticks([])
ax9[0].set_xlabel([])


plt.show()

#serie = serie.drop(columns='index')
serie = serie.set_index('date')
#serie.to_csv('./dynamic_omegas_duck_2019_daily.csv', sep=' ', header=False, columns=['year', 'morphodyn'])
serie.to_csv('./dynamic_omegas_duck_2019_daily_phibar.csv', sep=' ', header=False, columns=['year', 'morphodyn'])

print('\nEnd of analysis')
print('--- {} seconds ---'.format(time.time() - start_time))
