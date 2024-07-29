# /usr/bin/python
# -*- coding: utf-8 -*-

"""
Some statistics on wave time series
time Hs Tp Dir Level

WARNING: search 'to_csv' and change the name of the file for your needs
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

start_time = time.time()

##############################################################################

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


def masselink_plot(ax):
    """
       Background for Figure 4 of Masselink et al (1993)
       to plot the repartition of Omegas and beach states
    """
    from matplotlib.ticker import LogLocator, FixedLocator, FixedFormatter
    ax.vlines(x=2, ymin=0, ymax=15, color='gray', ls='--')
    ax.vlines(x=5, ymin=0, ymax=7, color='gray', ls='--')
    ax.hlines(y=3, xmin=0, xmax=8, color='gray', ls='--')
    ax.hlines(y=7, xmin=0, xmax=8, color='gray', ls='--')
    ax.hlines(y=15, xmin=0, xmax=8, color='gray', ls='--')
    ax.set_xticks(ticks=[0, 2, 5, 8])
    #ax.set_yticks(ticks=[0, 3, 7, 15])
    ax.set_xlim(0, 8)
    ax.set_yscale('log')
    ax.invert_yaxis()
    ax.set_ylim(16, 0)
    ticks = [0.1, 3, 5, 15]  # Using 0.1 instead of 0 because log scale can't have zero
    ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.yaxis.set_major_formatter(FixedFormatter(['0', '3', '5', '15']))
    ax.xaxis.tick_top()
    ax.set_xlabel('DIMENSIONLESS FALL VELOCITY OMEGA')
    ax.xaxis.set_label_position('top')
    ax.set_ylabel('RELATIVE TIDE RANGE RTR = MSR/Hb')
    ax.text(1, -0.1, 'REFLECTIVE', horizontalalignment='center')
    ax.text(3.5, -0.1, 'INTERMEDIATE', horizontalalignment='center')
    ax.text(6.5, -0.1, 'DISSIPATIVE', horizontalalignment='center')
    ax.text(1, 1.5, 'reflective', horizontalalignment='center')
    ax.text(3.5, 1.5, 'barred', horizontalalignment='center')
    ax.text(6.5, 1.5, 'barred dissipative', horizontalalignment='center')
    ax.text(1, 5, 'low tide terrace + rip', horizontalalignment='center')
    ax.text(3.5, 5, 'low tide bar/rip', horizontalalignment='center')
    ax.text(6.5, 5, 'non-barred dissipative', horizontalalignment='center')
    ax.text(1, 11, 'low tide terrace', horizontalalignment='center')
    ax.text(5, 11, 'ultra dissipative', horizontalalignment='center')
    ax.text(4, 15.5, 'transition to tide-dominated tidal flats',
            horizontalalignment='center', verticalalignment='center')


def calculate_return(df, colname):
    """
        Calculate the return period and exceedance probability of a variable
    """
    sorted_data = df.sort_values(by=colname)
    n = sorted_data.shape[0]
    sorted_data.insert(0, 'rank', range(1, 1+n))
    sorted_data['exceedance_probability'] = ((n - sorted_data['rank'] + 1) / (n))
    sorted_data['return_period'] = (1 / sorted_data['exceedance_probability'])
    return sorted_data


def hilbert(imf, sample_rate):
    """
        Calculate the variables needed to plot the Hilbert transform
        based on the imfs of a variable
    """
    _, IF, IA = emd.spectra.frequency_transform(imf, sample_rate, 'nht')
    freq_edges, freq_centres = emd.spectra.define_hist_bins(0, 1e-5, 5000, 'linear')
    _, spec_weighted = emd.spectra.hilberthuang(IF, IA, freq_edges, sum_imfs=False)
    _, spec_unweighted = emd.spectra.hilberthuang(IF, np.ones_like(IA), freq_edges, sum_imfs=False)
    return freq_centres, spec_unweighted, spec_weighted


def search_and_replace_NaN(var):
    """
        Find NaN and replace by 0
    """
    count = serie[var].isna().sum()
    serie[var] = serie[var].fillna(0)
    if count > 0:
        if count == 1:
            print('WARNING: {} value of {} was NaN and has been replaced by 0'.format(count, var))
        else:
            print('WARNING: {} values of {} were NaN and have been replaced by 0'.format(count,
                                                                                         var))


###################################################################################################
# Configuration choice
# Choose study and load data from the configuration file

BEACH_CONFIG = input("Enter the yaml file name (without _configuration.yaml):")
with open(r"./" + BEACH_CONFIG + "_configuration.yaml") as file:
    cl = yaml.full_load(file)
print("beach config choose: ", BEACH_CONFIG)

sns.set_theme()
sns.set_context("talk")

plt.rcParams.update({'font.size': 18,#26,
                     'axes.titlesize': 18,#26,
                     'axes.labelsize': 18,#26,
                     'xtick.labelsize': 18,#26,
                     'ytick.labelsize': 18,#26,
                     'legend.fontsize': 16})#24})

wave_file = cl["waves_data_path"]["input_waves"]
#wave_file = './forcing_duck_daily_8m_1990-2022.dat'
extension = pathlib.Path(wave_file).suffix

# Read data
if extension == '.nc':
    nc = xr.open_dataset(wave_file)
    serie = nc.to_dataframe()
elif extension in ('.txt', '.dat', '.csv'):  # TESTER CSV
    serie = pandas.read_csv(wave_file, sep='\s+')# delim_whitespace=True)
    if len(serie.columns) == 3:
        serie.columns = ['date', 'Hs', 'Tp']
    elif len(serie.columns) == 4:
        serie.columns = ['date', 'Hs', 'Tp', 'Dir']
    elif len(serie.columns) == 5:
        serie.columns = ['date', 'Hs', 'Tp', 'Dir', 'Level']

search_and_replace_NaN('Hs')
search_and_replace_NaN('Tp')
serie['Hs'] = np.where(serie['Hs'] < 0, 0, serie['Hs'])

serie['year'] = serie['date'].apply(_from_ordinal)

serie_prob_Hs = calculate_return(serie, 'Hs')
serie_prob_Hs['return_period'] = serie_prob_Hs['return_period'] / 24   # return period in days
serie_prob_Hs['exceedance_probability'] = serie_prob_Hs['exceedance_probability']  # * 100
H12_y = serie_prob_Hs.loc[(serie_prob_Hs['exceedance_probability'] >= 0.0012) &
                          (serie_prob_Hs['exceedance_probability'] <= 0.0014), 'Hs'].iloc[0]

serie_prob_Tp = calculate_return(serie, 'Tp')
serie_prob_Tp['return_period'] = serie_prob_Tp['return_period'] / 24  # return period in days
serie_prob_Tp['exceedance_probability'] = serie_prob_Tp['exceedance_probability']  # * 100
T12_y = serie_prob_Tp.loc[(serie_prob_Tp['exceedance_probability'] >= 0.0012) &
                          (serie_prob_Tp['exceedance_probability'] <= 0.0014), 'Tp'].iloc[0]


H_mean = serie['Hs'].mean()
T_mean = serie['Tp'].mean()

print('H_mean=', H_mean, 'm; H12_y=', H12_y, 'm')
print('T_mean=', T_mean, 's; T12_y=', T12_y, 's')

with open(r"./" + BEACH_CONFIG + "_configuration.yaml", 'r') as file:
    get_all = file.readlines()

with open(r"./" + BEACH_CONFIG + "_configuration.yaml", 'w') as file:
    for i, line in enumerate(get_all, 1):
        if 'T12_Y' in line:
            file.writelines('    T12_Y: {}\n'.format(T12_y))
        elif 'H12_Y' in line:
            file.writelines('    H12_Y: {}\n'.format(H12_y))
        elif 'H_mean' in line:
            file.writelines('    H_mean: {}\n'.format(H_mean))
        elif 'T_mean' in line:
            file.writelines('    T_mean: {}\n'.format(T_mean))
        else:
            file.writelines(line)

#############################PLOT SERIES HS AND TP#############################
# Figure 1

rep = 'nan'
while rep != 'y' or rep != 'n':
    rep = input('Figure 1: Plot timeseries of Hs and Tp? y/n \t')
    if rep in ('y', 'n'):
        break

if rep == 'y':
    fig1, ax1 = plt.subplots(nrows=2, ncols=1)
    ax1[0].plot(serie['year'], serie['Hs'], linewidth=1, color='black')
    ax1[0].set(ylabel='Hs (m)')
    ax1[1].plot(serie['year'], serie['Tp'], linewidth=1, color='black')
    ax1[1].set(xlabel='Year', ylabel='Tp (s)')

    max_annual_hs = pandas.DataFrame({'Year': [], 'max_hs': []})
    min_annual_hs = pandas.DataFrame({'Year': [], 'min_hs': []})
    max_annual_tp = pandas.DataFrame({'Year': [], 'max_tp': []})
    min_annual_tp = pandas.DataFrame({'Year': [], 'min_tp': []})

    id_max_hs = serie.groupby(lambda x: serie['year'][x].year)['Hs'].idxmax()
    max_annual_hs['max_hs'] = serie['Hs'].loc[id_max_hs]
    max_annual_hs['Year'] = serie['year'].loc[id_max_hs]

    id_min_hs = serie.groupby(lambda x: serie['year'][x].year)['Hs'].idxmin()
    min_annual_hs['min_hs'] = serie['Hs'].loc[id_min_hs]
    min_annual_hs['Year'] = serie['year'].loc[id_min_hs]

    id_max_tp = serie.groupby(lambda x: serie['year'][x].year)['Tp'].idxmax()
    max_annual_tp['max_tp'] = serie['Tp'].loc[id_max_tp]
    max_annual_tp['Year'] = serie['year'].loc[id_max_tp]

    id_min_tp = serie.groupby(lambda x: serie['year'][x].year)['Tp'].idxmin()
    min_annual_tp['min_tp'] = serie['Tp'].loc[id_min_tp]
    min_annual_tp['Year'] = serie['year'].loc[id_min_tp]

    #m1 = plt.get_current_fig_manager()
    #m1.window.showMaximized()

else:
    pass

##########################PLOT WIND ROSE OF DIRECTION##########################
# Figure 2

rep = 'nan'
while rep != 'y' or rep != 'n':
    rep = input('Figure 2: Plot wind rose of wave direction? y/n \t')
    if rep in ('y', 'n'):
        break

if rep == 'y':
    fig2 = plt.figure()
    ax2_0 = fig2.add_subplot(1, 2, 1, projection='windrose')
    ax2_0.bar(serie['Dir'], serie['Hs'], normed=True, bins=np.array([0, 1, 2, 3, 4, 5]),
              opening=0.8, linewidth=0)
    ax2_0.set_title('Windrose of Hs', weight='bold')
    ax2_0.set_legend(loc='lower left')
    ax2_1 = fig2.add_subplot(1, 2, 2, projection='windrose')
    ax2_1.bar(serie['Dir'], serie['Tp'], normed=True, bins=np.array([0, 6, 9, 12, 15, 18]),
              opening=0.8, linewidth=0)
    ax2_1.set_title('Windrose of Tp', weight='bold')
    ax2_1.set_legend(loc='lower right')

    #m2 = plt.get_current_fig_manager()
    #m2.window.showMaximized()
else:
    pass

################MULTIPLOT DENSITY EXCEEDANCE AND RETURN PERIOD#################
# Figure 3

rep = 'nan'
while rep != 'y' or rep != 'n':
    rep = input('Figure 3: Plot density,\
                probability of exceedance and return periods for Hs and Tp? y/n \t')
    if rep in ('y', 'n'):
        break

if rep == 'y':
    fig3, ax3 = plt.subplots(nrows=2, ncols=2)
    sns.kdeplot(data=serie, x='Hs', multiple='stack', ax=ax3[0, 0])
    sns.kdeplot(data=serie, x='Tp', multiple='stack', ax=ax3[0, 1])
    ax3[0, 0].set(xlabel='Hs (m)')
    ax3[0, 1].set(xlabel='Tp (s)')

    with sns.axes_style("dark"):
        sns.lineplot(data=serie_prob_Hs, x='Hs', y='exceedance_probability', ax=ax3[1, 0])
        ax3_1 = ax3[1, 0].twinx()
        ax3[1, 0].scatter(H12_y, 0.00136986, s=100, c='magenta',
                          label='Height exceeded 12 hours per year')
        sns.lineplot(data=serie_prob_Hs, x='Hs', y='return_period', ax=ax3_1, color='r')
        ax3_1.set(yscale='log')
        ax3[1, 0].set(xlabel='Hs (m)', ylabel='Probability of exceedance')
        ax3_1.set(ylabel='')
        sns.lineplot(data=serie_prob_Tp, x='Tp', y='exceedance_probability', ax=ax3[1, 1])
        ax3_2 = ax3[1, 1].twinx()
        ax3[1, 1].scatter(T12_y, 0.00136986, s=100, c='magenta',
                          label='Period exceeded 12 hours per year')
        sns.lineplot(data=serie_prob_Tp, x='Tp', y='return_period', ax=ax3_2, color='r')
        ax3_2.set(yscale='log')
        ax3[1, 1].set(xlabel='Tp (s)', ylabel='')
        ax3_2.set(ylabel='Return period (days)')

    #m3 = plt.get_current_fig_manager()
    #m3.window.showMaximized()

else:
    pass

###############################JOINT PROBABILITY###############################
# Figure 4

rep = 'nan'
while rep != 'y' or rep != 'n':
    rep = input('Figure 4: Plot joint repartition of (Hs, Tp) and levels of density? y/n \t')
    if rep in ('y', 'n'):
        break

fig4, ax4 = plt.subplots()
sns.scatterplot(data=serie, x='Tp', y='Hs', s=1, ax=ax4)
kde = sns.kdeplot(data=serie, x='Tp', y='Hs', color='r', ax=ax4,
                  cbar=True, cmap='tab10',
                  fill=False, levels=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.98])
ax4.set(xlabel='Tp (s)', ylabel='Hs (m)')

contours = []
for i in kde.get_children():
    if i.__class__.__name__ == 'QuadContourSet':
        contours.append(i.get_paths())

points_inside = []
#for k in range(len(contours[9])):
inside = []
for i in range(len(serie['Hs'])):
    inside.append(contours[0][9].contains_point((serie['Tp'][i], serie['Hs'][i])))
for i, x in enumerate(inside):
    if x:
        points_inside.append(i)

sns.scatterplot(data=serie.iloc[points_inside], x='Tp', y='Hs', color='red', s=1, ax=ax4)
#m4 = plt.get_current_fig_manager()
#m4.window.showMaximized()

if rep == 'y':
    print('Mean values of most likely couple (Hs, Tp):',
          serie.iloc[points_inside]['Hs'].mean(), 'm,',
          serie.iloc[points_inside]['Tp'].mean(), 's')

else:
    plt.close(fig4)

###########################OMEGA CALCULATION AND STATS##########################
MSR = float(cl["tide"]["MSR"])
D_50 = float(cl["sand_size"]["D50"])
RHOW = 1024
RHOS = 2650
VISK = 1.3e-3 / RHOW  # m2/s
D = D_50 * (9.81 * (RHOS / RHOW - 1) / (VISK ** 2)) ** (1.0 / 3.0)
WSED = VISK * (np.sqrt(10.36 ** 2 + 1.049 * D ** 3) - 10.36) / D_50  # m/s
serie['Hb'] = 0.39*9.81**(1/5)*(serie['Tp'].values*serie['Hs'].values*serie['Hs'].values)**(2/5)
serie['omega_hs'] = serie['Hs'].values / (serie['Tp'].values * WSED)
serie['omega'] = serie['Hb'].values / (serie['Tp'].values * WSED)

# Classification of omega in reflective, intermediate and dissipative
conditions = [
    (serie['omega'] <= 2),
    (serie['omega'] > 2) & (serie['omega'] <= 5),
    (serie['omega'] > 5)
    ]
state = ['reflective', 'intermediate', 'dissipative']
serie['state'] = np.select(conditions, state)

# Calculate season
serie['month'] = serie['year'].dt.month_name()

# Classification of the months into seasons, winter or summer
serie['season'] = np.where(serie['month'].str.contains('October|November|December|January|February|March'),
                                                       'Winter', 'Summer')

# Calculate the cumulative time from the beginning
serie['total_duration'] = pandas.to_timedelta(serie['year']
                                              - serie['year'][0]).astype('timedelta64[s]')
# Recalculate year from date because the column disappears
serie['year'] = serie['date'].apply(_from_ordinal)
# Calculate delta t between each row
serie['delta_t'] = serie['total_duration'].diff()
delta_t = serie['delta_t'].mean()

# Calculate the cumulative time for a beach state until it changes
serie['state_duration'] = serie['total_duration']
serie = serie.assign(k=[1]*len(serie.index))
for i in range(1, len(serie.index)):
    if serie.loc[i, 'state'] == serie.loc[i-1, 'state']:
        serie.loc[i, 'k'] = serie.loc[i-1, 'k'] + 1
    serie.loc[i-1,
              'state_duration'] = serie.loc[i-1,
                                            'total_duration'] - serie.loc[i-serie.loc[i-1, 'k'],
                                                                          'total_duration']

serie.iloc[-1, serie.columns.get_loc('state_duration')] = serie.iloc[-1, serie.columns.get_loc('total_duration')] -\
                                                          serie.iloc[-serie.iloc[-1, serie.columns.get_loc('k')],
                                                                     serie.columns.get_loc('total_duration')]
serie['state_duration'] = serie['state_duration'] / 3600 / 24  # Converts into days
serie['real_duration'] = serie['state_duration']

for i in range(len(serie.index)-1, 1, -1):  # question du zero non traitée
    if serie.loc[i, 'state_duration'] > serie.loc[i-1, 'state_duration']:
        serie.loc[i-1, 'real_duration'] = serie.loc[i, 'real_duration']

# Keep the state duration of the last element for each cluster
serie['real_duration'] = serie['real_duration'].mask(serie['real_duration'].shift(-1) == serie['real_duration'])

# Classification of the duration of each state into four classes
# You can modify the classes and labels
serie['duration_levels'] = pandas.cut(serie['real_duration'], [0, 1, 7, 15, 31],
                                      labels=['less than 1 day', 'between 1 and 7 days',
                                              'between 7 and 15 days', 'between 15 and 31 days'])

# Plot different Omega repartitions on Figure 4 of Masselink and Short (1992)
serie['RTR'] = MSR / serie['Hb']

rep = 'nan'
while rep != 'y' or rep != 'n':
    rep = input('Figure 5: Plot Omegas and highlight the most likely points? y/n \t')
    if rep in ('y', 'n'):
        break

if rep == 'y':
    # Figure 5: Plot omegas for the most likely couple (Hs, Tp)
    fig5, ax5 = plt.subplots()
    sns.scatterplot(data=serie, x='omega', y='RTR', s=1, ax=ax5)
    if 'points_inside' in locals() or 'points_inside' in globals():
        sns.scatterplot(data=serie.loc[points_inside], x='omega', y='RTR', s=1, color='red', ax=ax5)
    else:
        print('WARNING: Figure 6 was not plotted so the most likely points are unknown')
    masselink_plot(ax5)
else:
    pass

rep = 'nan'
while rep != 'y' or rep != 'n':
    rep = input('Figure 6: Plot comparison between Omega(Hs) and Omega(Hb)? y/n \t')
    if rep in ('y', 'n'):
        break

if rep == 'y':
    # Figure 6: Use omega(Hs) instead of omega(Hb)
    fig6, ax6 = plt.subplots()
    ax6.scatter(serie['omega_hs'], serie['RTR'], s=1, color='orange', marker='.')
    ax6.scatter(serie['omega'], serie['RTR'], s=1, marker='.')
    if 'points_inside' in locals() or 'points_inside' in globals():
        sns.scatterplot(data=serie.loc[points_inside], x='omega', y='RTR', s=1, color='red', ax=ax6)
    else:
        print('WARNING: Figure 6 was not plotted so the most likely points are unknown')
    masselink_plot(ax6)
else:
    pass

rep = 'nan'
while rep != 'y' or rep != 'n':
    rep = input('Figure 7: Plot Omegas repartition by season? y/n \t')
    if rep in ('y', 'n'):
        break

if rep == 'y':
    # Figure 7: Plot omegas by season (winter/summer)
    fig7, ax7 = plt.subplots()
    sns.scatterplot(data=serie, x='omega', y='RTR', s=10, hue='season', ax=ax7)
    if 'points_inside' in locals() or 'points_inside' in globals():
        sns.scatterplot(data=serie.loc[points_inside], x='omega', y='RTR', s=10, color='red', ax=ax7)
    else:
        print('WARNING: Figure 6 was not plotted so the most likely points are unknown')
    masselink_plot(ax7)
    ax7.legend(loc='lower right')
else:
    pass

rep = 'nan'
while rep != 'y' or rep != 'n':
    rep = input('Figure 8: Plot Omegas repartition by consecutive time? y/n \t')
    if rep in ('y', 'n'):
        break

if rep == 'y':
    # Figure 8: Plot omegas by consecutive time
    fig8, ax8 = plt.subplots()
    sns.scatterplot(data=serie, x='omega', y='RTR', hue='duration_levels',
                    palette='Set2', style='duration_levels', s=5, ax=ax8, size='duration_levels',
                    sizes={'less than 1 day': 5,
                           'between 1 and 7 days': 10,
                           'between 7 and 15 days': 30,
                           'between 15 and 31 days': 50})
    if 'points_inside' in locals() or 'points_inside' in globals():
        sns.scatterplot(data=serie.loc[points_inside], x='omega', y='RTR', s=1, color='red', ax=ax8)
    else:
        print('WARNING: Figure 6 was not plotted so the most likely points are unknown')
    masselink_plot(ax8)
else:
    pass

serie.to_csv('./stats_wavesduck_daily_8m_1990-2022.csv', sep=',', index=False)

rep = 'nan'
while rep != 'y' or rep != 'n':
    rep = input('Figure 9: Plot Omegas timeseries and monthly means of Omegas? y/n \t')
    if rep in ('y', 'n'):
        break

if rep == 'y':
    # Plot Omegas, Omegas(month)
    fig9, ax9 = plt.subplots(2)
    #ax9[0].plot(serie['year'], serie['omega'], linewidth=0.5)
    sns.scatterplot(data=serie, x='year', y='omega', hue='state', ax=ax9[0], s=10,
                    hue_order=['reflective', 'intermediate', 'dissipative'])
    legend_labels, _ = ax9[0].get_legend_handles_labels()
    ax9[0].legend(legend_labels, ['reflective',
                                  'intermediate',
                                  'dissipative'], bbox_to_anchor=(1, 1))
    sns.barplot(data=serie, x='month', y='omega', ax=ax9[1]) #hue='duration_levels',
    #            hue_order=['less than 1 day', 'between 1 and 7 days',
    #                       'between 7 and 15 days', 'between 15 and 31 days'],
    #            ci=None, linewidth=0)
    sns.barplot(data=serie, x='month', y='omega', ax=ax9[1],
                linewidth=1, linestyle='--', edgecolor='black', facecolor=(0,0,0,0), ci=None)
    legend_labels, _ = ax9[1].get_legend_handles_labels()
    #ax9[1].legend(legend_labels, ['less than 1 day',
    #                              'between 1 and 7 days',
    #                              'between 7 and 15 days',
    #                              'between 15 and 31 days'], bbox_to_anchor=(1, 1))
    ax9[0].set(xlabel='Year', ylabel='Dean number')
    ax9[1].set(xlabel='Month', ylabel='Dean number')
else:
    pass

plt.show()
print('End of analysis')
print('--- {} seconds ---'.format(time.time() - start_time))
