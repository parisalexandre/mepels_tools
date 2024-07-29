#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis of the morphodynamic evolution of a beach
over several months

@author: Alexandre Paris
"""

import glob
import sys
import time
import math
import datetime
import datefinder
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import yaml
from scipy.signal import find_peaks
from scipy.ndimage.filters import gaussian_filter
from scipy.optimize import curve_fit
import seaborn as sns
import pandas

sys.path.append('./libs/')
from indata_function import indata
from rotation_function import rotation
from interpolation_function import interpolation
from replacementValues_function import replacementValues

start_time = time.time()
sns.set_theme()
sns.set_context('talk')

###############################################################################


def _from_ordinal(serial_dn):
    """
    This function converts a serial date number into a datetime object

    The reference for datetime.fromordinal is 0001/01/01

    EXample:
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


def netcdf_open(elevation_netcdf_file):
    """
    Open and read the DEMs in NetCDF format
    """
    nc = xr.open_dataset(elevation_netcdf_file, drop_variables='project')
    xyz = nc.to_dataframe()
    nc.close()
    xyz = xyz.reset_index()
    return xyz


def func(xdata, coefa, coefb, coefc):
    """
    Function used by 'curve_fit' to fit the average profile
    ax^b + c

    Can be modified in Dean profile
    """
    return coefa * np.sign(xdata)*(np.abs(xdata))**coefb + coefc


def swap(bars, i, u, v):
    """
    Swap columns that contain bar u and bar v
    """
    (bars['crest_{}'.format(u)][i],
     bars['crest_{}'.format(v)][i]) = (bars['crest_{}'.format(v)][i],
                                       bars['crest_{}'.format(u)][i])
    (bars['trough_{}'.format(u)][i],
     bars['trough_{}'.format(v)][i]) = (bars['trough_{}'.format(v)][i],
                                        bars['trough_{}'.format(u)][i])
    (bars['z_crest_{}'.format(u)][i],
     bars['z_crest_{}'.format(v)][i]) = (bars['z_crest_{}'.format(v)][i],
                                         bars['z_crest_{}'.format(u)][i])
    (bars['z_trough_{}'.format(u)][i],
     bars['z_trough_{}'.format(v)][i]) = (bars['z_trough_{}'.format(v)][i],
                                          bars['z_trough_{}'.format(u)][i])
    (bars['h_bar_{}'.format(u)][i],
     bars['h_bar_{}'.format(v)][i]) = (bars['h_bar_{}'.format(v)][i],
                                       bars['h_bar_{}'.format(u)][i])
    (bars['w_bar_{}'.format(u)][i],
     bars['w_bar_{}'.format(v)][i]) = (bars['w_bar_{}'.format(v)][i],
                                       bars['w_bar_{}'.format(u)][i])
    return bars


def find_roots(x, y):
    s = np.abs(np.diff(np.sign(y))).astype(bool)
    return x[:-1][s] + np.diff(x)[s]/(np.abs(y[1:][s]/y[:-1][s])+1)


####################################################################################################
def analyse_morpho(X, Y, h, formula, date_bathy, Xthorn, Ebb, Hrms):
    """
    0. Prepare a plot for final figure (fig)
    0bis. Plot the results of this function in separated figures
        fig_profile and fig_bars
    1. Depth of closure
    2. Longshore-averaged cross-shore profile
    3. Slope on fitted profile at the depth of closure
    4. Depth of breaking (Thorton)
    5. Bar(s) detection on the averaged profile
    6. Bar(s) detection in 2D for each cross-shore profile
    7. Prepare and export bar data for future analysis
    """
    ##########
    # Plot initialization
    fig_profile, ax_profile = plt.subplots()
    fig = plt.figure(constrained_layout=True)
    #fig_profile.set_size_inches(12.8, 9.6, forward=True)
    fig.set_size_inches(25.6, 13.64, forward=True)
    gs = GridSpec(3, 2, figure=fig)
    ax4 = fig.add_subplot(gs[2, 1])

    ####################
    # Calculate the depth of closure and plot it
    if formula == 0:
        depth_closure = -(2.28 * H12_Y - 68.5 * ((H12_Y**2)
                                                 /
                                                 (G * (T12_Y**2))
                                                 )
                          )
        #ax4.axhline(depth_closure, label='DoC: Hallermeier')
        #ax_profile.axhline(depth_closure, label='DoC: Hallermeier')
    elif formula == 1:
        depth_closure = -(1.75 * H12_Y - 57.9 * ((H12_Y**2) / (G * (T12_Y**2))))
        ax4.axhline(depth_closure, label='DoC: Birkemeier')
        ax_profile.axhline(depth_closure, label='DoC: Birkemeier')
    elif formula == 2:
        depth_closure = -6.75 * H_mean
        ax4.axhline(depth_closure, label='DoC: Houston')
        ax_profile.axhline(depth_closure, label='DoC: Houston')
    elif formula == 3:
        depth_closure = -8.9 * H_mean
        ax4.axhline(depth_closure, label='DoC: Hallermeier approximated')
        ax_profile.axhline(depth_closure, label='DoC: Hallermeier approximated')
    else:
        print('ERROR: the depth of closure formula is not recognized')
        print('Please modify in the configuration file')
        sys.exit()

    ##########
    # Longshore-averaged cross-shore profile
    xpos = X[0, :]
    h = -h
    k = 0
    meanProfile = np.nanmean(h, axis=0)
    while k < np.shape(X)[0]:
        profile = h[k, :]
        if k == 1:
            ax4.plot(xpos, profile, color='0.8', linestyle='--',
                     label='cross-shore profile')
            ax_profile.plot(xpos, profile, color='0.8', linestyle='--',
                            label='cross-shore profile')
        else:
            ax4.plot(xpos, profile, color='0.8', linestyle='--')
            ax_profile.plot(xpos, profile, color='0.8', linestyle='--')
        k += 1

    if meanProfile[0] > depth_closure:
        depth_closure = meanProfile[0]

    ####################
    # Plot the cross-shore profiles, mean profile and different levels (HTL, LTL, water level)
    ax4.axhline(z_max, color='red', linestyle='-.', label='HTL')
    ax_profile.axhline(z_max, color='red', linestyle='-.', label='HTL')
    ax4.axhline(LTL_val, color='blue', linestyle='-.', label='LTL')
    ax_profile.axhline(LTL_val, color='blue', linestyle='-.', label='LTL')
    ax4.axhline(serie_eq.loc[serie_eq['year'] == date_bathy, 'Level'].values,
                color='orange', linestyle='-', label='water level')
    ax_profile.axhline(serie_eq.loc[serie_eq['year'] == date_bathy, 'Level'].values,
                       color='orange', linestyle='-', label='water level')
    ax4.plot(X[0, :], meanProfile, label='mean cross-shore profile', color='black')
    ax_profile.plot(X[0, :], meanProfile, label='mean cross-shore profile', color='black')

    ##########
    # Definition of the studied zone, between depth of closure and z_max
    k = 0
    while meanProfile[k] <= depth_closure:
        k += 1
    k_x_min = k
    if np.max(meanProfile) < z_max:
        k_x_max = np.argmax(meanProfile)
    else:
        k = 0
        while meanProfile[k] <= z_max:
            k += 1
        k_x_max = k

    profile_short = meanProfile[k_x_min:k_x_max]
    xpos_short = xpos[k_x_min:k_x_max]

    ####################
    # Fitted profile to calculate the slope at depth of closure
    popt, pcov = curve_fit(func, xpos_short, profile_short, check_finite=False)
    Beta = np.abs(popt[0]*popt[1]
                  *
                  (np.sign(xpos_short[0])*(np.abs(xpos_short[0]))**(popt[1] - 1)))

    ####################
    # Depth of breaking Thornton and Guza
    if math.isnan(h_breaking) is False:
        #ax4.axhline(-h_breaking, color='magenta', linestyle='-', label='breaking Thornton')
        ax4.axvline(x_breaking, color='magenta', linestyle='-')
        #ax_profile.axhline(-h_breaking, color='magenta', linestyle='-', label='breaking Thornton')
        ax_profile.axvline(x_breaking, color='magenta', linestyle='-', label='breaking Thornton')
        ax_bis_profile = ax_profile.twinx()
        ax_profile.plot(Xthorn, Ebb, label=r'$\epsilon_b$')
        #ax_bis_profile.plot(Xthorn, Hrms)
        ax_bis_profile.set_ylabel(r'Wave energy bore dissipation $\epsilon_b$')
        ax_bis_profile.set_ylim([-7.5, 7.5])
        ax_bis_profile.grid(False)

    ##########
    # Shoreline detection at z_max
    if date_bathy in serie_eq['year'].astype(str).unique() is False:
        shoreline = float('NaN')
    else:
        fr = find_roots(xpos, meanProfile - z_max)
        shoreline = fr[0]
    
    serie_eq.loc[serie_eq['year'] == date_bathy, 'shoreline'] = shoreline

    ####################
    # Bar(s) detection and plot on averaged profile
    peak, _ = find_peaks(profile_short, prominence=0.04)
    trough, _ = find_peaks(-profile_short, prominence=0.04)
    first_crests = ax4.plot(xpos_short[peak], profile_short[peak],
                            'o', color='r', label='bar crest')
    ax_profile.plot(xpos_short[peak], profile_short[peak],
                    'o', color='r', label='bar crest')
    first_troughs = ax4.plot(xpos_short[trough], profile_short[trough],
                             'o', color='b', label='bar trough')
    ax_profile.plot(xpos_short[trough], profile_short[trough],
                    'o', color='b', label='bar trough')
    ax4.set_xlabel('Cross-shore (m)')
    ax4.set_ylabel('Elevation (m)')
    ax4.set_xlim([-500, 100])  # WARNING
    ax4.set_ylim([-7.5, 7.5])  # WARNING
    ax4.grid(True)
    ax4.legend(loc='upper left')
    ax_profile.set_xlabel('Cross-shore (m)')
    ax_profile.set_ylabel('Elevation (m)')
    ax_profile.set_xlim([-500, 100])  # WARNING
    ax_profile.set_ylim([-7.5, 7.5])  # WARNING
    ax_profile.grid(True)
    ax_profile.legend(loc='upper left')
    ax_profile.set_title(date_bathy)
    #fig_profile.suptitle(date_bathy)
    #fig_profile.savefig('./figures_2000-2022_0-180/profile_bars_'+date_bathy+'.png', dpi=100)
    #plt.show()
    plt.close(fig_profile)

    ##########
    # Bar(s) in 2D for each cross-shore profile and plot
    coords_crests = []
    coords_troughs = []
    long_shore = Y[:, 0]
    for i in range(np.shape(X)[0]):  # One cross-shore profile by longshore coordinate
        profile = h[i, :]
        profile_short = profile[k_x_min:k_x_max]
        cross_short = xpos[k_x_min:k_x_max]
        # Profiles are filtered to avoid too much detection
        profile_short_filter = gaussian_filter(profile_short,
                                               sigma,
                                               mode='nearest')
        crests_bathy, _ = find_peaks(profile_short_filter,
                                     prominence=0.04)  # WARNING
        for _, n in enumerate(crests_bathy):
            coords_crests.append([long_shore[i],
                                  cross_short[n],
                                  profile_short_filter[n]])
        troughs_bathy, _ = find_peaks(-profile_short_filter,
                                      prominence=0.04)
        for _, n in enumerate(troughs_bathy):
            coords_troughs.append([long_shore[i],
                                   cross_short[n],
                                   profile_short_filter[n]])

    coords_crests = np.asarray(coords_crests)
    # Coords_crests is organized as follows:
    # first column: longshore coordinate of the profile
    # second column: cross-shore coordinate of the bar crest if detected
    # third column: elevation of the bar crest
    if len(coords_crests) > 0:
        # First sort on the elevation then on the longshore coordinate
        coords_crests = coords_crests[coords_crests[:, 2].argsort()]
        coords_crests = coords_crests[coords_crests[:, 0].argsort(kind='mergesort')]
    # Same thing for the troughs
    coords_troughs = np.asarray(coords_troughs)
    if len(coords_troughs) > 0:
        coords_troughs = coords_troughs[(coords_troughs[:, 2]).argsort()]
        coords_troughs = coords_troughs[coords_troughs[:, 0].argsort(kind='mergesort')]
    zminPlot = np.nanmin(h)
    zmaxPlot = np.nanmax(h)
    zlevels = np.linspace(-7.5, 7.5, 100)  # WARNING
    zclevels = np.arange(round(zminPlot), 1, 1)

    fig_bars, ax_bars = plt.subplots()
    ax = fig.add_subplot(gs[1, 0])
    fig_bars.set_size_inches(12.8, 9.6, forward=True)
    surf = ax.contourf(Y, X, h, cmap=cmapBathy, levels=zlevels)
    surf_bars = ax_bars.contourf(Y, X, h, cmap=cmapBathy, levels=zlevels)
    cbar = plt.colorbar(surf, spacing='uniform', orientation='vertical',
                        pad=0, ax=ax_bars,
                        ticks=[-7.5, -5, -2.5, 0, 2.5, 5, 7.5])
    cbar.ax.tick_params(labelsize=10)
    cs = ax.contour(Y, X, h, levels=zclevels, colors='black', linewidths=1)
    cs_bars = ax_bars.contour(Y, X, h, levels=zclevels, colors='black', linewidths=1)
    cs.clabel(levels=zclevels, colors='black', fontsize='xx-small')
    cs_bars.clabel(levels=zclevels, colors='black', fontsize='xx-small')
    
    # Thornton
    if math.isnan(h_breaking) is False:
        ax.contour(Y, X, h,
                   levels=[-h_breaking],
                   colors='magenta', linewidths=2)
        ax_bars.contour(Y, X, h,
                        levels=[-h_breaking],
                        colors='magenta', linewidths=2)

    fig_bars.suptitle(date_bathy)
    ax.set_ylabel('Cross-shore (m)')
    ax.set_xlabel('Long-shore (m)')
    ax.set_ylim([np.max(X), np.min(X)])
    ax.set_xlim([np.min(Y), np.max(Y)])
    ax_bars.set_ylabel('Cross-shore (m)')
    ax_bars.set_xlabel('Long-shore (m)')
    ax_bars.set_ylim([np.max(X), np.min(X)])
    ax_bars.set_xlim([np.min(Y), np.max(Y)])

    # Plot crests and troughs on 2D bathymetry if detected
    if len(coords_crests) > 0:
        for i in range(len(coords_crests)):
            ax.plot(coords_crests[i, 0], coords_crests[i, 1],
                    'x', color='red')
            ax_bars.plot(coords_crests[i, 0], coords_crests[i, 1],
                         'x', color='red')
        ax.plot(coords_crests[-1, 0], coords_crests[-1, 1], 'x', color='red',
                label='Bar crests')
        ax_bars.plot(coords_crests[-1, 0], coords_crests[-1, 1], 'x', color='red',
                     label='Bar crests')

    if len(coords_troughs) > 0:
        for i in range(len(coords_troughs)):
            ax.plot(coords_troughs[i, 0], coords_troughs[i, 1],
                    'x', color='cyan', markersize=5)
            ax_bars.plot(coords_troughs[i, 0], coords_troughs[i, 1],
                         'x', color='cyan', markersize=5)
        ax.plot(coords_troughs[-1, 0], coords_troughs[-1, 1], 'x', color='cyan',
                label='Bar troughs', markersize=5)
        ax_bars.plot(coords_troughs[-1, 0], coords_troughs[-1, 1], 'x', color='cyan',
                     label='Bar troughs', markersize=5)

    fig_bars.legend()
    #fig_bars.savefig('./figures_2000-2022_0-180/bathy_bars_'+date_bathy+'.png', dpi=100)
    #plt.show()
    plt.close(fig_bars)

    ####################
    # Clean the bars, calculate the dimensions and the number of bars
    # See bars_dimensions comments for more information
    if len(coords_troughs) > 0 and len(coords_crests) > 0:
        bars, nb_bar, fig, gs = bars_dimensions(coords_crests, coords_troughs, long_shore, fig, gs)
    else:
        bars = []
        nb_bar = 0

    return bars, nb_bar, fig, gs


###############################################################################
def bars_dimensions(coords_crests, coords_troughs, long_shore, fig, gs):
    """
    1. Create a pandas dataFrame with the results of analyse_morpho
    2. Sort and clean bars
    """
    ##########
    # Put the results of analyse_morpho in a pandas dataFrame
    serie_crests = pandas.DataFrame({'longshore': coords_crests[:, 0],
                                     'crest': coords_crests[:, 1],
                                     'z_crest': coords_crests[:, 2]})
    serie_troughs = pandas.DataFrame({'longshore': coords_troughs[:, 0],
                                      'trough': coords_troughs[:, 1],
                                      'z_trough': coords_troughs[:, 2]})
    serie_crests = serie_crests.set_index(['longshore',
                                          (serie_crests.groupby('longshore')
                                           .cumcount().add(1).astype(str))]).unstack()
    serie_troughs = serie_troughs.set_index(['longshore',
                                            (serie_troughs.groupby('longshore')
                                             .cumcount().add(1).astype(str))]).unstack()
    serie_crests.columns = serie_crests.columns.map('_'.join)
    serie_troughs.columns = serie_troughs.columns.map('_'.join)
    serie_crests = serie_crests.reset_index()
    serie_troughs = serie_troughs.reset_index()
    longshore = pandas.DataFrame({'longshore': long_shore})
    crest_trough = serie_crests.merge(serie_troughs, on='longshore', how='left')
    bars = longshore.merge(crest_trough, on='longshore', how='left')

    ####################
    # Initialize the number of bars
    if 'crest_1' in bars:
        nb_bar = 1
    else:
        nb_bar = 0

    # Calculate the height and the width for each bar and update the number of bars
    bars, nb_bar = bars_heights_widths(bars, nb_bar)

    ################################################################################################
    # From here, if there are two bars, bar 1 is supposed to be the farthest one
    # but many problems remain and are being addressed now
    ################################################################################################

    ##########
    # If only one bar, eliminate points too far (tolerance of 120 m) from the mean position
    # and between 0 and 50 m from the shore
    if nb_bar == 1:
        for i in range(len(bars)):
            if (math.isnan(bars['crest_1'][i]) is False and
                    math.isclose(bars['crest_1'][i],
                                 bars['crest_1'].mean(),
                                 abs_tol=120) is False and  # WARNING
                    bars['crest_1'][i] > -50):
                bars['crest_1'][i] = float('NaN')

    ####################
    # If two bars
    if nb_bar == 2:
        k = 0
        for i in range(len(bars)):  # Count how many times bar 1 is farther than bar 2
            if (math.isnan(bars['crest_1'][i]) is False and
                    math.isnan(bars['crest_2'][i]) is False and
                    bars['crest_1'][i] < bars['crest_2'][i]):
                k = k + 1
        if k < 29:  # WARNING
            # If bar 1 is not often enough farther than bar 2, they are swapped
            for i in range(len(bars)):
                if (math.isnan(bars['crest_1'][i]) is False and
                        math.isnan(bars['crest_2'][i]) is False and
                        bars['crest_1'][i] < bars['crest_2'][i] and
                        bars['crest_1'].isnull().sum() > (len(bars) - 10)):
                    # Points where bar 1 is farther than bar 2
                    # while there is too much NaN value in bar 1
                    # i.e bar 1 is not really a bar
                    bars['crest_1'][i] = float('NaN')
                    bars = swap(bars, i, 1, 2)
                elif (math.isnan(bars['crest_1'][i]) is False and
                        math.isnan(bars['crest_2'][i]) is False and
                        bars['crest_1'][i] < bars['crest_2'][i] and
                        bars['crest_2'].isnull().sum() > (len(bars) - 10)):
                    # Points where bar 1 is farther than bar 2
                    # while there is too much NaN value in bar 2
                    bars = swap(bars, i, 1, 2)

        else:
            k = 0
            for i in range(len(bars)):  # Same as before but now when bar 2 is farther than bar 1
                if (math.isnan(bars['crest_1'][i]) is False and
                        math.isnan(bars['crest_2'][i]) is False and
                        bars['crest_1'][i] > bars['crest_2'][i]):
                    k = k + 1
                if k < 20:  # WARNING
                    for i in range(len(bars)):
                        if (math.isnan(bars['crest_1'][i]) is False and
                                math.isnan(bars['crest_2'][i]) is False and
                                bars['crest_1'][i] > bars['crest_2'][i] and
                                bars['crest_2'].isnull().sum() > (len(bars) - 10)):
                            bars = swap(bars, i, 1, 2)

        # Bar 2 is deleted if is NaN or contains to much NaN
        if (bars['crest_2'].isnull().sum() == len(bars) or
                bars['crest_2'].isnull().sum() > (len(bars) - 25)):
            bars = bars.drop(columns=['crest_2', 'z_crest_2', 'trough_2',
                                      'z_trough_2', 'h_bar_2', 'w_bar_2'])
            nb_bar = 1

    ##########
    # If one bar (detected before or remaining after first previous processing)
    # delete points too far from the mean, with a tolerance of 100 m
    if nb_bar == 1:
        mean_crest_1 = bars['crest_1'].mean()
        for i in range(len(bars)):
            if (math.isnan(bars['crest_1'][i]) is False and
                    math.isclose(bars['crest_1'][i],
                                 mean_crest_1, abs_tol=100) is False):  # WARNING
                bars['crest_1'][i] = float('NaN')

    ####################
    # If two bars
    if nb_bar == 2:
        if 'trough_3' in bars:
            # Sometimes, 'trough_3' can remain whereas there are only two bars left
            bars = bars.drop(columns=['trough_3', 'z_trough_3'])
        if bars['crest_1'].mean() > bars['crest_2'].mean():
            # Swap bars if mean position of bar 1 is closer to the shore than bar 2
            for i in range(len(bars)):
                bars = swap(bars, i, 1, 2)
        for i in range(len(bars)):
            if bars['crest_1'][i] > bars['crest_2'][i]:
                # For each point, swap if bar 1 closer to the shore than bar 2
                bars = swap(bars, i, 1, 2)

        mean_crest_1 = bars['crest_1'].mean()
        mean_crest_2 = bars['crest_2'].mean()
        for i in range(len(bars)):
            # Check all points of bar 1 when points of bar 2 are NaN
            if math.isnan(bars['crest_2'][i]):
                if math.isclose(bars['crest_1'][i], mean_crest_1, abs_tol=50) is False:
                    # Swap if point of bar 1 is too far from mean position of bar 1
                    bars = swap(bars, i, 1, 2)
                elif math.isclose(bars['crest_1'][i], mean_crest_2, abs_tol=50) is True:
                    # Swap if point of bar 1 is close to mean position of bar 2
                    bars = swap(bars, i, 1, 2)
                elif math.isclose(bars['crest_1'][i], bars['crest_1'].min(), abs_tol=50) is False:
                    # Swap if point of bar 1 is too far from minimum of bar 1
                    bars = swap(bars, i, 1, 2)
                elif math.isclose(bars['crest_1'][i], bars['crest_2'].max(), abs_tol=50) is True:
                    # Swap if point of bar 1 is close to maximum of bar 2
                    bars = swap(bars, i, 1, 2)
            # Update mean values to take into account the swaps
            mean_crest_1 = bars['crest_1'].mean()
            mean_crest_2 = bars['crest_2'].mean()

        for i in range(len(bars)):
            # Same as previous for loop but now for points of bar 2 when points of bar 1 are NaN
            if math.isnan(bars['crest_1'][i]):
                if math.isclose(bars['crest_2'][i], mean_crest_2, abs_tol=50) is False:
                    bars = swap(bars, i, 1, 2)
                elif math.isclose(bars['crest_2'][i], mean_crest_1, abs_tol=50) is True:
                    bars = swap(bars, i, 1, 2)
                elif (math.isclose(bars['crest_2'][i],
                      bars['crest_2'].max(), abs_tol=50) is False and
                        math.isclose(bars['crest_2'][i], mean_crest_1, abs_tol=50) is True):
                    # Swap if point of bar 2 is too far from max of bar 2 and close to mean of bar 1
                    bars = swap(bars, i, 1, 2)
                elif math.isclose(bars['crest_2'][i], bars['crest_1'].min(), abs_tol=50) is True:
                    bars = swap(bars, i, 1, 2)
            mean_crest_1 = bars['crest_1'].mean()
            mean_crest_2 = bars['crest_2'].mean()

        # Two previous for loops are repeated one time each with some variations
        for i in range(len(bars)):
            if math.isnan(bars['crest_2'][i]):
                if math.isclose(bars['crest_1'][i], mean_crest_1, abs_tol=50) is False:
                    bars = swap(bars, i, 1, 2)
                elif math.isclose(bars['crest_1'][i], mean_crest_2, abs_tol=50) is True:
                    bars = swap(bars, i, 1, 2)
                elif math.isclose(bars['crest_1'][i], bars['crest_2'].max(), abs_tol=50) is True:
                    bars = swap(bars, i, 1, 2)
            mean_crest_1 = bars['crest_1'].mean()
            mean_crest_2 = bars['crest_2'].mean()

        for i in range(len(bars)):
            if math.isnan(bars['crest_1'][i]):
                if math.isclose(bars['crest_2'][i], mean_crest_2, abs_tol=50) is False:
                    bars = swap(bars, i, 1, 2)
                elif math.isclose(bars['crest_2'][i], mean_crest_1, abs_tol=50) is True:
                    bars = swap(bars, i, 1, 2)
                elif math.isclose(bars['crest_2'][i], bars['crest_2'].max(), abs_tol=50) is False:
                    bars = swap(bars, i, 1, 2)
            mean_crest_1 = bars['crest_1'].mean()
            mean_crest_2 = bars['crest_2'].mean()

        # Other comparisons
        for i in range(len(bars)):
            if (math.isnan(bars['crest_2'][i]) and
                    abs(bars['crest_1'][i] - mean_crest_1) >
                    abs(bars['crest_1'][i] - mean_crest_2)):
                # Swap if point of bar 1 is closer to mean of bar 2 than mean of bar 1
                bars = swap(bars, i, 1, 2)
            elif (math.isnan(bars['crest_2'][i]) and
                    (abs(bars['crest_1'][i] - bars['crest_1'].min()) >
                     abs(bars['crest_1'][i] - bars['crest_2'].min())) and
                    bars['crest_1'].min() > -300 and
                    math.isclose(bars['crest_1'][i], mean_crest_1, abs_tol=30) is False):  # WARNING
                # Swap if point of bar 1 is closer to min of bar 2 than min of bar 1
                # and min of bar 1 is located between 0 and 300 m
                # and point of bar 1 is too far from mean of bar 1
                bars = swap(bars, i, 1, 2)
            elif (math.isnan(bars['crest_1'][i]) and
                    abs(bars['crest_2'][i] - mean_crest_2) >
                    abs(bars['crest_2'][i] - mean_crest_1)):
                # Swap if point of bar 2 is closer to mean of bar 1 than mean of bar 2
                bars = swap(bars, i, 1, 2)
            elif (math.isnan(bars['crest_1'][i]) and
                    abs(bars['crest_2'][i] - bars['crest_2'].min()) >
                    abs(bars['crest_2'][i] - bars['crest_1'][i].min())):
                # Swap if point of bar 2 is closer to min of bar 1 than min of bar 2
                bars = swap(bars, i, 1, 2)
            mean_crest_1 = bars['crest_1'].mean()
            mean_crest_2 = bars['crest_2'].mean()

    ##########
    # The two bars are sorted again after all previous processings
    if nb_bar == 2:
        if bars['crest_1'].mean() > bars['crest_2'].mean():
            bars.columns = ['longshore', 'crest_2', 'crest_1', 'z_crest_2', 'z_crest_1',
                            'trough_2', 'trough_1', 'z_trough_2', 'z_trough_1',
                            'h_bar_2', 'h_bar_1', 'w_bar_2', 'w_bar_1']

    ####################
    # Bars too small are removed and number of bars is updated
    if nb_bar > 0 and nb_bar != 3:
        for i in range(1, nb_bar+1):
            if bars['crest_{}'.format(i)].count() < (len(bars)/2):  # WARNING
                bars = bars.drop(columns=['crest_{}'.format(i),
                                          'trough_{}'.format(i),
                                          'z_crest_{}'.format(i),
                                          'z_trough_{}'.format(i),
                                          'h_bar_{}'.format(i),
                                          'w_bar_{}'.format(i)])
                nb_bar = nb_bar - 1

        if ('crest_1' in bars) is False and nb_bar == 2 and ('crest_2' in bars):
            # If bars were swapped and the only one remaining is bar 2, it is renamed as bar 1
            bars.columns = ['longshore', 'crest_1', 'z_crest_1', 'trough_1',
                            'z_trough_1', 'h_bar_1', 'w_bar_1']

    ##########
    # Cases with three bars
    if nb_bar == 3:
        for i in range(len(bars)):
            if (math.isnan(bars['crest_1'][i]) is False and
                math.isnan(bars['crest_2'][i]) is False and
                    math.isnan(bars['crest_3'][i]) is False):
                # Bar 1 is supposed to be the farthest from the shore
                if bars['crest_1'][i] > bars['crest_3'][i]:
                    bars = swap(bars, i, 1, 3)
                    for j in range(len(bars)):
                        if bars['crest_1'][j] > bars['crest_2'][j]:
                            bars = swap(bars, j, 1, 2)
                            for k in range(len(bars)):
                                if bars['crest_2'][k] > bars['crest_3'][k]:
                                    bars = swap(bars, k, 2, 3)
            elif (math.isnan(bars['crest_1'][i]) is False and
                    math.isnan(bars['crest_2'][i]) is False and
                    math.isnan(bars['crest_3'][i])):
                if bars['crest_1'][i] > bars['crest_2'][i]:
                    bars = swap(bars, i, 1, 2)
        mean_crest_1 = bars['crest_1'].mean()
        mean_crest_2 = bars['crest_2'].mean()

        for i in range(len(bars)):
            if (math.isnan(bars['crest_1'][i]) is False and
                    math.isnan(bars['crest_2'][i]) and
                    math.isnan(bars['crest_3'][i]) and
                    math.isclose(bars['crest_1'][i],
                                 mean_crest_1,
                                 abs_tol=50) is False):
                # Swap with bar 2 if only point of bar 1 is non NaN and too far from mean of bar 1
                bars = swap(bars, i, 1, 2)
            mean_crest_1 = bars['crest_1'].mean()
            mean_crest_2 = bars['crest_2'].mean()
        for i in range(len(bars)):
            if (math.isnan(bars['crest_1'][i]) and
                    math.isnan(bars['crest_2'][i]) is False and
                    math.isnan(bars['crest_3'][i]) and
                    math.isclose(bars['crest_2'][i],
                                 bars['crest_2'].max(),
                                 abs_tol=100) is False):
                # Swap with bar 1 if only point of bar 2 is non NaN and is too far from max of bar 2
                bars = swap(bars, i, 1, 2)
            mean_crest_1 = bars['crest_1'].mean()
            mean_crest_2 = bars['crest_2'].mean()

    if nb_bar == 3:
        # If bar 3 is too small, it is removed and number of bars is updated
        if bars['crest_3'].isnull().sum() > (len(bars) - 10):
            bars = bars.drop(columns=['crest_3', 'z_crest_3', 'trough_3',
                                      'z_trough_3', 'h_bar_3', 'w_bar_3'])
            nb_bar = 2

    #####################################
    # KEEP ALL BARS BETWEEN 0 AND 180 M #
    #####################################
    #if nb_bar == 3 and 'crest_1' in bars and 'crest_2' in bars and 'crest_3' in bars:
    #    bars = bars.drop(columns=['crest_1', 'trough_1', 'z_crest_1',
    #                              'z_trough_1', 'h_bar_1', 'w_bar_1',
    #                              'crest_3', 'trough_3', 'z_crest_3',
    #                              'z_trough_3', 'h_bar_3', 'w_bar_3' ])
    #    bars.columns = ['longshore', 'crest_1', 'z_crest_1', 'trough_1',
    #                    'z_trough_1', 'h_bar_1', 'w_bar_1']
    #    nb_bar = 1
    #
    #if nb_bar == 2 and 'crest_1' in bars and 'crest_2' in bars:
    #    bars = bars.drop(columns=['crest_1', 'trough_1', 'z_crest_1',
    #                              'z_trough_1', 'h_bar_1', 'w_bar_1'])
    #    bars.columns = ['longshore', 'crest_1', 'z_crest_1', 'trough_1',
    #                    'z_trough_1', 'h_bar_1', 'w_bar_1']
    #    nb_bar = 1

    #if nb_bar == 1 and 'crest_1' in bars:
    #    for i in range(len(bars)):
    #        if bars['crest_1'][i] < -180:
    #            bars['crest_1'][i] = 'NaN'


    #if 'crest_1' in bars and bars['crest_1'].count() < (len(bars)/5):
    #    bars = bars.drop(columns=['crest_1', 'trough_1', 'z_crest_1',
    #                              'z_trough_1', 'h_bar_1', 'w_bar_1'])

    # Remove all bars between 0 and 180 m
    #if nb_bar == 3 and 'crest_1' in bars and 'crest_2' in bars and 'crest_3' in bars:
    #    bars = bars.drop(columns=['crest_2', 'trough_2', 'z_crest_2',
    #                              'z_trough_2', 'h_bar_2', 'w_bar_2',
    #                              'crest_3', 'trough_3', 'z_crest_3',
    #                              'z_trough_3', 'h_bar_3', 'w_bar_3' ])
    #    bars.columns = ['longshore', 'crest_1', 'z_crest_1', 'trough_1',
    #                    'z_trough_1', 'h_bar_1', 'w_bar_1']
    #    nb_bar = 1
    #
    #if nb_bar == 2 and 'crest_1' in bars and 'crest_2' in bars:
    #    bars = bars.drop(columns=['crest_2', 'trough_2', 'z_crest_2',
    #                              'z_trough_2', 'h_bar_2', 'w_bar_2'])
    #    bars.columns = ['longshore', 'crest_1', 'z_crest_1', 'trough_1',
    #                    'z_trough_1', 'h_bar_1', 'w_bar_1']
    #    nb_bar = 1

    #if nb_bar == 1 and 'crest_1' in bars:
    #    for i in range(len(bars)):
    #        if bars['crest_1'][i] > -180:
    #            bars['crest_1'][i] = 'NaN'

    #if 'crest_1' in bars and bars['crest_1'].count() < (len(bars)/5):
    #    bars = bars.drop(columns=['crest_1', 'trough_1', 'z_crest_1',
    #                              'z_trough_1', 'h_bar_1', 'w_bar_1'])
    
    #################################################################
    # KEEP BAR CORRESPONDING TO THORNTON AND GUZA DEPTH OF BREAKING #
    # AND SET crest_1 = x_breaking ##################################
    #################################################################
    if nb_bar > 1:
        if 'crest_3' in bars:
            if math.isclose(bars['crest_3'].mean(), x_breaking, abs_tol=20):
                bars = bars.drop(columns=['crest_1', 'trough_1', 'z_crest_1',
                                          'z_trough_1', 'h_bar_1', 'w_bar_1'
                                          'crest_2', 'trough_2', 'z_crest_2',
                                          'z_trough_2', 'h_bar_2', 'w_bar_2'])
                bars = bars.rename(columns={'crest_3': 'crest_1', 'trough_3': 'trough_1',
                                            'z_crest_3': 'z_crest_1', 'z_trough_3': 'z_trough_1',
                                            'h_bar_3': 'h_bar_1', 'w_bar_3': 'w_bar_1'})
                bars['crest_1'] = x_breaking
            elif math.isclose(bars['crest_2'].mean(), x_breaking, abs_tol=20):
                bars = bars.drop(columns=['crest_1', 'trough_1', 'z_crest_1',
                                          'z_trough_1', 'h_bar_1', 'w_bar_1'
                                          'crest_3', 'trough_3', 'z_crest_3',
                                          'z_trough_3', 'h_bar_3', 'w_bar_3'])
                bars = bars.rename(columns={'crest_2': 'crest_1', 'trough_2': 'trough_1',
                                            'z_crest_2': 'z_crest_1', 'z_trough_2': 'z_trough_1',
                                            'h_bar_2': 'h_bar_1', 'w_bar_2': 'w_bar_1'})
                bars['crest_1'] = x_breaking
            elif math.isclose(bars['crest_1'].mean(), x_breaking, abs_tol=20):
                bars = bars.drop(columns=['crest_2', 'trough_2', 'z_crest_2',
                                          'z_trough_2', 'h_bar_2', 'w_bar_2'
                                          'crest_3', 'trough_3', 'z_crest_3',
                                          'z_trough_3', 'h_bar_3', 'w_bar_3'])
                bars['crest_1'] = x_breaking
        elif 'crest_2' in bars:
            if math.isclose(bars['crest_2'].mean(), x_breaking, abs_tol=20):
                bars = bars.drop(columns=['crest_1', 'trough_1', 'z_crest_1',
                                          'z_trough_1', 'h_bar_1', 'w_bar_1'])
                bars = bars.rename(columns={'crest_2': 'crest_1', 'trough_2': 'trough_1',
                                            'z_crest_2': 'z_crest_1', 'z_trough_2': 'z_trough_1',
                                            'h_bar_2': 'h_bar_1', 'w_bar_2': 'w_bar_1'})
                bars['crest_1'] = x_breaking
            elif math.isclose(bars['crest_1'].mean(), x_breaking, abs_tol=20):
                bars = bars.drop(columns=['crest_2', 'trough_2', 'z_crest_2',
                                          'z_trough_2', 'h_bar_2', 'w_bar_2'])
                bars['crest_1'] = x_breaking
        elif math.isclose(bars['crest_1'].mean(), x_breaking, abs_tol=20):
            bars['crest_1'] = x_breaking
        else:
            bars = bars.drop(columns=['crest_1', 'trough_1', 'z_crest_1',
                                      'z_trough_1', 'h_bar_1', 'w_bar_1'])
    
    ############################################
    # Set bar = Thornton. If positive -> no bar
    #if x_breaking < 0:
    #    bars['crest_1'] = x_breaking
    #    nb_bar = 1
    #else:
    #    if 'crest_1' in bars:
    #        bars = bars.drop(columns=['crest_1', 'trough_1', 'z_crest_1',
    #                                  'z_trough_1', 'h_bar_1', 'w_bar_1'])
    #    nb_bar = 0
    
    ###########################
    # PLOT DIMENSIONS OF BARS #
    ###########################
    fig_dim, ax_dim = plt.subplots()
    fig_dim.set_size_inches(12.8, 9.6, forward=True)
    ax1 = fig.add_subplot(gs[1, 1])
    ax2 = ax1.twinx()
    ax2_dim = ax_dim.twinx()
    if (nb_bar > 0) and 'h_bar_1' in bars:
        ax1.plot(bars['longshore'], bars['h_bar_1'],
                 label='outer bar height',
                 linestyle='--')
        ax2.plot(bars['longshore'], bars['w_bar_1'],
                 label='outer bar width')
        ax_dim.plot(bars['longshore'], bars['h_bar_1'],
                    label='outer bar height',
                    linestyle='--')
        ax2_dim.plot(bars['longshore'], bars['w_bar_1'],
                     label='outer bar width')
        if (nb_bar == 2) and 'h_bar_2' in bars:
            ax1.plot(bars['longshore'], bars['h_bar_2'],
                     label='inner bar height',
                     linestyle='--')
            ax2.plot(bars['longshore'], bars['w_bar_2'],
                     label='inner bar width')
            ax_dim.plot(bars['longshore'], bars['h_bar_2'],
                        label='inner bar height',
                        linestyle='--')
            ax2_dim.plot(bars['longshore'], bars['w_bar_2'],
                         label='inner bar width')
        else:
            for i in range(2, nb_bar):
                if 'h_bar_{}'.format(i) in bars:
                    ax1.plot(bars['longshore'], bars['h_bar_{}'.format(i)],
                             label='inner bar {} height'.format(i-1))
                    ax2.plot(bars['longshore'], bars['w_bar_{}'.format(i)],
                             linestyle='--', color='r',
                             label='inner bar {} width'.format(i-1))
                    ax_dim.plot(bars['longshore'], bars['h_bar_{}'.format(i)],
                                label='inner bar {} height'.format(i-1))
                    ax2_dim.plot(bars['longshore'], bars['w_bar_{}'.format(i)],
                                 linestyle='--', color='r',
                                 label='inner bar {} width'.format(i-1))

    ax1.set_xlabel('Long-shore coordinate (m)')
    ax1.set_ylabel('Bar height (m)')
    ax2.set_ylabel('Bar width (m)')
    ax1.set_ylim([0, 2])
    ax2.set_ylim([0, 200])
    ax1.set_xlim([610, 1070])
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')

    ax_dim.set_xlabel('Long-shore coordinate (m)')
    ax_dim.set_ylabel('Bar height (m)')
    ax2_dim.set_ylabel('Bar width (m)')
    ax_dim.set_ylim([0, 2])
    ax2_dim.set_ylim([0, 200])
    ax_dim.set_xlim([610, 1070])
    ax_dim.legend(loc='upper left')
    ax2_dim.legend(loc='upper right')
    fig_dim.suptitle(date_bathy)
    #fig_dim.savefig('./figures_2000-2022_0-180/bars_dimensions_'+date_bathy+'.png')
    plt.close(fig_dim)

    ####################
    # Plot remaining bars with RMS on bar position
    fig_position, ax_position = plt.subplots()
    fig_position.set_size_inches(12.8, 9.6, forward=True)
    ax3 = fig.add_subplot(gs[2, 0])
    if (nb_bar > 0) and 'crest_1' in bars:
        rms = np.sqrt(np.nansum((bars['crest_1'] - bars['crest_1'].mean())**2) / len(bars))
        ax3.axhline(y=bars['crest_1'].mean())
        ax3.axhline(y=(bars['crest_1'].mean() + rms), linestyle='--')
        ax3.axhline(y=(bars['crest_1'].mean() - rms), linestyle='--')
        ax3.scatter(bars['longshore'], bars['crest_1'],
                    label='outer bar')
        ax_position.axhline(y=bars['crest_1'].mean())
        ax_position.axhline(y=(bars['crest_1'].mean() + rms), linestyle='--')
        ax_position.axhline(y=(bars['crest_1'].mean() - rms), linestyle='--')
        ax_position.scatter(bars['longshore'], bars['crest_1'],
                            label='outer bar')
        if (nb_bar == 2) and 'crest_2' in bars:
            rms = np.sqrt(np.nansum((bars['crest_2'] - bars['crest_2'].mean())**2) / len(bars))
            ax3.axhline(y=bars['crest_2'].mean(), color='r')
            ax3.axhline(y=(bars['crest_2'].mean() + rms), color='r', linestyle='--')
            ax3.axhline(y=(bars['crest_2'].mean() - rms), color='r', linestyle='--')
            ax3.scatter(bars['longshore'], bars['crest_2'], color='r',
                        label='inner bar')
            ax_position.axhline(y=bars['crest_2'].mean(), color='r')
            ax_position.axhline(y=(bars['crest_2'].mean() + rms), color='r', linestyle='--')
            ax_position.axhline(y=(bars['crest_2'].mean() - rms), color='r', linestyle='--')
            ax_position.scatter(bars['longshore'], bars['crest_2'], color='r',
                                label='inner bar')
        if (nb_bar == 3) and 'crest_3' in bars:
            rms = np.sqrt(np.nansum((bars['crest_2'] - bars['crest_2'].mean())**2) / len(bars))
            ax3.axhline(y=bars['crest_2'].mean(), color='r')
            ax3.axhline(y=(bars['crest_2'].mean() + rms), color='r', linestyle='--')
            ax3.axhline(y=(bars['crest_2'].mean() - rms), color='r', linestyle='--')
            ax3.scatter(bars['longshore'], bars['crest_2'], color='r',
                        label='inner bar')
            rms = np.sqrt(np.nansum((bars['crest_3'] - bars['crest_3'].mean())**2) / len(bars))
            ax3.axhline(y=bars['crest_3'].mean(), color='orange')
            ax3.axhline(y=(bars['crest_3'].mean() + rms), color='orange', linestyle='--')
            ax3.axhline(y=(bars['crest_3'].mean() - rms), color='orange', linestyle='--')
            ax3.scatter(bars['longshore'], bars['crest_3'], color='orange',
                        label='inner bar 2')

            ax_position.axhline(y=bars['crest_2'].mean(), color='r')
            ax_position.axhline(y=(bars['crest_2'].mean() + rms), color='r', linestyle='--')
            ax_position.axhline(y=(bars['crest_2'].mean() - rms), color='r', linestyle='--')
            ax_position.scatter(bars['longshore'], bars['crest_2'], color='r',
                                label='inner bar')
            ax_position.axhline(y=bars['crest_3'].mean(), color='orange')
            ax_position.axhline(y=(bars['crest_3'].mean() + rms), color='orange', linestyle='--')
            ax_position.axhline(y=(bars['crest_3'].mean() - rms), color='orange', linestyle='--')
            ax_position.scatter(bars['longshore'], bars['crest_3'], color='orange',
                                label='inner bar 2')
    ax3.legend(loc='upper right')
    ax3.set_ylabel('Cross-shore (m)')
    ax3.set_xlabel('Long-shore (m)')
    ax3.set_ylim([np.max(X), np.min(X)])
    ax3.set_xlim([610, 1070])  # WARNING
    ax_position.legend(loc='upper right')
    ax_position.set_ylabel('Cross-shore (m)')
    ax_position.set_xlabel('Long-shore (m)')
    ax_position.set_ylim([np.max(X), np.min(X)])
    ax_position.set_xlim([610, 1070])  # WARNING
    fig_position.suptitle(date_bathy)
    #fig_position.savefig('./figures_2000-2022_0-180/bars_positions_'+date_bathy+'.png')
    plt.close(fig_position)

    #######################
    # FRAGMENTATION INDEX #
    #######################
    if nb_bar > 0 and 'crest_1' in bars:
        # Identify groups of cumulative NaNs
        na_groups = bars['crest_1'].notna().cumsum()[bars['crest_1'].isna()]
        # Length of NaNs groups over total length
        lengths_consecutive_na = na_groups.groupby(na_groups).agg(len)
        longest_na_ga = lengths_consecutive_na.max()
        bars['frag_1'] = len(na_groups) / len(bars)
        if nb_bar == 2 and 'crest_2' in bars:
            na_groups = bars['crest_2'].notna().cumsum()[bars['crest_2'].isna()]
            lengths_consecutive_na = na_groups.groupby(na_groups).agg(len)
            longest_na_ga = lengths_consecutive_na.max()
            bars['frag_2'] = len(na_groups) / len(bars)
    else:
        bars['frag_1'] = float('NaN')
        bars['frag_2'] = float('NaN')

    return bars, nb_bar, fig, gs


###############################################################################
def bars_heights_widths(bars, nb_bar):
    while ('crest_{}'.format(nb_bar) in bars) and ('trough_{}'.format(nb_bar) in bars):
        bars['h_bar_{}'.format(nb_bar)] = np.abs(bars['z_trough_{}'.format(nb_bar)] -
                                                 bars['z_crest_{}'.format(nb_bar)])
        bars['w_bar_{}'.format(nb_bar)] = np.abs(bars['crest_{}'.format(nb_bar)] -
                                                 bars['trough_{}'.format(nb_bar)])
        nb_bar += 1

    nb_bar -= 1
    return bars, nb_bar


###############################################################################
def omega_eq(serie):
    deltat = round(serie['delta_t'].mean())
    num_per_day = round(86400 / deltat)
    #phi = int(input('Enter the value of phi calculated by ShoreFor: '))
    phi = 250
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


###############################################################################
def airy(hi, k, omega):
    while omega**2 - G * k * np.tanh(k * hi) > 0.0001:
        k = k + 0.0001
    cel = math.sqrt(G * np.tanh(k * hi) / k)
    return cel, k


def E(Hrms_i):
    return RHOW*G*(Hrms_i**2) / 8


def Cg(hi, k, omega, theta):
    cel, k = airy(hi, k, omega)
    Cg = cel*np.cos(theta)*(1 + 2 * k/(np.sinh(2 * k * hi)))/2
    return Cg


def Eb(hi, Hrms, gamma, B, f):
    par = 1 - 1 / (1 + (Hrms / (gamma * hi))**2)**(5/2)
    return par * (3*math.sqrt(np.pi)/16)*RHOW*G*(B**3)*f*(Hrms**5)/((gamma**2) * (hi**3))


###############################################################################
def thornton(serie, X, h):
    '''
    Calculate the depth of breaking following Thornton and Guza 1983
    '''
    h = np.nanmean(h, axis=0)
    h = h[h >= 0]
    X = X[0, 0:len(h)]
    B = 0.8
    theta = math.radians(serie.loc[serie['year'] == date_bathy, 'Dir'].values.mean())
    f = 1/serie.loc[serie['year'] == date_bathy, 'Tp'].values.mean()
    omega = 2 * np.pi * f
    Hs = serie.loc[serie['year'] == date_bathy, 'Hs'].values.mean()
    Hrms_0 = Hs / 1.42

    gamma = 0.42
    
    k = omega**2 / G
    Hrms = np.zeros(len(X))
    Hrms[0] = Hrms_0
    i = 1
    while Hrms[i-1] < 1 and i < len(X):
        E_i_1 = E(Hrms[i-1])
        Cg_i_1 = Cg(h[i-1], k, omega, theta)
        Eb_i_1 = Eb(h[i-1], Hrms[i-1], gamma, B, f)
        try:
            Hrms[i] = math.sqrt(8 * ((E_i_1*Cg_i_1) +
                                     (np.abs(X[i] - X[i-1]))*Eb_i_1)/(RHOW*G*Cg(h[i], k, omega, theta)))
        except ValueError:
            Hrms[i] = 0

        i = i + 1

    #fig = plt.figure()
    #plt.plot(X, Hrms)

    Eb_bathy = np.zeros(len(Hrms))
    for i in range(len(h)):
        Eb_bathy[i] = Eb(h[i], Hrms[i], gamma, B, f)
    
    #plt.plot(X, Eb_bathy)
    peak, _ = find_peaks(Eb_bathy)#, prominence=0.01)
    #plt.plot(X[peak[0]], Eb_bathy[peak[0]], 'o')
    #plt.plot(X[peak], Eb_bathy[peak], 'o')
    #plt.show()
    if len(peak) > 0:
        x_breaking = np.max(X[peak]) 
        h_breaking = h[peak][np.argmax(X[peak])]
        #x_breaking = X[peak[0]]
        #h_breaking = h[peak[0]]
        #if x_breaking > -180:
        #    x_breaking = float('NaN')
        #    h_breaking = float('NaN')
    else:
        x_breaking = float('NaN')
        h_breaking = float('NaN')

    return serie, x_breaking, h_breaking, X, Eb_bathy, Hrms


###############################################################################
def preprocessing(readfile):
    x_raw, y_raw, z_raw = indata(readfile, 'nc', '\t')
    Xorig = 70
    # North:
    Yorig = 840
    Ly = 450
    # Full bathy:
    #Yorig = 500
    #Ly = 1140

    Lx = 520
    A = np.zeros(2)
    B = np.zeros(2)
    A[0] = Xorig - Lx
    A[1] = Yorig - Ly / 2
    B[0] = Xorig
    B[1] = Yorig + Ly / 2
    Arot = A.copy()
    Arot[0], Arot[1] = rotation(A[0], A[1], 180, Xorig, Yorig, 1)
    gridsize = 5
    Nx = int(abs(Lx / float(gridsize)))
    Ny = int(Ly / float(gridsize))
    dx = float(gridsize)
    dy = float(gridsize)
    Xmin = A[0] - float(gridsize)
    Ymin = A[1] - float(gridsize)
    Xmax = Xmin + Nx * dx + 2 * float(gridsize)
    Ymax = Ymin + Ny * dy + 2 * float(gridsize)
    val_lim_y_raw = Ymax
    if (Nx % 2) != 0:
        Nx += 1
    if (Ny % 2) != 0:
        Ny += 1
    x = np.linspace(Xmin, Xmax, Nx)
    y = np.linspace(Ymin, Ymax, Ny)
    [X, Y] = np.meshgrid(x, y)
    X, Y = rotation(X, Y, 180, Xorig, Yorig, 1)
    x = X[0, :]
    y = Y[:, 0]
    hraw = np.zeros((Ny, Nx))
    h_interp = interpolation(hraw, 2, Nx, Ny, 'nc', x_raw, y_raw, z_raw,
                             'linear', X, Y)
    h_replaced = replacementValues(0, h_interp, ' ', val_lim_y_raw, 100000,
                                   2, Nx, Ny, 'nc', 'nearest', X, Y, 4, '\t',
                                   'False')
    X, Y = rotation(X, Y, 180, Xorig, Yorig, 0)
    h = h_replaced
    h = -h
    for i in range(np.shape(h)[0]):
        for j in range(np.shape(h)[1]):
            if h[i, j] > 20:
                h[i, j] = float('NaN')

    return X, Y, h


###############################################################################
# Choice and reading of the configuration file
# ONLY DUCK HERE
G = 9.81
RHOW = 1024
VISK = 1.3e-3 / RHOW
with open(r'./duck_configuration.yaml') as file:
    cl = yaml.full_load(file)
    cmapBathy = cl['parameter']['colormaps']
    T12_Y = cl['storm']['T12_Y']
    H12_Y = cl['storm']['H12_Y']
    T_mean = cl['mean']['T_mean']
    H_mean = cl['mean']['H_mean']
    D_50 = float(cl['sand_size']['D50'])
    SRHO = float(cl['sand_size']['rhos'])
    D = D_50 * (G * (SRHO/RHOW - 1) / (VISK**2))**(1.0/3.0)
    WSED = VISK * (np.sqrt(10.36**2 + 1.049*D**3) - 10.36) / D_50
    formula = cl['doc']['formula']
    #HTL_val = cl['sealevel']['HTL']
    #z_max = HTL_val + cl['sealevel']['diff_to_zero_hydro']
    #LTL_val = cl['sealevel']['LTL'] + cl['sealevel']['diff_to_zero_hydro']
    sigma = cl['fit_dissipative']['sigma']

###############################################################################
# Read the forcing serie
serie = pandas.read_csv('./stats_waves_duck_daily_8m_1990-2022.csv')
serie_eq, serie = omega_eq(serie)
# Adapt to the year or years you want
serie_eq['annee'] = serie_eq['year'].dt.year
serie_eq.drop(serie_eq[serie_eq['annee'] < 2000].index, inplace=True)
#serie_eq.drop(serie_eq[serie_eq['annee'] != 2022].index, inplace=True)
serie_eq['year'] = serie_eq['date'].apply(_from_ordinal)
serie_eq = serie_eq.drop(columns=['annee'])
serie_eq['diff'] = (serie_eq['omega'] - serie_eq['omega_eq']) / serie['omega_eq']
serie_eq['mean'] = serie_eq['diff'].rolling(720).mean()  # monthly
#serie_eq = storm_detection(serie_eq)

# Resample to one day
serie_eq = serie_eq.set_index('year')
serie_eq = serie_eq.drop(columns='date')
#serie_eq = serie_eq.resample('1D').mean()
serie_eq = serie_eq.reset_index()
serie_eq['date'] = serie_eq['year'].apply(str)
serie_eq['date'] = serie_eq['date'].apply(datenum)

# Breaking
serie_eq['L0'] = G * serie['Tp']**2 / (2 * np.pi)

# Niveaux de marée
HTL_val = serie['Level'].max()
LTL_val = serie['Level'].min()
z_max = HTL_val

###############################################################################
# ADAPT WITH THE DIRECTORY WITH YOUR BATHYS ###################################
###############################################################################
all_elevation_files = sorted(glob.glob('../bathy/*.nc'))
# Remove problematic bathymetries
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20020407.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20021003.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20021010.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20050917.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20061101.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20080519.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20081109.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20100906.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20100114.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20150710.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20150807.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20150915.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20151008.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20151019.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20151023.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20151030.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20160527.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20181004.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20200417.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20200504.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20200515.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20200527.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20200805.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20200821.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20201228.nc')
all_elevation_files.remove('../bathy/FRF_geomorphology_DEMs_surveyDEM_20210209.nc')

bars = {}
mean_dimensions = pandas.DataFrame(columns=['date', 'mean_height', 'mean_width',
                                            'mean_depth', 'mean_position', 'mean_height_2',
                                            'mean_width_2', 'mean_depth_2', 'mean_position_2',
                                            'frag_1', 'frag_2'])

for elevation_file in all_elevation_files:
    matches = datefinder.find_dates(elevation_file)
    for match in matches:
        date_bathy = match.strftime('%Y-%m-%d')
    print(date_bathy)
    X, Y, h = preprocessing(elevation_file)
    serie_eq, x_breaking, h_breaking, Xthorn, Ebb, Hrms = thornton(serie_eq, X, h)
    bars[date_bathy], nb_bar, fig, gs, = analyse_morpho(X, Y, h, formula, date_bathy, Xthorn, Ebb, Hrms)
    
    if nb_bar > 0 and 'h_bar_1' in bars[date_bathy]:
        mean_height = bars[date_bathy]['h_bar_1'].mean()
        mean_width = bars[date_bathy]['w_bar_1'].mean()
        mean_depth = bars[date_bathy]['z_crest_1'].mean()
        mean_position = bars[date_bathy]['crest_1'].mean()
        mean_frag = bars[date_bathy]['frag_1'].mean()
    else:
        mean_height = float('NaN')
        mean_width = float('NaN')
        mean_depth = float('NaN')
        mean_position = float('NaN')
        mean_frag = float('NaN')
    if nb_bar == 2 and 'h_bar_2' in bars[date_bathy]:
        mean_height_2 = bars[date_bathy]['h_bar_2'].mean()
        mean_width_2 = bars[date_bathy]['w_bar_2'].mean()
        mean_depth_2 = bars[date_bathy]['z_crest_2'].mean()
        mean_position_2 = bars[date_bathy]['crest_2'].mean()
        mean_frag_2 = bars[date_bathy]['frag_2'].mean()
    else:
        mean_height_2 = float('NaN')
        mean_width_2 = float('NaN')
        mean_depth_2 = float('NaN')
        mean_position_2 = float('NaN')
        mean_frag_2 = float('NaN')
    new_row = pandas.Series({'date': date_bathy,
                             'mean_height': mean_height,
                             'mean_width': mean_width,
                             'mean_depth': mean_depth,
                             'mean_position': mean_position,
                             'mean_frag': mean_frag,
                             'mean_height_2': mean_height_2,
                             'mean_width_2': mean_width_2,
                             'mean_depth_2': mean_depth_2,
                             'mean_position_2': mean_position_2,
                             'mean_frag_2': mean_frag_2})
    mean_dimensions = pandas.concat([mean_dimensions,
                                     new_row.to_frame().T],
                                    ignore_index=True)

    if (serie_eq['year'] == date_bathy).any():
        conditions = [(serie_eq['year'] == date_bathy), (serie_eq['year'] != date_bathy)]
        state = ['yes', 'no']
        serie_eq['date_bathy'] = np.select(conditions, state)
        ax5 = fig.add_subplot(gs[0, :])
        sns.scatterplot(data=serie_eq, x='year', y='Hs', ax=ax5, hue='date_bathy',
                        size='date_bathy', sizes=(20, 500), size_order=['yes', 'no'],
                        legend=False)
        ax5.axvline(serie_eq.loc[serie_eq['year'] == date_bathy, 'year'].values,
                    color='orange', linewidth=3)
        plt.title(date_bathy)
        ax5.set_xlabel('Year')
        ax5.set_ylabel('Hs (m)')
        #plt.savefig('./figures_2000-2022_0-180/figure_'+date_bathy+'.png')
        #plt.show()
    plt.close()


mean_dimensions['date'] = mean_dimensions['date'].apply(datenum)
mean_dimensions['time'] = mean_dimensions['date'].apply(_from_ordinal)
mean_dimensions = mean_dimensions.set_index('time')
mean_dimensions = mean_dimensions.drop(columns='date')
#mean_dimensions = mean_dimensions.resample('1H').mean()
mean_dimensions = mean_dimensions.reset_index()
mean_dimensions['date'] = mean_dimensions['time'].apply(str)
mean_dimensions['date'] = mean_dimensions['date'].apply(datenum)
serie_final = serie_eq.merge(mean_dimensions, on='date', how='left')
serie_final['time'] = serie_final['date'].apply(_from_ordinal)


#############
# SAVE DATA #
#############
# TO SAVE THE SHORELINE POSITION, UNCOMMENT shoreline_evolution LINES
# AND MODIFY NAME OF CSV
shoreline_evolution = serie_final[['date', 'shoreline']].copy()
shoreline_evolution = shoreline_evolution[shoreline_evolution['shoreline'].notna()].reset_index()
shoreline_evolution = shoreline_evolution.drop(columns='index')
#shoreline_evolution['diff_shoreline'] = shoreline_evolution.loc[1:, 'shoreline'] - shoreline_evolution.at[0, 'shoreline']
shoreline_evolution = shoreline_evolution.set_index('date')
#shoreline_evolution.to_csv('./shoreline_duck_2000-2022.csv', sep=' ', header=False, columns=['shoreline'])

# TO SAVE THE BAR POSITION, UNCOMMENT evolution_position LINES
# HEIGHT evolution_height
# WIDTH evolution_width
# DEPTH evolution_depth
evolution_height = serie_final[['date', 'mean_height']].copy()
evolution_width = serie_final[['date', 'mean_width']].copy()
#evolution_depth = serie_final[['date', 'mean_depth']].copy()
evolution_position = serie_final[['date', 'mean_position', 'mean_frag']].copy()
evolution_height = evolution_height[evolution_height['mean_height'].notna()].reset_index()
evolution_height = evolution_height.drop(columns='index')
#evolution_height['diff_height'] = evolution_height.loc[1:, 'mean_height'] - evolution_height.at[0, 'mean_height']
evolution_width = evolution_width[evolution_width['mean_width'].notna()].reset_index()
evolution_width = evolution_width.drop(columns='index')
#evolution_width['diff_width'] = evolution_width.loc[1:, 'mean_width'] - evolution_width.at[0, 'mean_width']
#evolution_depth = evolution_depth[evolution_depth['mean_depth'].notna()].reset_index()
#evolution_depth = evolution_depth.drop(columns='index')
#evolution_depth['diff_depth'] = evolution_depth.loc[1:, 'mean_depth'] - evolution_depth.at[0, 'mean_depth']
evolution_position = evolution_position[evolution_position['mean_position'].notna()].reset_index()
evolution_position = evolution_position.drop(columns='index')
evolution_position['diff_position'] = evolution_position.loc[1:, 'mean_position'] - evolution_position.at[0, 'mean_position']
evolution_height = evolution_height.set_index('date')
evolution_width = evolution_width.set_index('date')
#evolution_depth = evolution_depth.set_index('date')
evolution_position = evolution_position.set_index('date')

#evolution_height.to_csv('./evolution_bars_height.csv', header=False, sep=' ', columns=['mean_height'])
#evolution_width.to_csv('./evolution_bars_width.csv', header=False, sep=' ', columns=['mean_width'])
#evolution_depth.to_csv('./evolution_bars_depth.csv', header=False, sep=' ', columns=['diff_depth'])
#evolution_position.to_csv('./evolution_bars_frag_thornton_2000-2022.csv', header=False, sep=' ', columns=['mean_position', 'mean_frag'])
