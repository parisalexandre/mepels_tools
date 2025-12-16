#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Reconstruction of a cross-shore profile

@author: Alexandre Paris
"""

import time
import seaborn as sns
import numpy as np
import pandas
from matplotlib import pyplot as plt
from scipy import interpolate
from scipy.signal import find_peaks

start_time = time.time()
sns.set_theme()
sns.set_context('talk')

###############################################################################
def gaussian(x, mu, sigma, a):
    return a * np.exp((-1 / 2) * ((x - mu) / sigma) ** 2)


def fill_nan(a):
    inds = np.arange(a.shape[0])
    good = np.where(np.isfinite(a))
    f = interpolate.interp1d(inds[good], a[good], bounds_error=False)
    b = np.where(np.isfinite(a), a, f(inds))
    return b


###############################################################################
mean_profile = np.loadtxt('mean_profile_duck_2000-2022.txt')
x = mean_profile[:,0]
z_mean = mean_profile[:,1]
z_mean = z_mean[np.where(~np.isnan(z_mean))]
x = x[np.where(~np.isnan(z_mean))]
deltax = np.mean(np.diff(x))

#HS DEPENDING ON MONTH AT DUCK
hs_month = {'January': 0.89,
            'February': 0.91,
            'March': 0.94,
            'April': 0.83,
            'May': 0.75,
            'June': 0.59,
            'July': 0.52,
            'August': 0.64,
            'September': 1.01,
            'October': 0.97,
            'November': 0.95,
            'December': 0.89}

# Get values
values = np.array(list(hs_month.values()))
# Find min, max, and interval
vmin, vmax = values.min(), values.max()
interval = (vmax - vmin) / 3
# Define bins
low_max = vmin + interval
med_max = vmin + 2*interval

# READ BARS HEIGHT AND WIDTH
heights = pandas.read_csv('/Users/alexandreparis/SHOREFOR/data/barres/evolution_bars_height.csv',
                          names=['time', 'height'], delimiter=' ')
widths = pandas.read_csv('/Users/alexandreparis/SHOREFOR/data/barres/evolution_bars_width.csv',
                         names=['time', 'width'], delimiter=' ')

h_winter = heights['height'].quantile(0.9)
h_summer = heights['height'].quantile(0.1)
h_spring = heights['height'].mean() #same for fall
w_winter = widths['width'].quantile(0.9)
w_summer = widths['width'].quantile(0.1)
w_spring = widths['width'].mean() #same for fall

#DICTIONARY WITH SHOREFOR RESULTS
results = {'date': ['02-15', '03-25', '04-11', '04-17', '05-20', '06-25', '07-16', '08-15', '09-03', '09-09', '09-17', '09-23', '10-15', '10-25', '11-14', '11-22', '12-06'],
           'hs_7': [0.70, 0.93, 0.66, 0.68, 0.63, 0.41, 0.47, 0.52, 0.85, 1.1, 0.94, 1.5, 1.8, 0.87, 1.1, 1.9, 0.96],
#           'height': [h_winter, h_winter, h_spring, h_spring, h_spring, h_summer, h_summer, h_summer, h_winter, h_winter, h_winter, h_winter, h_winter, h_winter, h_winter, h_winter, h_winter],
#           'width': [w_winter, w_winter, w_spring, w_spring, w_spring, w_summer, w_summer, w_summer, w_winter, w_winter, w_winter, w_winter, w_winter, w_winter, w_winter, w_winter, w_winter],
           'shoreline': [-61, -61, -61, -61, -60, -59, -58, -57, -57, -58, -58, -59, -61, -61, -61, -63, -63],
           'pos': [-152, -149, -155, -150, -127, -97, -77, -48, -48, -60, -65, -78, -112, -109, -115, -156, -153]}

## BARS DIMENSIONS
#height = 0.97
#pos = -154 # positive values going offshore
#width = 90
## SHORELINE POSITION
#shoreline = -63 # positive values going offshore

# ADAPT THE HTL VALUE
HTL = 1.87

for i, d in enumerate(results['date']):
    #i = i + 13
    if results['hs_7'][i] < low_max:
        height = h_summer
        width = w_summer
    elif low_max <= results['hs_7'][i] < med_max:
        height = h_spring
        width = w_spring
    elif results['hs_7'][i] >= med_max:
        height = h_winter
        width = w_winter  
    #height = results['height'][i]
    #width = results['width'][i]
    pos = results['pos'][i]
    shoreline = results['shoreline'][i]
    day = results['date'][i][0:2]
    month = results['date'][i][3:5]
    newdate = f'{day}{month}2019'

    a = height*2
    b = width
    sigma_bar = width / 3
    
    z_mean_htl = z_mean[z_mean < HTL]
    x_htl = x[z_mean < HTL]
    
    # Bar creation
    bar = np.zeros(len(x_htl))
    bar[:] = gaussian(x_htl, pos, sigma_bar, a/2) + gaussian(x_htl, pos+b, sigma_bar, -a/2)
    bar_profile = z_mean_htl + bar
    
    # Height adjusment, increasing the Gaussian height
    peaks, _ = find_peaks(bar_profile, prominence=0.1)
    troughs, _ = find_peaks(-bar_profile, prominence=0.1)
    while bar_profile[peaks]-bar_profile[troughs] < height:
        a = a + 0.05
        bar[:] = gaussian(x_htl, pos, sigma_bar, a/2) + gaussian(x_htl, pos+b, sigma_bar, -a/2)
        bar_profile = z_mean_htl + bar
        peaks, _ = find_peaks(bar_profile, prominence=0.1)
        troughs, _ = find_peaks(-bar_profile, prominence=0.1)
    
    # Same for width
    #if x_htl[troughs] > 0 or np.abs(x_htl[troughs]) < width/2:
    if np.abs(x_htl[troughs]) < width/2:
        b = width / 2
        bar[:] = gaussian(x_htl, pos, sigma_bar, a/2) + gaussian(x_htl, pos+b, sigma_bar, -a/2)
        bar_profile = z_mean_htl + bar
    else:        
        while np.abs(x_htl[peaks] - x_htl[troughs]) < width:
            b = b + 0.1
            bar[:] = gaussian(x_htl, pos, sigma_bar, a/2) + gaussian(x_htl, pos+b, sigma_bar, -a/2)
            bar_profile = z_mean_htl + bar
            peaks, _ = find_peaks(bar_profile, prominence=0.1)
            troughs, _ = find_peaks(-bar_profile, prominence=0.1)
    
    ## Definition of new x
    # Look for x>0 (arbitrary)
    new_x_pos = x_htl[x_htl >= 0]
    # Shortest distance between the mean profile last point under HTL and the new shoreline
    diff = -x[z_mean < HTL][-1] - shoreline
    # Extract the part of the new profile which is after x>0 (arbitrary)
    bar_prof = bar_profile[x_htl >= 0]
    # HTL is imposed at the end of the new profile
    bar_prof[-1] = HTL
    # Recalculate the x
    for i, v in enumerate(new_x_pos):
        if i != len(new_x_pos)-1 and i!= 0:
            new_x_pos[i] = float('NaN')
            bar_prof[i] = float('NaN')
        else:
            if diff > 0:
                new_x_pos[i] = v + diff / (len(new_x_pos) - i)
            else:
                new_x_pos[i] = v + diff / (len(new_x_pos) - i)
    
    # Interpolation between 0 et HTL
    new_x_pos = fill_nan(new_x_pos)
    bar_prof = fill_nan(bar_prof)
    bar_profile = np.append(bar_profile[x_htl < 0], bar_prof)
    
    new_x = np.append(x_htl[x_htl < 0], new_x_pos)
    delta_new_x = np.mean(np.diff(new_x))
    bar_profile = bar_profile[0:len(new_x)]
    
    # Save the new profile
    np.savetxt(f'./bed_w{width}_{newdate}.dep', bar_profile, fmt='%1.8f', newline=' ')
    np.savetxt(f'./x_{newdate}.dep', new_x, fmt='%1.8f', newline=' ')
    
    # Load measured profile
    zmeas = np.loadtxt(f'./profiles/north_{month}{day}.dep')
    
    # Plot
    fig, ax = plt.subplots()
    ax.plot(x, z_mean, label='mean profile', linestyle=':', linewidth=4, color='green')
    ax.plot(x, zmeas, label='measured profile', color='black', linestyle='--')
    ax.plot(new_x, bar_profile, label='reconstructed profile', color='red')
    ax.set_ylim([-7.5, 4])
    ax.set_xlabel('Cross-shore (m)')
    ax.set_ylabel('Elevation (m)')
    ax.set_title(f'2019-{day}-{month}', y=0.9)
    ax.axhline(HTL, color='blue', linestyle='--', label='HTL')
    ax.legend(loc='upper left')
    #plt.tight_layout()
    fig.set_tight_layout(True)
    plt.show()
