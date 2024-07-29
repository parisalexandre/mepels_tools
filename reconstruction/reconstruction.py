#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Reconstruction of a cross-shore profile

@author: Alexandre Paris
"""

import time
import seaborn as sns
import numpy as np
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

# BARS DIMENSIONS
height = 0.97
pos = -154 # positive values going offshore
width = 90
# SHORELINE POSITION
shoreline = -63 # positive values going offshore
# ADAPT THE HTL VALUE
HTL = 1.87

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
peaks, _ = find_peaks(bar_profile, prominence=0.01)
troughs, _ = find_peaks(-bar_profile, prominence=0.01)
while bar_profile[peaks]-bar_profile[troughs] < height:
    a = a + 0.05
    bar[:] = gaussian(x_htl, pos, sigma_bar, a/2) + gaussian(x_htl, pos+b, sigma_bar, -a/2)
    bar_profile = z_mean_htl + bar
    peaks, _ = find_peaks(bar_profile, prominence=0.01)
    troughs, _ = find_peaks(-bar_profile, prominence=0.01)

# Same for width
while np.abs(x_htl[peaks]-x_htl[troughs]) < width:
    b = b + 0.1
    bar[:] = gaussian(x_htl, pos, sigma_bar, a/2) + gaussian(x_htl, pos+b, sigma_bar, -a/2)
    bar_profile = z_mean_htl + bar
    peaks, _ = find_peaks(bar_profile, prominence=0.01)
    troughs, _ = find_peaks(-bar_profile, prominence=0.01)

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
np.savetxt('./bed_w{}_06122019.dep'.format(width), bar_profile, fmt='%1.8f', newline=' ')

# Load measured profile
zmeas = np.loadtxt('./north_06122019.dep')

# Plot
fig, ax = plt.subplots()
ax.plot(x, z_mean, label='mean profile', linestyle=':', linewidth=4, color='green')
ax.plot(x, zmeas, label='measured profile', color='black', linestyle='--')
ax.plot(new_x, bar_profile, label='reconstructed profile', color='red')
ax.set_ylim([-7.5, 7.5])
ax.set_xlabel('Cross-shore (m)')
ax.set_ylabel('Elevation (m)')
ax.set_title('2019-12-06')
ax.axhline(HTL, color='blue', linestyle='--', label='HTL')
ax.legend(loc='upper left')
plt.tight_layout()
plt.show()
