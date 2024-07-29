#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This function can replace the missing values with two options:
    1 - Using another bathymetry
    2 - By approximation of the point

With a big number of missing values, prefer to use other bathymetry data
"""

import sys
from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt
from interpolation_function import interpolation
from indata_function import indata
#sys.path.append("./")


def replacementValues(
    miss_val_param,
    h_bathy,
    missing_bathy_path,
    val_lim_y_raw,
    miss_more,
    partition,
    N_x,
    N_y,
    extension,
    interp_method,
    X,
    Y,
    n,
    delim_ponct_miss_val,
    web,
):

    """
    Parameters
    ----------
    miss_val_param : int
        defined the method to refill the hole on bathymetry.
    h_bathy : array
        bathymetry data.
    missing_bathy_path : str
        path to load the bathymetry used for replaced missing values.
    val_lim_y_raw : int
        inferior limite under we don't calculate/modify values.
    miss_more : int
        to define the method of fill when miss_val_param is choice to use an
        other bathymetry.
    partition : int
        define the number of partitial area to reduce the RAM use.
    N_x : float
        meshgrid number according to x axis.
    N_y : float
        meshgrid number according to y axis.
    extension : str
        extension of the initial file .txt/.xyz/.grd/.nc.
    interp_method : str
        determine the method used for interpolate.
    X : array
        matrix of point position according to x axis
    Y : array
        matrix of point position according to y axis
    n : int
        determine the number of point used in 1 direction to approximate
        the values of the points.
    delim_ponct_miss_val : str
        give the argement between two data row in the file
    web:
        if the bathymetry is on the web
    Returns
    -------
    h_bathy : array
        new bathymetry with values on each point of the grid.
    """

    miss_val = np.argwhere(np.isnan(h_bathy))
    if miss_val_param == 1:
        print("miss_val :", len(miss_val))
        # if missing values exist then complete with large scale bathymetry
        if np.size(miss_val) >= 1:
            print("At least 1 missing value")

            x_raw2, y_raw2, z_raw2 = indata(
                missing_bathy_path, extension, delim_ponct_miss_val, web
            )
            fig = plt.figure(
                num=4, figsize=(10, 10), dpi=100, facecolor="w", edgecolor="w"
            )
            ax = fig.add_subplot()
            if extension == "grd":
                surf = ax.contourf(y_raw2, x_raw2, z_raw2, cmap="gist_earth")
            else:
                surf = ax.scatter(
                    y_raw2,
                    x_raw2,
                    c=z_raw2,
                    vmin=np.min(z_raw2),
                    vmax=0,
                    cmap="gist_earth",
                )

                plt.grid()
            plt.ylim([np.min(x_raw2), np.max(x_raw2)])
            ax.invert_yaxis()
            plt.ylabel("Cross-shore  [m]", size=12)
            plt.xlabel("Long-shore  [m]", size=12)
            plt.colorbar(surf, spacing="uniform", orientation="horizontal")
            plt.title("missing values contour")
            plt.show(block=False)
            z_raw2[np.where(y_raw2 < val_lim_y_raw)] = np.NaN

            if len(miss_val) >= miss_more:
                print("Replacement values by an other Bathymetry:",
                      missing_bathy_path)
                X = x_raw2.reshape(N_y, N_x)
                Y = y_raw2.reshape(N_y, N_x)
                h_bathy = interpolation(
                    np.zeros((N_y, N_x)),
                    partition,
                    N_x,
                    N_y,
                    extension,
                    x_raw2,
                    y_raw2,
                    z_raw2,
                    interp_method,
                    X,
                    Y,
                )
            else:
                print("Replacement ponctual values by other bathymetry:",
                      missing_bathy_path)
                h_bathy[miss_val] = griddata(
                    (x_raw2, y_raw2),
                    z_raw2,
                    (X[miss_val], Y[miss_val]),
                    interp_method,
                )
            h_bathy[h_bathy == 0] = np.NaN
            miss_val = np.argwhere(np.isnan(h_bathy))
            print("missing values stay:", miss_val)
    if miss_val_param == 2:
        print("miss_val :", len(miss_val))
        # # rplc_val == replace value
        rplc_val = 0
        rplc_val_x = 0
        rplc_val_y = 0
        passage = 0
        while len(miss_val) >= 1:
            h_bathy = np.nan_to_num(h_bathy)
            for k in range(len(miss_val)):
                if miss_val[k, 0] > n:
                    if h_bathy[miss_val[k, 0] - (1), miss_val[k, 1]] != 0:
                        for kk in range(n):
                            rplc_val = (
                                rplc_val
                                + h_bathy[miss_val[k, 0] - (1 + kk), miss_val[k, 1]]
                                - h_bathy[miss_val[k, 0] - (2 + kk), miss_val[k, 1]]
                            )

                        rplc_val_x = (
                            rplc_val / n
                            + h_bathy[miss_val[k, 0] - 1, miss_val[k, 1]]
                        )
                        rplc_val = 0
                    elif (
                        miss_val[k, 0] + (1) < len(h_bathy[:, 0])
                        and h_bathy[miss_val[k, 0] + (1), miss_val[k, 1]] != 0
                    ):
                        for kk in range(n):
                            rplc_val = (
                                rplc_val
                                + h_bathy[miss_val[k, 0] + (1 + kk), miss_val[k, 1]]
                                - h_bathy[miss_val[k, 0] + (2 + kk), miss_val[k, 1]]
                            )

                        rplc_val_x = (
                            rplc_val / n
                            + h_bathy[miss_val[k, 0] - 1, miss_val[k, 1]]
                        )
                        rplc_val = 0
                    else:
                        rplc_val_x = "nan"
                        rplc_val = 0
                else:
                    if h_bathy[miss_val[k, 0] + (1), miss_val[k, 1]] != 0:
                        for kk in range(n):
                            rplc_val = (
                                rplc_val
                                + h_bathy[miss_val[k, 0] + (1 + kk), miss_val[k, 1]]
                                - h_bathy[miss_val[k, 0] + (2 + kk), miss_val[k, 1]]
                            )

                        rplc_val_x = (
                            rplc_val / n
                            + h_bathy[miss_val[k, 0] + 1, miss_val[k, 1]]
                        )
                        rplc_val = 0
                    else:
                        rplc_val_x = "nan"
                        rplc_val = 0
                if miss_val[k, 1] > n:
                    if h_bathy[miss_val[k, 0], miss_val[k, 1] - (1)] != 0:
                        for kk in range(n):
                            rplc_val = rplc_val + (
                                h_bathy[miss_val[k, 0], miss_val[k, 1] - (1 + kk)]
                                - h_bathy[miss_val[k, 0], miss_val[k, 1] - (2 + kk)]
                            )
                        rplc_val_y = (
                            rplc_val / n
                            + h_bathy[miss_val[k, 0], miss_val[k, 1] - 1]
                        )
                        rplc_val = 0
                    elif (
                        miss_val[k, 1] + (1) < len(h_bathy[0, :])
                        and h_bathy[miss_val[k, 0], miss_val[k, 1] + (1)] != 0
                    ):
                        for kk in range(n):
                            rplc_val = rplc_val + (
                                h_bathy[miss_val[k, 0], miss_val[k, 1] + (1 + kk)]
                                - h_bathy[miss_val[k, 0], miss_val[k, 1] + (2 + kk)]
                            )

                        rplc_val_y = (
                            rplc_val / n
                            + h_bathy[miss_val[k, 0], miss_val[k, 1] - 1]
                        )
                        rplc_val = 0
                    else:
                        rplc_val_y = "nan"
                        rplc_val = 0
                else:
                    if h_bathy[miss_val[k, 0], miss_val[k, 1] + (1)] != 0:
                        for kk in range(n):
                            rplc_val = rplc_val + (
                                h_bathy[miss_val[k, 0], miss_val[k, 1] + (1 + kk)]
                                - h_bathy[miss_val[k, 0], miss_val[k, 1] + (2 + kk)]
                            )
                        rplc_val_y = (
                            rplc_val / n
                            + h_bathy[miss_val[k, 0], miss_val[k, 1] + 1]
                        )
                        rplc_val = 0
                    else:
                        rplc_val_y = "nan"
                if rplc_val_x == "nan" and rplc_val_y != "nan":
                    rplc_val_f = 0
                elif rplc_val_x != "nan" and rplc_val_y == "nan":
                    rplc_val_f = 0
                elif rplc_val_x != "nan" and rplc_val_y != "nan":
                    rplc_val_f = (rplc_val_x + rplc_val_y) / 2
                elif rplc_val_x == "nan" and rplc_val_y == "nan":
                    rplc_val_f = 0
                h_bathy[miss_val[k, 0], miss_val[k, 1]] = (
                    h_bathy[miss_val[k, 0], miss_val[k, 1]] + rplc_val_f
                )
            rplc_val = 0
            rplc_val_x = 0
            rplc_val_y = 0
            miss_val = np.zeros(())
            h_bathy[h_bathy == 0] = np.NaN
            miss_val = np.argwhere(np.isnan(h_bathy))
            passage += 1
        print(passage)

    return h_bathy
