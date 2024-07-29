#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This function can partition the maps without limit (warning: the real
partition is partition^2) and uses the z_raw row to make interpolation
on meshgrid
"""

from scipy.interpolate import griddata


def interpolation(
    hraw, partition, n_x, n_y, extension, x_raw, y_raw, z_raw, methoduse, X, Y
):

    """
    Parameters
    ----------
    hraw : array
        matrix size (n_x,n_y) empty.
    partition : int
        1, 2, 3 or 4.
    n_x: int
        number of point present on meshgrid follow x axis
    n_y: int
        number of point present on meshgrid follow y axis
    extension: str
        .txt/.xyz/.grd
    x_raw: array
        Vector or Matrix for position according to x axis
    y_raw: array
        Vector or Matrix for position according to y axis
    z_raw: array
        Vector or Matrix for position according to z axis
    methoduse: str
        nearest, linear, quadratic or cubique
    X: array
        meshgrid matrix acording to x axis
    Y: array
        meshgrid matrix according to y axis

    Returns
    -------
    hraw : array
        matrix hraw fill by value (bathymetry) on each point of the mesh
        defined by X_raw and Y_raw
    """

    fillx = n_x // partition
    filly = n_y // partition
    endfillx = n_x - (fillx * partition)
    endfilly = n_y - (filly * partition)
    # fillx: calcul the number of point present in one subdivisions on x axis
    # filly: same as fillx but for y axis
    # endfillx: Using "//" made a round on the inferior value, so we need to
    # re-add some value to obtaint the good dimension of mesh at the end.
    # The adding values are used on the last partition on x axis
    # endfilly: same than endfillx but for y axis

    if extension in ("grd", "nc"):
        x_raw = x_raw.flatten()
        y_raw = y_raw.flatten()
        z_raw = z_raw.flatten()

    print("Interpolate bathymetry on grid")
    k = 0
    m = 0
    block = 1
    while k < (partition):
        while m < (partition):
            # print("Block number", block, '/', partition**2)
            # for X meshgrid matrix
            # TypeError: only integer scalar arrays can be converted to a scalar index
            if k < (partition - 1) and m < (partition - 1):
                a = int(m * filly)
                b = int((m + 1) * filly)
                c = int(k * fillx)
                d = int((k + 1) * fillx)
                X1 = X[
                    a:b,
                    c:d
                    # int(m * filly) : int((m + 1) * filly), int(k * fillx) : int((k + 1) * fillx)
                ]

            elif k < (partition - 1) and m == (partition - 1):
                X1 = X[
                    m * filly: (m + 1) * filly + endfilly,
                    k * fillx: (k + 1) * fillx,
                ]

            elif k == (partition - 1) and m < (partition - 1):
                X1 = X[
                    m * filly: (m + 1) * filly,
                    k * fillx: (k + 1) * fillx + endfillx,
                ]

            elif k == (partition - 1) and m == (partition - 1):
                X1 = X[
                    m * filly: (m + 1) * filly + endfilly,
                    k * fillx: (k + 1) * fillx + endfillx,
                ]
            # for Y meshgrid matrix
            Y1 = []
            if k < (partition - 1) and m < (partition - 1):
                Y1 = Y[m * filly: (m + 1) * filly, k * fillx: (k + 1) * fillx]

            elif k < (partition - 1) and m == (partition - 1):
                Y1 = Y[
                    m * filly: (m + 1) * filly + endfilly,
                    k * fillx: (k + 1) * fillx,
                ]

            elif k == (partition - 1) and m < (partition - 1):
                Y1 = Y[
                    m * filly: (m + 1) * filly,
                    k * fillx: (k + 1) * fillx + endfillx,
                ]

            elif k == (partition - 1) and m == (partition - 1):
                Y1 = Y[
                    m * filly: (m + 1) * filly + endfilly,
                    k * fillx: (k + 1) * fillx + endfillx,
                ]

            hraw1 = griddata((x_raw, y_raw), z_raw, (X1, Y1), method=methoduse)

            # to fill the new matix h
            if m < (partition - 1) and k < (partition - 1):
                hraw[
                    m * (filly): (m + 1) * (filly),
                    k * (fillx): (k + 1) * (fillx),
                ] = (
                    hraw[
                        m * (filly): (m + 1) * (filly),
                        k * (fillx): (k + 1) * (fillx),
                    ]
                    + hraw1
                )
            elif m == (partition - 1) and k < (partition - 1):
                hraw[
                    m * (filly): (m + 1) * (filly) + endfilly,
                    k * (fillx): (k + 1) * (fillx),
                ] = (
                    hraw[
                        m * (filly): (m + 1) * (filly) + endfilly,
                        k * (fillx): (k + 1) * (fillx),
                    ]
                    + hraw1
                )
            if m < (partition - 1) and k == (partition - 1):
                hraw[
                    m * (filly): (m + 1) * (filly),
                    k * (fillx): (k + 1) * (fillx) + endfillx,
                ] = (
                    hraw[
                        m * (filly): (m + 1) * (filly),
                        k * (fillx): (k + 1) * (fillx + endfillx),
                    ]
                    + hraw1
                )
            elif m == (partition - 1) and k == (partition - 1):
                hraw[
                    m * (filly): (m + 1) * (filly) + endfilly,
                    k * (fillx): (k + 1) * (fillx) + endfillx,
                ] = (
                    hraw[
                        m * (filly): (m + 1) * (filly) + endfilly,
                        k * (fillx): (k + 1) * (fillx) + endfillx,
                    ]
                    + hraw1
                )
            m += 1
            # print("interpolation on block", block, "done")
            block += 1
        k += 1
        m = 0

    print("Interpolation finished with success!")
    return hraw
