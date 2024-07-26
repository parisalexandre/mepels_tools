#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This function induces rotation to modify position of one point or a vector
of points. The function appears two times, in the beginning and in the end.
"""

import numpy as np


def rotation(X, Y, axisrotation, x_orig, y_orig, trigo_sense):
    """
    Parameters
    ----------
    X : array
        matrix of point position according to x axis
        or just a x point position
    Y : array
        matrix of point position according to y axis
        or just a y point position
    axisrotation : int
        defined the rotation angle in degres.
    x_orig : array
        x position used to induce rotation
    y_orig : array
        y position used to induce rotation
    trigo_sense : int
        define the sense of rotation (trigonometric def)

    Returns
    -------
    X : array
        new rotation matrix/point : position of points according to x axis
    Y : array
        new rotation matrix/point : position of points according to y axis
    """

    x_old = X
    y_old = Y
    if trigo_sense == 0:
        print(
            "Apply rotation in trigonometric direction with angle :",
            axisrotation,
            " degrees",
        )
        axisrotation_r = axisrotation * np.pi / 180.0
        X = (
            x_orig
            + (x_old - x_orig) * np.cos(axisrotation_r)
            - (y_old - y_orig) * np.sin(axisrotation_r)
        )
        Y = (
            y_orig
            + (x_old - x_orig) * np.sin(axisrotation_r)
            + (y_old - y_orig) * np.cos(axisrotation_r)
        )

    else:
        print(
            "Apply rotation in anti-trigonometric direction with angle :",
            axisrotation,
            " degrees",
        )
        axisrotation_r = axisrotation * np.pi / 180.0
        X = (
            x_orig
            + (x_old - x_orig) * np.cos(axisrotation_r)
            + (y_old - y_orig) * np.sin(axisrotation_r)
        )
        Y = (
            y_orig
            - (x_old - x_orig) * np.sin(axisrotation_r)
            + (y_old - y_orig) * np.cos(axisrotation_r)
        )

    return X, Y
