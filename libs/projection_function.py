#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This function transforms one geographic coordinates system into another
(Ex :WGS to Lambbert-93). The transformation is possible two times, at the
beginning and at the end.
"""

from pyproj import CRS
from pyproj import Transformer


def projection(var_x, var_y, coor_in, coor_out):

    """
    Parameters
    ----------
    var_x : array
        matrix of point position according to x axis
    var_y : array
        matrix of point position according to y axis
    coor_in: str
        EPSG of bathymetry input
    coor_out: str
        EPSG of bathymetry output


    Returns
    -------
    rot_var_X : array
        new projection matrix of point position according to x axis
    rot_var_Y : array
        new projection matrix of point position according to y axis

    """

    crs_input = CRS(coor_in)
    crs_output = CRS(coor_out)
    print("Original projection system :", crs_input.name)
    print("Output projection system :", crs_output.name)

    transformer = Transformer.from_crs(crs_input, crs_output, always_xy=True)
    x_old = var_x
    y_old = var_y
    rot_var_x, rot_var_y = transformer.transform(x_old, y_old, radians=False)

    return rot_var_x, rot_var_y
