#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Read the input data
"""

from netCDF4 import Dataset
import numpy as np


def indata(input_path, extension, delimponct):
    """
    Args:
        input_path: str
            file location
        extension: .txt/.xyz/.grd/.nc
            Bathymetry file extension
        delimponct: str
            ponctuation between column (.txt or .xyz)

    Returns:
        x_raw: Array of points according to x direction\n
        y_raw: Array of points according to y direction\n
        z_raw: Array of values for the mesh defined by x_raw and y_raw\n
    """

    # load variable function of extension
    if extension in ("txt", "xyz"):
        x_raw, y_raw, z_raw = np.loadtxt(
            input_path,
            delimiter=delimponct,
            usecols=(0, 1, 2),
            unpack=True,
        )

    if extension == "grd":
        nc_fil = Dataset(input_path)
        x_read = nc_fil.variables["x"][:]
        y_read = nc_fil.variables["y"][:]
        z_raw = nc_fil.variables["z"][:, :]
        [x_raw, y_raw] = np.meshgrid(x_read.data, y_read.data)

    if extension == "nc":
        nc_fil = Dataset(input_path)
        try:
            x_read = nc_fil.variables["xFRF"][:]
            y_read = nc_fil.variables["yFRF"][:]
            z_raw = nc_fil.variables["elevation"][:, :]
            [x_raw, y_raw] = np.meshgrid(x_read.data, y_read.data)
        except KeyError:
            try:
                x_raw = nc_fil.variables["longitude"][:, :]
            except KeyError:
                x_raw = nc_fil.variables["X"][:, :]
            try:
                y_raw = nc_fil.variables["latitude"][:, :]
            except KeyError:
                y_raw = nc_fil.variables["Y"][:, :]
            try:
                z_raw = nc_fil.variables["elevation"][:, :]
            except KeyError:
                z_raw = nc_fil.variables["h"][:, :]

    return x_raw, y_raw, z_raw
