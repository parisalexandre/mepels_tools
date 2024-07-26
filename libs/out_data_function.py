#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This module contains the function out_data which save bathymetric data in
netCDF format.
"""

import os
import time
import numpy as np
import yaml
from netCDF4 import Dataset


def out_data(
    X,
    Y,
    Nx,
    Ny,
    h,
    filterBathy,
    save_Plot_Bathy,
    fig1,
    fig2,
    fig3,
    configname,
    date_bathy,
):

    """
    This function is for recovers calculated parameters and saved the
    data and figures in the files choice for each.

    Args:
        X : array
            matrix of point position according to x axis
        Y : array
            matrix of point position according to y axis
        Nx: int
            number of mesh on x axis
        Ny: int
            number of mesh on y axis
        h: array
            matrix size (Ny,Nx) with z value for each point on the mesh
        filterBathy: int
            0 or 1
        save_Plot_Bathy: int
            0, 1 or 2
        fig1:
            figure or np.nan
        fig2:
            figure or np.nan
        fig3:
            figure or np.nan
        configname: str
            name of the configation choice
        date_bathy: str
            date of the bathymetry used
    Return:
        Nothing on the code but save all different data
    """

    with open(r"./configurations/" + configname +
              "_configuration.yaml") as file:
        config_file = yaml.full_load(file)

    if filterBathy == 0:
        output_path = (
            config_file["data_path"]["savePath_preproc"]
            + configname
            + date_bathy
            + "_nonFilter"
        )
    else:
        output_path = (
            config_file["data_path"]["savePath_preproc"]
            + configname
            + date_bathy
            + "_Filter"
        )

    print(output_path)

    if save_Plot_Bathy in (0, 2):
        if fig2 != "nan":
            fig2.savefig(output_path + "_gradh.png", dpi=100)
        if fig1 != "nan":
            fig1.savefig(
                output_path + "_bathy.png",
                format="png",
                dpi=100,
            )
        if fig3 != "nan":
            fig3.savefig(output_path + "_bathyLarge.png", dpi=100)

    if save_Plot_Bathy != 0:

        file_name = output_path + "_bathy.nc"

        x = X[0, 0: int(Nx)]
        y = Y[0: int(Ny), 0]

        # open netcdf file
        rootgrp = Dataset(file_name, "w", format="NETCDF4")

        rootgrp.createDimension("x", Nx)
        rootgrp.createDimension("y", Ny)

        var_x = rootgrp.createVariable('X', 'f8', ('y', 'x',))
        var_x.long_name = 'X cross-shore position'
        var_x.units = 'meters'

        var_y = rootgrp.createVariable('Y', 'f8', ('y', 'x'))
        var_y.long_name = 'Y long-shore position'
        var_y.units = 'meters'

        h_var = rootgrp.createVariable('h', 'f8', ('y', 'x',))
        h_var.long_name = 'Bathymetry'
        h_var.units = 'meters'

        var_x[:, :] = X
        var_y[:, :] = Y
        h_var[:, :] = h

        # attributs
        rootgrp.description = "Positionning and remodeling raw bathymetry"
        if config_file["changeStudyArea"]["changeArea"] == 1:
            rootgrp.Xorig_seaward = config_file["changeStudyArea"]["Xorig"]
            rootgrp.Yorig_seaward = config_file["changeStudyArea"]["Yorig"]
            rootgrp.cross_shore_length = config_file["changeStudyArea"]["Lx"]
            rootgrp.long_shore_length = config_file["changeStudyArea"]["Ly"]
            rootgrp.orignie_side = config_file["changeStudyArea"]["origstudy"]
        else:
            print('To save the area definition (Lx, Ly, Xorig, Yorig and \
                  origin side) the parameter "changeArea" need to be \
                  equal to 1 on the configuration file')
        rootgrp.rotation_degrees = config_file["changeRotation"]["axisrotation"]
        rootgrp.rotation_sense_def = "0 : trigo sense / 1 : anti-trigo sense"
        rootgrp.orientation_rotation = config_file["changeRotation"]["trigosense"]

        if config_file["projection"]["before"] == 1:
            rootgrp.info_projection1 = "projection make before interpolation"
            rootgrp.Coordinate_initialy = config_file["projection"]["coorIn_B"]
            rootgrp.Coordinate_finally = config_file["projection"]["coorOut_B"]
        else:
            rootgrp.info_projection1 = "no projection before interpolation"

        if config_file["projection"]["after"] == 1:
            rootgrp.info_projection2 = "projection make after data replacement"
            rootgrp.Coordinate_initialy = config_file["projection"]["coorIn_A"]
            rootgrp.Coordinate_finally = config_file["projection"]["coorOut_A"]
        else:
            rootgrp.info_projection2 = "no projection after data replacement"

        rootgrp.local_landmark_def = (
            "0 : no local coordinate / 1 : local coordinate"
        )
        rootgrp.local_landmark = config_file["projection"]["changeCoordGeo"]

        rootgrp.output_path = config_file["data_path"]["savePath_preproc"]

        rootgrp.history = "Created " + time.ctime(time.time())

        rootgrp.close()

        np.savetxt(
            output_path + ".txt",
            np.column_stack(
                (
                    np.reshape(X, Nx * Ny),
                    np.reshape(Y, Nx * Ny),
                    np.reshape(h, Nx * Ny),
                )
            ),
            delimiter=",",
            newline=os.linesep,
            fmt="%10.5f",
        )

    print("bathymetry save into:", config_file["data_path"]["savePath_preproc"])
