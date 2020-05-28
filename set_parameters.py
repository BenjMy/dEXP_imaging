# -*- coding: utf-8 -*-
"""
Created on Thu May 28 07:57:02 2020

@author: Benjamin
"""


def set_par(shape,max_elevation,**kwargs):
    # # -------------------------------  Graphical parameters
    scaled=0 
   
    for key, value in kwargs.items():
        if key == 'SI': # data interpolation
            SI = value        
    shape = shape # data interpolation
    qorder = 0 # derivative order of the continuated field

    # # -------------------------------  Imaging parameters
    SI = 2 # default structural index
    zp=0  # initial depth (conditionned upward or downward)

    # ---- z-discretisation - Upward continuation parameters parameters
    max_elevation=max_elevation
    nlay = 25 # discretisation of upward continuation
    
    # ----------- ridges analysis
    minAlt_ridge = max_elevation*0.25
    maxAlt_ridge = max_elevation*0.75
    
    for key, value in kwargs.items():
        if key == 'minAlt_ridge': # data interpolation
            minAlt_ridge = value  
        if key == 'minAlt_ridge': # data interpolation
            minAlt_ridge = value  


    
    return scaled, shape, SI, zp, qorder, max_elevation, nlay, minAlt_ridge, maxAlt_ridge