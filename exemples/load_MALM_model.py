# -*- coding: utf-8 -*-
"""
Created on Thu May 28 08:01:49 2020

@author: Benjamin
"""
import os
import numpy as np 
from fatiando.vis import mpl #, myv
from fatiando import gridder, mesher, utils
from fatiando.gravmag import prism, imaging, transform
from fatiando.vis.mpl import square
import dEXP as dEXP

def load_MALM(ZZ=-3.75,shape=(30,30)):

    # ------------------------------  Model parametes
    ZZ = ZZ # depth of the synthetic anomaly
    x1, x2, y1, y2, z1, z2 = -5,5,-5,5,ZZ-2.5/2,ZZ+2.5/2
    Rsoil = 1000
    
    # ------------------------------- Load data
    filename = 'MSoilR' + str(Rsoil) + 'AnoR1Z'+ str(ZZ) +'L5h2.5'
    MainPath='E:/Padova/Simulation/MALM_SensitivityAnalysis/Sens3d/' + filename + '/Data/'
    #MainPath='E:/Padova/Simulation/MALM_SensitivityAnalysis/Sens3d/' + filename + '/Data/'
    os.chdir(MainPath)
    x,y,z,gz=np.loadtxt('3d_SENS_FAT.txt', unpack=True)
    
    # ------------------------------- remove B return electrode effects 
    B = [-104.16666666666667, -104.16666666666667, 0]
    U = dEXP.cor_field_B(x,y,z,gz,B,rho=100)
    
    # U_cor = U
    xp,yp,U = gridder.interp(x,y,U,shape)
    # xp,yp,gz_cor= gridder.interp(x,y,gz_cor,shape)
    
    return xp, yp, z, U