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


# my own functions
import dEXP as dEXP
import plot_dEXP as pEXP

# icsd functions
from icsd3d.importers.read import load_obs, load_geom

import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 15



def load_MALM_synthetic(ZZ=-13.75,shape=(30,30)):

    # ------------------------------  Model parametes
    ZZ = ZZ # depth of the synthetic anomaly
    x1, x2, y1, y2, z1, z2 = -5,5,-5,5,ZZ-2.5/2,ZZ+2.5/2
    Rsoil = 1000
    maxdepth=20

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
    
    return xp, yp, z, U, maxdepth, shape



#%% Porto Marghera dataset 

# landfill geometry
# utm coordinates 

def definep1p2(path):

    RemLineNb, Injection, coordE, pointsE =  load_geom(path) # find the geom file in the folder path
    nb, x , y, z = np.array(coordE[:-3]).T
    
    coords_liner = [(284046.43,	5036328.39),
              (284132.24,	5036277.32),
              (284146,	5036297),
              (284281.07,	5036214.66),
              (284317.46,	5036281.81),
              (284245,	5036313),
              (284097.55,	5036411.08)]
    coords_liner= np.asarray(coords_liner)
    

    
    # ------------- define p1 and p2   -----------------------------------------
    B = coordE[RemLineNb+2,1:] # position of the return electrode
    plt.figure(figsize=(20,10))
    ax = plt.subplot()
    plt.plot(coords_liner[:,0],coords_liner[:,1],'*-')
    for i in range(len(coords_liner[:,0])):
        ax.annotate(str(i), (coords_liner[i,0], coords_liner[i,1]))
    plt.scatter(B[0],B[1],color='red')
    ptsliner = [2,3]
    
    # p2, p1 = [[coords_liner[ptsliner[0],0],coords_liner[ptsliner[0],1]],
    #           [coords_liner[ptsliner[1],0],coords_liner[ptsliner[1],1]]]
    
    # p1 = [coords_liner[2,0], 5.03629e6] #coords_liner[1,1]]
    # p2 = [coords_liner[4,0],5.03618e6]

    p1 = [284146.0994252465, 5036293.660496061] #coords_liner[1,1]]
    p2 = [284369.22176640795,5036158.52301766]
    
    plt.plot([p1[0],p2[0]],[p1[1],p2[1]],color='red')
    
    return coords_liner, p1, p2, B
    
def creategrid(coords_liner, B, shape):
    
    # extrapolation borders points
    xnew_min = coords_liner[2,0]
    ynew_min = B[1] 
    xnew_max = coords_liner[4,0] 
    ynew_max = coords_liner[2,1] 
    
    # create a regular grid
    xnew, ynew = gridder.regular((xnew_min, xnew_max, ynew_min, ynew_max), shape)
    # plt.figure()
    plt.scatter(xnew,ynew, color='green')
    plt.show()
    
    return xnew, ynew
    
    
def load_field_u_LandfillPorto(path, filename):

    U_field = np.genfromtxt(path+filename) # load observation data
    x, y, z = np.loadtxt(path + 'coord_zerolevel.dat', unpack=True) # load observation data
    # mesh_xyz = np.array(mesh_xyz)




    # print(mesh_xyz)
    return x, y, z , U_field

    
def load_MALM_LandfillPorto(path, filename, shape=None, field=True, interp=True):

    if field==True:
        # x , y, z, U_raw  = load_field_u_LandfillPorto(path, filename + '_uz0.dat')
        print('load' + str(path + filename) + '_uz0_grid.dat')
        x, y, z, U_raw = np.loadtxt(path + filename + '_uz0_grid.dat', unpack=True)
        # x , y, z = mesh_xyz
        # ------------------------------- Plot the data
        plt.figure()
        plt.scatter(x, y,c=U_raw, cmap='viridis',vmin=None, vmax=1)
        plt.colorbar()
        plt.axis('square')
        plt.show()
    
    else:         
        U_raw = load_obs(path, filename + '.txt') # load observation data
        RemLineNb, Injection, coordE, pointsE =  load_geom(path) # find the geom file in the folder path
        nb, x , y, z = np.array(coordE[:-3]).T

    # path2files="example_2add_later/Landfill_3d/Ano_0_E13/" # Test elec outside landfill
    # # path2files="example_2add_later/Landfill_3d/Ano_0_E89/" # Test 2 elec outside landfill
    # # path2files="example_2add_later/Landfill_3d/Ano_0_EA/"  # Real A position without anomaly
    # # path2files="example_2add_later/Landfill_3d/Ano_1_BH_EA/"  # Real A position with big hole anomaly
    # path2files="example_2add_later/Landfill_3d/Ano_1_SH_EA/"  # Real A position  with small hole anomaly
    # # path2files="examples/Landfill_3d/RealData_Ano_0/"  # Real data / without anomaly
    # # path2files="examples/Landfill_3d/RealData_Ano_1_SH/"  # Real data / with (small) anomaly
    # -----------------------------------------------------------------------------------------
    
    if shape is None:
        shape = (31,31)
    else:
        shape = shape
    
    coords_liner, p1, p2, B = definep1p2(path)
    xnew, ynew = creategrid(coords_liner, B, shape)
    
    # ------------- Raw data   -----------------------------------------
    
    pEXP.plot_line(x, y, U_raw,p1,p2, title='U_raw', interp=interp) #, vmin=0.01, vmax=0.1, 
    
    # ------------- Interpolation  -----------------------------------------
    U = gridder.interp_at(x, y, U_raw, xnew, ynew, algorithm='cubic', extrapolate=True)
    xp,yp,U_int = gridder.interp(xnew,ynew,U,shape)
    pEXP.plot_line(xp, yp, U_int,p1,p2,title='U_int', interp=interp)
    
    
    # ------------- correction of B  + interpolation  ----------------------
    # remove influence of B
    if field==False:
        Ucor = dEXP.cor_field_B(x,y,z, U_raw, B,rho=10)
    else:
        Ucor = U_raw
    
    Ucorint_tmp = gridder.interp_at(x, y, Ucor, xnew, ynew, algorithm='cubic', extrapolate=True)
    xp,yp,Ucor_int = gridder.interp(xnew,ynew,Ucorint_tmp,shape)
    pEXP.plot_line(xp, yp, Ucor_int,p1,p2, title='U_cor&int', interp=interp)
    
    # --------------- plot 2d maps -----------------------------------------
    plt.figure(figsize=(20,10))

    plt.subplot(2,2,1)
    plt.tricontourf(x, y, U_raw, 50, cmap='viridis', vmin = 0.01, vmax = 1)
    cbar = plt.colorbar(orientation='vertical', aspect=50)
    cbar.set_label('$u_{B}$ (V)')
    plt.plot(coords_liner[:,0],coords_liner[:,1],'*-')
    plt.axis('square')
    
    plt.subplot(2,2,2)
    plt.tricontourf(x, y, Ucor, 50, cmap='viridis', vmin = 0.01, vmax = 1)
    cbar = plt.colorbar(orientation='vertical', aspect=50)
    cbar.set_label('$u_{Bcor}$ (V)')
    plt.plot(coords_liner[:,0],coords_liner[:,1],'*-')
    plt.axis('square')
    
    plt.subplot(2,2,3)
    plt.tricontourf(xnew, ynew, U_int, 50, cmap='viridis', vmin = 0.01, vmax = 1)
    cbar = plt.colorbar(orientation='vertical', aspect=50)
    cbar.set_label('$u_{Binterp}$ (V)')
    plt.plot(coords_liner[:,0],coords_liner[:,1],'*-')
    plt.axis('square')
    
    plt.subplot(2,2,4)
    plt.tricontourf(xp, yp, Ucor_int, 50, cmap='viridis', vmin = 0.01, vmax = 1)
    cbar = plt.colorbar(orientation='vertical', aspect=50)
    cbar.set_label('$u_{Bint&cor}$ (V)')
    plt.plot(coords_liner[:,0],coords_liner[:,1],'*-')
    plt.axis('square')
    
    # ---------------------------------------
    zp = 0
    z = zp
    # ---------------------------------------
    
    p = [p1,p2]
    U = [U_raw, Ucor, U_int, Ucor_int]
    coord_xyz = [x,y,z]
    coord_xyz_int = [xp,yp,zp]
    
    max_elevation = 20

    return coord_xyz, coord_xyz_int, U, coords_liner, shape, max_elevation, p

