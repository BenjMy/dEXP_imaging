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
import utils_dEXP as uEXP

# icsd functions
from icsd3d.importers.read import load_obs, load_geom

import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 15

import pickle


def load_MALM_sens3d(filename=None):

    file = open(filename, 'rb')
    u = pickle._Unpickler(file)
    u.encoding = 'latin1'
    data = u.load()
    
       
    # dstruct = { "SoilR" : SoilR, "AnoR" : AnoR , "HWD" : [HAno,widthAno,depthAno], "XYU" : uz0_grid,
    #            "shape" : shape, "p12": [p1,p2]}
    # afile = open(SimName + '.pkl', 'wb')
    # pickle.dump(dstruct, afile)
    # afile.close()

    SimName='M' + 'SoilR' + str(data['SoilR']) + 'AnoR' + str(data['AnoR']) + 'Z' + str(data['HWD'][0]) + 'W' + str(data['HWD'][1]) +  'D' + str(data['HWD'][2])

    maxdepth = data['HWD'][2] * 1.5
    shape = data['shape']
    p1 = data['p12'][0]
    p2 = data['p12'][1]
    xyzu = data['XYU']
    xp, yp, z, U  = xyzu[:,0],  xyzu[:,1], xyzu[:,2],  xyzu[:,3]
    
    # ano_prop = data

    return xp, yp, z, U, maxdepth, shape, p1, p2, SimName, data


def load_MALM_synthetic(ZZ=-13.75,shape=(30,30),field=False):

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
    if field==True:
         xp,yp, zp, U  = np.loadtxt(filename + '_uz0_grid.dat', unpack=True)
    else:
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

def definep1p2(path,radius,AnoPosp1p2):
    # print(AnoPosXY)
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
    # plt.figure(figsize=(20,10))
    # ax = plt.subplot()
    # plt.plot(coords_liner[:,0],coords_liner[:,1],'*-')
    # for i in range(len(coords_liner[:,0])):
    #     ax.annotate(str(i), (coords_liner[i,0], coords_liner[i,1]))
    # plt.scatter(B[0],B[1],color='red')
    ptsliner = [2,3]
    
    # p2, p1 = [[coords_liner[ptsliner[0],0],coords_liner[ptsliner[0],1]],
    #           [coords_liner[ptsliner[1],0],coords_liner[ptsliner[1],1]]]
    
    # p1 = [coords_liner[2,0], 5.03629e6] #coords_liner[1,1]]
    # p2 = [coords_liner[4,0],5.03618e6]

    # p1 = [284146.0994252465, 5036293.660496061] #coords_liner[1,1]]
    # p2 = [284369.22176640795,5036158.52301766]
    x1, y1, x2, y2 = squaremat(AnoPosp1p2,r=radius)
    p2 = [x1,y1] #coords_liner[1,1]]
    p1 = [x2,y2]
    
    # print(p1)
    # print(p2)

    # p1 = [284146.0994252465, 5036293.660496061] #coords_liner[1,1]]
    # p2 = [coords_liner[4,0],5036190.0]

    # print(p1)
    # print(p2)


    plt.plot([p1[0],p2[0]],[p1[1],p2[1]],color='red')
    
    return coords_liner, p1, p2, B
    
def creategrid(coords_liner, B, shape, offset=500):
    
    # extrapolation borders points
    xnew_min = coords_liner[2,0] - offset
    ynew_min = B[1] - offset
    xnew_max = coords_liner[4,0] + offset
    ynew_max = coords_liner[2,1] + offset
    
    # create a regular grid
    xnew, ynew = gridder.regular((xnew_min, xnew_max, ynew_min, ynew_max), shape)
    # plt.figure()
    # plt.scatter(xnew,ynew, color='green')
    # plt.show()
    
    return xnew, ynew
    
    
def load_field_u_LandfillPorto(path, filename):

    U_field = np.genfromtxt(path+filename) # load observation data
    x, y, z = np.loadtxt(path + 'coord_zerolevel.dat', unpack=True) # load observation data
    # mesh_xyz = np.array(mesh_xyz)




    # print(mesh_xyz)
    return x, y, z , U_field

    
def load_MALM_LandfillPorto(path, filename, shape=None, field=True, interp=True, radius=20):

    if field==True:
        
        file = open(path + 'PortoM.pkl', 'rb')
        u = pickle._Unpickler(file)
        u.encoding = 'latin1'
        data_struct = u.load()
        
        data_struct['HZ']
        xA = (data_struct['HZ'][0][0]+data_struct['HZ'][0][1])/2
        yA = data_struct['HZ'][0][2]
        AnoPos = [xA,yA]

    
        # x , y, z, U_raw  = load_field_u_LandfillPorto(path, filename + '_uz0.dat')
        print('load' + str(path + filename) + '_uz0_grid.dat')
        x, y, z, U_raw = np.loadtxt(path + filename + '_uz0_grid.dat', unpack=True)
        # x , y, z = mesh_xyz
        # ------------------------------- Plot the data
        # plt.figure()
        # plt.scatter(x, y,c=U_raw, cmap='viridis',vmin=None, vmax=1)
        # plt.colorbar()
        # plt.axis('square')
        # plt.show()
    
    else:         
        U_raw = load_obs(path, filename + '.txt') # load observation data
        RemLineNb, Injection, coordE, pointsE =  load_geom(path) # find the geom file in the folder path
        nb, x , y, z = np.array(coordE[:-3]).T
        AnoPos = []
        data_struct = []

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
    
    # print(AnoPosxAyA)
    coords_liner, p1, p2, B = definep1p2(path,radius,AnoPos)
    xnew, ynew = creategrid(coords_liner, B, shape)
    
    # ------------- Raw data   -----------------------------------------
    
    # pEXP.plot_line(x, y, U_raw,p1,p2, title='U_raw', interp=interp) #, vmin=0.01, vmax=0.1, 
    print(len(x))
    print(len(U_raw))
    print(len(xnew))

    # ------------- Interpolation  -----------------------------------------
    U_int = gridder.interp_at(x, y, U_raw, xnew, ynew, algorithm='cubic', extrapolate=True)
    # xp,yp,U_int = gridder.interp(xnew,ynew,U,shape)
    
    # pEXP.plot_line(xnew, ynew, U_int,p1,p2,title='U_int', interp=interp)
    
    
    # ------------- correction of B  + interpolation  ----------------------
    # remove influence of B
    if field==False:
        Ucor = uEXP.cor_field_B(x,y,z, U_raw, B,rho=10)
    else:
        Ucor = U_raw
    
    Ucor_int = gridder.interp_at(x, y, Ucor, xnew, ynew, algorithm='cubic', extrapolate=True)
    # xp,yp,Ucor_int = gridder.interp(xnew,ynew,Ucorint_tmp,shape)
    # pEXP.plot_line(xnew, ynew, Ucor_int,p1,p2, title='U_cor&int', interp=interp)
    
    # --------------- plot 2d maps -----------------------------------------
    # plt.figure(figsize=(20,10))

    # plt.subplot(2,2,1)
    # plt.tricontourf(x, y, U_raw, 50, cmap='viridis', vmin = 0.01, vmax = 1)
    # cbar = plt.colorbar(orientation='vertical', aspect=50)
    # cbar.set_label('$u_{B}$ (V)')
    # plt.plot(coords_liner[:,0],coords_liner[:,1],'*-')
    # plt.axis('square')
    
    # plt.subplot(2,2,2)
    # plt.tricontourf(x, y, Ucor, 50, cmap='viridis', vmin = 0.01, vmax = 1)
    # cbar = plt.colorbar(orientation='vertical', aspect=50)
    # cbar.set_label('$u_{Bcor}$ (V)')
    # plt.plot(coords_liner[:,0],coords_liner[:,1],'*-')
    # plt.axis('square')
    
    # plt.subplot(2,2,3)
    # plt.tricontourf(xnew, ynew, U_int, 50, cmap='viridis', vmin = 0.01, vmax = 1)
    # cbar = plt.colorbar(orientation='vertical', aspect=50)
    # cbar.set_label('$u_{Binterp}$ (V)')
    # plt.plot(coords_liner[:,0],coords_liner[:,1],'*-')
    # plt.axis('square')
    
    # plt.subplot(2,2,4)
    # plt.tricontourf(xp, yp, Ucor_int, 50, cmap='viridis', vmin = 0.01, vmax = 1)
    # cbar = plt.colorbar(orientation='vertical', aspect=50)
    # cbar.set_label('$u_{Bint&cor}$ (V)')
    # plt.plot(coords_liner[:,0],coords_liner[:,1],'*-')
    # plt.axis('square')
    
    # ---------------------------------------
    zp = 0
    z = zp
    # ---------------------------------------
    
    p = [p1,p2]
    U = [U_raw, Ucor, U_int, Ucor_int]
    coord_xyz = [x,y,z]
    coord_xyz_int = [xnew, ynew,zp]
    
    max_elevation = 50
    
    return coord_xyz, coord_xyz_int, U, coords_liner, shape, max_elevation, p, data_struct



def mirrorImage( a, b, c, x1, y1): 
	temp = -2 * (a * x1 + b * y1 + c) /(a * a + b * b) 
	x = temp * a + x1 
	y = temp * b + y1 
	return (x, y) 


# %% check position of points against line

def isabove(xvcol, yvcol, Ucol, p1,p2):
    


    p = np.vstack([xvcol,yvcol]).T
    isabove = lambda p, p1,p2: np.cross(p-p1, p2-p1) < 0
    
    
    # xvcol = x.reshape([x.shape[0]*y.shape[0],1])
    # yvcol = y.reshape([x.shape[0]*y.shape[0],1])
    # Ucol = U.reshape([x.shape[0]*y.shape[0],1])
    
    c_above=isabove(p,p1,p2)
    
    xvcol_a = np.delete(xvcol , np.where(c_above == True))
    yvcol_a = np.delete(yvcol , np.where(c_above == True))
    U_a = np.delete(Ucol, np.where(c_above == True))
    
    xvcol_b = np.delete(xvcol , np.where(c_above == False))
    yvcol_b = np.delete(yvcol , np.where(c_above == False))
    U_b = np.delete(Ucol, np.where(c_above == False))
    
    p_a= np.vstack([xvcol_a,yvcol_a]).T
    p_a= np.vstack([xvcol_a,yvcol_a]).T

    p_b= np.vstack([xvcol_b,yvcol_b]).T
    p_b= np.vstack([xvcol_b,yvcol_b]).T
    
    plt.figure()
    plt.subplot(1,2,1)
    plt.scatter(xvcol, yvcol,c=Ucol, cmap='viridis', vmax=max(Ucol), vmin=0)
    # plt.scatter(xvcol[c_above==True], yvcol[c_above==True],c='red', cmap='viridis')
    # plt.scatter(xvcol[c_above==False], yvcol[c_above==False],c='black', cmap='viridis')

    p12x=[p1[0],p2[0]]
    p12y=[p1[1],p2[1]]
    
    plt.plot(p12x,p12y, marker="o", color="k")
    plt.axis('square')
    
    plt.subplot(1,2,2)
    plt.scatter(xvcol_a, yvcol_a,c=U_a, cmap='viridis', vmax=max(Ucol), vmin=0)
    plt.colorbar()
    plt.axis('square')

    return U_a, p_a, c_above, U_b, p_b


# %% calcul equation of the line using polyfit

def slope(p1,p2):

    p12x=[p1[0],p2[0]]
    p12y=[p1[1],p2[1]]
    
    # Calculate the coefficients. This line answers the initial question. 
    coefficients = np.polyfit(p12x, p12y, 1)
    
    # Print the findings
    # print ('a =', coefficients[0])
    # print ('b =', coefficients[1])
    
    a = -coefficients[0]
    b = 1
    c = -coefficients[1]
    

    
    return a, b, c

def mirrorU_alongLine(U,p,c_above,a,b,c):
    
   plt.figure()
    
   # Umirror = np.copy(U).tolist()
    # pmirror = np.copy(p).tolist()
   # pmirror = np.ones(p.shape).tolist()
   Umirror= []
   pmirror= []
   for i, bool_pi in enumerate(zip(c_above,p)):
        xmir, ymir = mirrorImage(a, b, c, bool_pi[1][0], bool_pi[1][1]); 
        # plt.scatter(xmir, ymir,c=U[i], cmap='viridis')
        # plt.annotate(str(i)+ '_m', [xmir, ymir])
        Umirror.append(U[i])
        pmirror.append([xmir, ymir])
        
   pmirror = np.vstack(pmirror)
   Umirror = np.array(Umirror)
   plt.scatter(pmirror[:,0],pmirror[:,1],c=Umirror, cmap='viridis',vmax=max(Umirror))
   plt.colorbar()
   plt.axis('square')
    
   return Umirror, pmirror


def squaremat(AnoPosXYmat,r=130):
    # MainPath= 'E:/Padova/Experiments/GC_2019_Landfill_Porto_MALM'
    # os.chdir(MainPath)
    # grid = pv.read('Ano_field_EA.vtk') 
    # grid.set_active_scalars("_Marker")
    # cpos = grid.plot()
    
    # bodies = grid.split_bodies()
    # for i, body in enumerate(bodies):
    #     print('Body {} volume: {:.3f}'.format(i, body.volume)) 
    
    
    # position of the anomaly
    ya =  5.036220e6
    xa =  284267 #(284272.0 + 284250.00)/2 -0.055
    
    # xa = AnoPosXYmat[0] - 0.100
    # ya = AnoPosXYmat[1]

    # landfill geometry
    # utm coordinates 
    coords_liner = [(284046.43,	5036328.39),
              (284132.24,	5036277.32),
              (284146,	5036297),
              (284281.07,	5036214.66),
              (284317.46,	5036281.81),
              (284245,	5036313),
              (284097.55,	5036411.08)]
    coords_liner= np.asarray(coords_liner)
    
    slope = (coords_liner[2,1]-coords_liner[3,1])/(coords_liner[2,0]-coords_liner[3,0]) # (yb-ya)/(xb-xa)
    slope = slope +0.06
    x1 = xa + r*np.cos(slope)
    y1 = ya + r*np.sin(slope)
    
    x2 = xa - r*np.cos(slope)
    y2 = ya - r*np.sin(slope)
    
    plt.figure(figsize=(20,10))
    ax = plt.subplot()
    plt.plot(coords_liner[:,0],coords_liner[:,1],'*-')
    for i in range(len(coords_liner[:,0])):
        ax.annotate(str(i), (coords_liner[i,0], coords_liner[i,1]))
        
    plt.scatter(xa,ya,c='red')
    plt.scatter(x1,y1,c='blue')
    plt.scatter(x2,y2,c='blue')
    plt.axis('square')

    return x1,y1,x2,y2

def zeros_sym(U,c_above,value=0):
    
    Uzeros = np.copy(U)
    ind_True = np.where(c_above==True)[0]
    Uzeros[ind_True] = value
    
    return Uzeros
    

def rotate_and_rescale_all(X_raw,Y_raw,coordE,p1,p2,coords_liner):
    
    # rotate all data
    origin=(max(X_raw), min(Y_raw))
    
    
    point_torotate = np.array([ X_raw, Y_raw])
    Xdraw, Ydraw = rotate_60(origin, point_torotate, angle_deg=rot, showfig=False)
    Xd = Xdraw-min(Xdraw)
    Yd = Ydraw-min(Ydraw)
    
    point_torotate = np.array([coordE[:,1],  coordE[:,2]])
    coordErx, coordEry = rotate_60(origin, point_torotate, angle_deg=rot, showfig=False)
    coordErx = coordErx-min(Xdraw)
    coordEry = coordEry-min(Ydraw)
    
    point_torotate = np.array([[p1[0],p2[0]],[p1[1],p2[1]]])
    px, py = rotate_60(origin,point_torotate,angle_deg=rot, showfig=False)
    px = px-min(Xdraw)
    py = py-min(Ydraw)
    
    p1 = [px[0],py[0]]
    p2 = [px[1],py[1]]
    
    point_torotate = np.array(coords_liner).T
    coords_linerx, coords_linery = rotate_60(origin,point_torotate,angle_deg=rot, showfig=False)
    coords_linerx = coords_linerx-min(Xdraw)
    coords_linery = coords_linery-min(Ydraw)
    
    coords_liner = np.array([coords_linerx, coords_linery]).T

    return Xd, Yd, coordErx,coordEry,coords_liner,p1,p2

    