# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 14:27:34 2020

@author: Benjamin
"""
import numpy as np
# icsd functions
#from importers.read import load_obs, load_geom
import matplotlib.pyplot as plt
from fatiando import gridder, mesher, utils
import lib.utils_dEXP as uEXP
import lib.plot_dEXP as pEXP
import lib.tools.rotate as rotate
import os

def load_geom(path):
    """ load the geometry of the acquisition (*geom file custum for Mise-à-la-masse data)
    
    Parameters
    ----------

    """
    geom_files = [f for f in os.listdir(path) if f.endswith('.geom')]
    if len(geom_files) != 1:
        raise ValueError('should be only one geom file in the current directory')
    
    fileNameElec = geom_files[0]  
    line_number = 0
    line_of_injection = []
    line_of_remotes = []
    # Open the file in read only mode
    with open(path + fileNameElec, 'r') as read_obj:
        # Read all lines in the file one by one
        for line in read_obj:
            # For each line, check if line contains the string
            line_number += 1
            if ('#Remote') in line:
                # If yes, then add the line number & line as a tuple in the list
                line_of_remotes.append((line_number))
            if ('#Injection') in line:
                line_of_injection.append((line_number))
    
    RemLineNb= int(line_of_remotes[0])-1
    Injection= int(line_of_injection[0])-1
    
    coordE = np.loadtxt(path+ fileNameElec)
    pointsE= np.vstack(coordE[:RemLineNb,1:4])
    
    return RemLineNb, Injection, coordE, pointsE


def rotate_60(origin,point_torotate,angle_deg=60, showfig=False):
    
    # X_raw, Y_raw, U_raw_r = np.genfromtxt('malm_models/XYObs_real_f_m3.txt', delimiter='', unpack=True)
        
    
    # point_torotate = np.array([X_raw, Y_raw])
    print(point_torotate)
    angle_rad = np.deg2rad(angle_deg)
    Xr, Yr = rotate.rotate_around_point_lowperf(point=point_torotate,
                                                radians=angle_rad,
                                                origin=origin)
    
    # p1r, p2r = rotate.rotate_around_point_lowperf(point=np.array(p1,p2),
    #                                             radians=angle_rad,
    #                                             origin=(max(X_raw), min(Y_raw)))
    
    if showfig == True:
        plt.figure()
        plt.subplot(1,2,1)
        
        plt.scatter(point_torotate[:,0], point_torotate[:,1])
        plt.grid()
        plt.legend()
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')
        
        plt.subplot(1,2,2)
        plt.scatter(Xr, Yr)
        plt.grid()
        plt.legend()
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')
    
    return Xr, Yr
    


def load_MALM_Porto_real(path, filename,shape,radius,rcor, showfig=False, rot=0):
    
    X_raw, Y_raw, U_raw = np.genfromtxt(filename, delimiter='', unpack=True)
    # U_raw = load_obs(path, filename + '.txt') # load observation data
    RemLineNb, Injection, coordE, pointsE =  load_geom(path) # find the geom file in the folder path
    nb, x , y, z = np.array(coordE[:-3]).T
    AnoPos = []
    
    if shape is None:
        shape = (31,31)
    else:
        shape = shape
    

    
    origin=(max(X_raw), min(Y_raw))

    point_torotate = np.array([ X_raw, Y_raw])
    # print(point_torotate)
    Xd, Yd = rotate_60(origin, point_torotate, angle_deg=rot, showfig=False)
    
    # print(AnoPosxAyA)
    coords_liner, p1, p2, B = definep1p2(path,radius,AnoPos)
    # xnew, ynew = creategrid(coords_liner, B, shape)
    
    point_torotate = np.array([[p1[0],p2[0]],[p1[1],p2[1]]])
    px, py = rotate_60(origin,point_torotate,angle_deg=rot, showfig=False)
    p1 = [px[0],py[0]]
    p2 = [px[1],py[1]]

    # print(coordE)
    point_torotate = np.array(coords_liner).T
    print(point_torotate)

    coords_linerx, coords_linery = rotate_60(origin,point_torotate,angle_deg=rot, showfig=False)
    coords_liner = np.array([coords_linerx, coords_linery]).T
    
    # Xint, Yint, U_int = np.genfromtxt(path + 'XYObs_real_f_m3_grd_surf.txt', delimiter='', unpack=True)
    # # print(len(Xint))

    # point_torotate = coordE
    # coordEx, coordEy = rotate_60(origin,point_torotate,angle_deg=rot, showfig=True)
    # coordE = np.array([coordEx, coordEy])
    
    
    # ------------- Raw data   -----------------------------------------
    
    
    # ------------- Interpolation  -----------------------------------------
    # U_int = gridder.interp_at(Xd, Yd, U_raw, xnew, ynew, algorithm='linear', extrapolate=True)
    # xp,yp,U_int = gridder.interp(xnew,ynew,U,shape)
    
    
    # ------------- correction of B  + interpolation  ----------------------
    # remove influence of B
    Zd = np.zeros(len(Xd))
    Ucor = uEXP.cor_field_B(Xd, Yd, Zd, U_raw, B,rho=rcor,
                            plt_2 = np.array([coords_liner[:,0],coords_liner[:,1]]))


    # Ucor_int = gridder.interp_at(Xd, Yd, Ucor, xnew, ynew, algorithm='linear', extrapolate=True)
    # xp,yp,Ucor_int = gridder.interp(xnew,ynew,Ucorint_tmp,shape)

    # Zd = np.zeros(len(Xint))
    # Uint_cor = uEXP.cor_field_B(Xint, Yint, Zd, U_int, B,rho=rcor)

    if showfig==True:
        pEXP.plot_line(Xd, Yd, U_raw,p1,p2, title='U_raw') #, vmin=0.01, vmax=0.1, 
        # pEXP.plot_line(Xint, Yint, U_int,p1,p2,title='U_int')
        pEXP.plot_line(Xd, Yd, Ucor,p1,p2, title='U_cor')
        # pEXP.plot_line(Xint, Yint, Uint_cor,p1,p2, title='U_cor&int')
        
        # --------------- plot 2d maps -----------------------------------------
        plt.figure(figsize=(20,10))
    
        plt.subplot(2,1,1)
        plt.scatter(Xd, Yd, c=U_raw, cmap='viridis')
        cbar = plt.colorbar(orientation='vertical', aspect=50)
        cbar.set_label('$u_{raw}$ (V)')
        plt.plot(coords_liner[:,0],coords_liner[:,1],'*-')
        plt.axis('square')
        
        plt.subplot(2,1,2)
        plt.scatter(Xd, Yd, c=Ucor, cmap='viridis')
        cbar = plt.colorbar(orientation='vertical', aspect=50)
        cbar.set_label('$u_{Bremoved}$ (V)')
        plt.plot(coords_liner[:,0],coords_liner[:,1],'*-')
        plt.axis('square')
        
        # plt.subplot(2,2,3)
        # plt.scatter(Xint, Yint, c=U_int, cmap='viridis')
        # cbar = plt.colorbar(orientation='vertical', aspect=50)
        # cbar.set_label('$u_{Binterp}$ (V)')
        # plt.plot(coords_liner[:,0],coords_liner[:,1],'*-')
        # plt.axis('square')
        
        # plt.subplot(2,2,4)
        # plt.scatter(Xint, Yint, c=Uint_cor, cmap='viridis')
        # cbar = plt.colorbar(orientation='vertical', aspect=50)
        # cbar.set_label('$u_{Bint&rmv}$ (V)')
        # plt.plot(coords_liner[:,0],coords_liner[:,1],'*-')
        # plt.axis('square')
    
    # ---------------------------------------
    zp = 0
    z = zp
    # ---------------------------------------
    
    p = [p1,p2]
    # U = [U_raw, Ucor, U_int, Uint_cor]
    U = [U_raw, Ucor]

    
    coord_xyz = [Xd, Yd, np.zeros(len(Yd))]
    # coord_xyz_int = [Xint, Yint,zp]
    
    max_elevation = 50

    # return coord_xyz, coord_xyz_int, U, coords_liner, shape, max_elevation, p
    return coord_xyz, coord_xyz, U, coords_liner, shape, max_elevation, p

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
    slope = np.deg2rad(60+90)
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