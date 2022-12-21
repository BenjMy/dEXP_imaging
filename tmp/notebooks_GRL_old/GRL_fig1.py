# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 14:51:17 2020

@author: Benjamin
"""
import os
# from importers.read import load_obs, load_geom
import numpy as np
#import lib.tools.rotate as rotate
import math
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})
from fatiando import gridder, mesher, utils

import lib.utils_dEXP as uEXP
import lib.plot_dEXP as pEXP

import notebooks_JGR.loadmalm as MALM


# MainPath= 'E:/Padova/Software/SourceInversion/Potential_field_imaging/dEXP_imaging/examples_in_prep/'
MainPath =  '/home/ben/OneDrive/Padova/Software/SourceInversion/Potential_field_imaging/dEXP_imaging/examples_in_prep/'

os.chdir(MainPath)

rot=60
radius=200
# position of the anomaly
ya =  5.036220e6
xa =  284267 #(284272.0 + 284250.00)/2 -0.055
offset = 255
# def plot_scatter3d():
#     return Xr, Yr
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


def rotate_around_point_lowperf(point, radians, origin=(0, 0)):
    """Rotate a point around a given point.
    
    I call this the "low performance" version since it's recalculating
    the same values more than once [cos(radians), sin(radians), x-ox, y-oy).
    It's more readable than the next function, though.
    """
    x, y = point
    ox, oy = origin

    qx = ox + math.cos(radians) * (x - ox) + math.sin(radians) * (y - oy)
    qy = oy + -math.sin(radians) * (x - ox) + math.cos(radians) * (y - oy)

    return qx, qy


def rotate_and_rescale_all(X_raw,Y_raw,coordE,p1,p2,coords_liner,B):
    
    # rotate all data
    origin=(max(X_raw), min(Y_raw))
    
    
    point_torotate = np.array([ X_raw, Y_raw])
    Xdraw, Ydraw = rotate_60(origin, point_torotate, angle_deg=rot, showfig=False)
    Xd = Xdraw-min(Xdraw)+offset
    Yd = Ydraw-min(Ydraw)+offset
    
    point_torotate = np.array([coordE[:,1],  coordE[:,2]])
    coordErx, coordEry = rotate_60(origin, point_torotate, angle_deg=rot, showfig=False)
    coordErx = coordErx-min(Xdraw) +offset
    coordEry = coordEry-min(Ydraw) +offset
    
    point_torotate = np.array([[p1[0],p2[0]],[p1[1],p2[1]]])
    px, py = rotate_60(origin,point_torotate,angle_deg=rot, showfig=False)
    px = px-min(Xdraw)
    py = py-min(Ydraw)
    
    p1 = [px[0]+offset,py[0]+offset]
    p2 = [px[1]+offset,py[1]+offset]
    
    point_torotate = np.array(coords_liner).T
    coords_linerx, coords_linery = rotate_60(origin,point_torotate,angle_deg=rot, showfig=True)
    coords_linerx = coords_linerx-min(Xdraw)
    coords_linery = coords_linery-min(Ydraw)

    point_torotate = B
    Bx, By = rotate_60(origin,point_torotate,angle_deg=rot, showfig=False)
    Bx = Bx-min(Xdraw)+offset
    By = By-min(Ydraw)+offset
    B = np.c_[Bx,By]
    
    coords_liner = np.array([coords_linerx, coords_linery]).T

    return Xd, Yd, coordErx,coordEry,coords_liner,p1,p2, B


def rotate_60(origin,point_torotate,angle_deg=60, showfig=False):
    
    angle_rad = np.deg2rad(angle_deg)
    Xr, Yr = rotate_around_point_lowperf(point=point_torotate,
                                                radians=angle_rad,
                                                origin=origin)    
    return Xr, Yr
    

def definep1p2(path,radius,AnoPosp1p2):

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
    ptsliner = [2,3]
    x1, y1, x2, y2 = squaremat(AnoPosp1p2,r=radius)
    p2 = [x1,y1] #coords_liner[1,1]]
    p1 = [x2,y2]
    plt.plot([p1[0],p2[0]],[p1[1],p2[1]],color='red')
    
    return coords_liner, p1, p2, B

def squaremat(AnoPosXYmat,r=130):
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

#-----------------------------------------------------------------------------
# read data

PathElecs = MainPath + '/malm_models/'
filename = './malm_models/XYObs_real_f_m3.txt'
X_raw, Y_raw, U_raw = np.genfromtxt(filename, delimiter='', unpack=True)

#E:\Padova\Redaction\Articles\1b_InversionUsingGravityPotMethod\notebooks\data\phNO

# path2files = 'data/ph_low_top_2m/'
path2files = '/home/ben/Documents/GitHub/BenjMy/dEXP_imaging/notebooks_JGR/data/phNO/'
filename= 'NoAno'
# # path2files = 'E:/Padova/Experiments/GC_2019_Landfill_Porto_MALM/Test/ph_low/'
#filename= 'Ano'
X_raw, Y_raw, z, U_raw = np.loadtxt(path2files + filename + '_uz0_grid.dat', unpack=True)


# read elecs to plot
RemLineNb, Injection, coordE, pointsE =  load_geom(PathElecs) # find the geom file in the folder path

# define 2d profile
AnoPos = []
coords_liner, p1, p2, B = definep1p2(PathElecs,radius,AnoPos)


X_raw = X_raw + offset 
Y_raw = Y_raw + offset 
coordE = coordE + offset 
Xd, Yd, coordErx,coordEry,coords_liner,p1,p2,B = rotate_and_rescale_all(X_raw,Y_raw,coordE,p1,p2,coords_liner,B[0:2])


# ------------------------------
fig = plt.figure()
ax=fig.gca(projection='3d')
sc=ax.scatter(coordErx[RemLineNb], coordEry[RemLineNb], -9 ,color='red', s=50, marker=(5, 2))
sc=ax.scatter(coordErx[RemLineNb+1], coordEry[RemLineNb+1], 0 ,color='green', s=50)
sc=ax.scatter(coordErx[Injection], coordEry[Injection],  0,color='green', s=50)
sc=ax.scatter(Xd , Yd, 0, color='black')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_zlim(-25,0)
ax.plot(coords_liner[2:5,0],coords_liner[2:5,1],'k')
ax.plot(coords_liner[2:5,0],coords_liner[2:5,1],-9,'b')
ax.plot(coords_liner[2:5,0],coords_liner[2:5,1],-9-7,'b',)
ax.plot(coords_liner[2:5,0],coords_liner[2:5,1],-25,'k')
l0 = np.c_[[coordErx[RemLineNb+1],  coordEry[RemLineNb+1],0],
           [coordErx[Injection],  coordEry[Injection],0]].T
ax.plot(l0[:,0],l0[:,1],l0[:,2],'g--')
ax.plot([p1[0],p2[0]],[p1[1],p2[1]],color='m')


#%%
fig, ax = plt.subplots()
sc=ax.scatter(coordErx[RemLineNb], coordEry[RemLineNb],
              color='red', 
              s=50, marker=(5, 2),label='MALM electrode')
sc=ax.scatter(coordErx[RemLineNb+1], coordEry[RemLineNb+1],
              color='green', s=50,label='remotes B and N')
sc=ax.scatter(coordErx[Injection], coordEry[Injection], color='green', s=50)
sc=ax.scatter(Xd , Yd, color='black',
              label='surface potential')
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.plot(coords_liner[2:5,0]+605,coords_liner[2:5,1]+165,'k',
              label='liner')
l0 = np.c_[[coordErx[RemLineNb+1],  coordEry[RemLineNb+1],0],
           [coordErx[Injection],  coordEry[Injection],0]].T
ax.plot(l0[:,0],l0[:,1],'g--')
ax.plot([p1[0]+349.71,p2[0]+349.71],[p1[1],p2[1]-200],color='m',
              label='symetry plan')
#ax.axis('equal')
# ax.set_xlim([200,420])
leg = ax.legend(loc='upper left', frameon=False);
# plt.savefig('fig2.png',
#          dpi=450)
plt.show()
np.savetxt('elecs_ref.txt',np.array([coordErx,coordEry]).T)

#%%
# ------------------------------

pEXP.plot_line(Xd, Yd, U_raw,p1,p2, title='U_raw') #, vmin=0.01, vmax=0.1, 
# pEXP.plot_line(Xint, Yint, U_int,p1,p2,title='U_int')
# pEXP.plot_line(Xd, Yd, Ucor,p1,p2, title='U_cor')
# pEXP.plot_line(Xint, Yint, Uint_cor,p1,p2, title='U_cor&int')

# --------------- plot 2d maps -----------------------------------------
plt.figure()
# plt.subplot(2,1,1)
plt.scatter(Xd, Yd, c=U_raw, cmap='viridis', vmin=0, vmax=0.005)
cbar = plt.colorbar(orientation='vertical', aspect=50)
cbar.set_label('$u_{raw}$ (V)')
plt.plot(coords_liner[2:5,0],coords_liner[2:5,1],'k')
# plt.
plt.xlim(min(Xd),max(Xd))
plt.ylim(min(Yd),max(Yd))
# plt.axis('square')

# Xg, Yg = np.meshgrid(Xd, Yd)
# Zg = griddata((Xd, Yd), U_raw, Xg, Yg, method='linear')

prl = 60
# shape = shape  (max(xp)-min(xp))/
shape = (150,150)

xint_scipy, yint_scipy = gridder.regular((min(Xd)-prl, max(Xd)+prl, 
                          min(Yd)-prl, max(Yd)+prl),shape=shape)

U_int_scipy = gridder.interp_at(Xd, Yd, U_raw, xint_scipy, yint_scipy, algorithm='linear', extrapolate=False)
InterpData = np.array([xint_scipy, yint_scipy, U_int_scipy]).T
where_are_NaNs = np.isnan(InterpData)
InterpData[where_are_NaNs] = 0.0074
xint_scipy, yint_scipy, U_int_scipy = InterpData.T

ax, plt = pEXP.plot_field(xint_scipy, yint_scipy,U_int_scipy, shape,Vminmax=[0,0.012])
# ax, plt = pEXP.plot_field(xint_scipy,yint_scipy,U_mirror_int, shape,Vminmax=[0,0.35])

ax.plot(coords_liner[2:5,0],coords_liner[2:5,1],'k')
plt.axis('square')
plt.xlim(min(xint_scipy),max(xint_scipy))
plt.ylim(min(yint_scipy),max(yint_scipy))
plt.xlim(300,500)
plt.ylim(300,500)
#%% Solution 2

rcor=1
B = np.hstack(B)
B = np.r_[B,0]
Zd = np.zeros(len(Xd))
Ucor = uEXP.cor_field_B(Xd, Yd, Zd, U_raw, B,rho=rcor,
                        plt_2 = np.array([coords_liner[:,0],coords_liner[:,1]]))


Extrapolate = True
U_int_scipy = gridder.interp_at(Xd, Yd, Ucor, xint_scipy, yint_scipy, 
                                algorithm='cubic', extrapolate=True)
ax, plt = pEXP.plot_field(xint_scipy,yint_scipy,U_int_scipy, shape,Vminmax=[0,0.009])
ax.plot(coords_liner[2:5,0],coords_liner[2:5,1],'k')
plt.axis('square')
plt.xlim(min(xint_scipy),max(xint_scipy))
plt.ylim(min(yint_scipy),max(yint_scipy))
# plt.xlim(300,500)
# plt.ylim(300,500)
# plt.savefig('fig_publi_Noano_synth.png', dpi=450)
plt.savefig('fig_publi_Cor_field.png', dpi=450)

#%%
# find point position with respect to line equation defined by p1 and p2
U_a, p_a, bool_above, U_b, p_b = MALM.isabove(xint_scipy,yint_scipy, U_int_scipy, 
                                  np.array(p1),np.array(p2))
# Compute p1 and p2 line equation ax + by + c = 0
a, b, c = MALM.slope(p1,p2)
# a, b, c = MALM.slope(np.array([p1[0],p2[0]]),np.array([p1[1],p2[1]]))

# Mirror points with respect to p1p2 line
Umirror, xy_mirror = MALM.mirrorU_alongLine(U_a,p_a,bool_above,a,b,c)

plt.scatter(xy_mirror[:,0],xy_mirror[:,1])
# plt.scatter(xy_mirror[:,0],xy_mirror[:,1])


U_a_int = gridder.interp_at(xy_mirror[:,0], xy_mirror[:,1], Umirror, xint_scipy,yint_scipy, algorithm='cubic', 
                        extrapolate=True)   
# U_mirror_int = np.copy(U_a_int)
U_mirror_int = np.copy(U_int_scipy)
U_mirror_int[np.where(bool_above == True)[0]]= U_a_int[np.where(bool_above == True)[0]]

# ax, plt = pEXP.plot_field(xint_scipy,yint_scipy,U_mirror_int, shape,Vminmax=[0,0.009])
ax, plt = pEXP.plot_field(xint_scipy,yint_scipy,U_mirror_int, shape,Vminmax=[0,0.35])
ax.plot(coords_liner[2:5,0],coords_liner[2:5,1],'k')
plt.axis('square')
plt.xlim(min(xint_scipy),max(xint_scipy))
plt.ylim(min(yint_scipy),max(yint_scipy))
plt.xlim(300,500)
plt.ylim(300,500)
# plt.savefig('fig_publi_Noano_synth.png', dpi=450)
plt.savefig('fig_publi_Cor_field.png', dpi=450)