# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 09:13:31 2020

@author: Benjamin
"""


import harmonica as hm
import verde as vd
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np


import dEXP as dEXP
from dEXP import _fit
import plot_dEXP as pEXP
import utils_dEXP as uEXP
from fatiando.gravmag import prism, imaging, transform


# Define the coordinates for two point masses
# easting = [5e3, 15e3]
# northing = [7e3, 13e3]
easting = [7e3]
northing = [10e3]

# The vertical coordinate is positive upward so negative numbers represent
# depth
# upward = [-0.5e3, -1e3]
upward = [-1e3]
points = [easting, northing, upward]
# We're using "negative" masses to represent a "mass deficit" since we assume
# measurements are gravity disturbances, not actual gravity values.
# masses = [3e11, -10e11]
masses = [3e11]
max_elevation = np.abs(upward)*3
nlay = 25
qorder= 1
zp = 0
interp = True
smooth = False
shape_grid = (200,200)
# Define computation points on a grid at 500m above the ground
coordinates = vd.grid_coordinates(
    region=[0, 20e3, 0, 20e3], shape=shape_grid, extra_coords=500
)

# Compute the downward component of the gravitational acceleration
gravity = hm.point_mass_gravity(
    coordinates, points, masses, field="g_z", coordinate_system="cartesian"
)
#%%


#%%
print(gravity)
# Plot the results on a map
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_aspect("equal")
# Get the maximum absolute value so we can center the colorbar on zero
maxabs = vd.maxabs(gravity)
img = ax.contourf(
    *coordinates[:2], gravity, 60, vmin=-maxabs, vmax=maxabs, cmap="seismic"
)
plt.colorbar(img, ax=ax, pad=0.04, shrink=0.73, label="mGal")
# Plot the point mass locations
ax.plot(easting, northing, "oy")
ax.set_title("Gravitational acceleration (downward)")
# Convert axes units to km
ax.set_xticklabels(ax.get_xticks() * 1e-3)
ax.set_yticklabels(ax.get_yticks() * 1e-3)
ax.set_xlabel("km")
ax.set_ylabel("km")
plt.tight_layout()
plt.show()

yp,xp = coordinates[:2]
xp = np.hstack(np.reshape(xp,[1,xp.shape[0]**2]))
yp = np.hstack(np.reshape(yp,[1,yp.shape[0]**2]))
gravity = np.hstack(np.reshape(gravity,[1,gravity.shape[0]**2]))
    
for slicedir in enumerate('x'):
    if slicedir[1]=='y':
        print(slicedir[1])
        x_axis = slicedir[1] # direction of p1p2 profil
        SI = 1.5
        p1 =[0,(easting[0])]
        p2 =[20e3,(easting[0])]
    else:
       p1 =[(northing[0]),0]
       p2 =[(northing[0]),20e3]
       x_axis = slicedir[1] # direction of p1p2 profil
       SI = 2
        
    pEXP.plot_line(xp, yp, gravity,p1,p2, interp=True, Xaxis=x_axis)
    pEXP.plot_field(xp, yp, gravity, shape_grid)

    #%%
    # Upward continuation of the field data with discretisation in altitude controlled by the number of layers (nlay) and the maximum elevation desired (max_elevation)
    mesh, label_prop = dEXP.upwc(xp, yp, zp, gravity, shape_grid, 
                     zmin=0, zmax=max_elevation, nlayers=nlay, 
                     qorder=qorder)

    plt, cmap = pEXP.plot_xy(mesh, label=label_prop, Xaxis=x_axis, p1p2=np.array([p1, p2]))
    plt.colorbar(cmap)

    #%%

    xderiv = transform.derivx(xp, yp, gravity, shape_grid,order=qorder)
    yderiv = transform.derivy(xp, yp, gravity, shape_grid,order=qorder)
    zderiv = transform.derivz(xp, yp, gravity, shape_grid,order=qorder)
    
    # # plt.plot(xderiv)
    # # plt.plot(yderiv)
    
    # # interp = True
    pEXP.plot_line(xp, yp, xderiv ,p1,p2,title='xderiv',savefig=False, interp=interp, Xaxis=x_axis)
    pEXP.plot_line(xp, yp, yderiv ,p1,p2,title='yderiv',savefig=False, interp=interp, Xaxis=x_axis)
    
    # # p1_perp,p2_perp = uEXP.perp_p1p2(p1,p2, offset=0)
    # # pEXP.plot_line(xp, yp, yderiv ,p1_perp,p2_perp,title='yderiv',savefig=False, interp=interp, Xaxis=x_axis)    
    pEXP.plot_line(xp, yp, zderiv ,p1,p2,title='zderiv',savefig=False, interp=interp, Xaxis=x_axis)


    #%%
    # Ridges identification: plot all extremas obtained via find_peaks function (numpy) for a given altitude
    dEXP.ridges_minmax_plot(xp, yp, mesh, p1, p2,
                                          label=label_prop,
                                          method_peak='find_peaks',
                                          showfig=True,
                                          interp=True,smooth=True,
                                          Xaxis=x_axis)  

    #%%
    # Ridges identification at all levels: plot extremas obtained via find_peaks function (numpy) for all 3 types of extremas familly RI, RII and RIII
    D  = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,
                                          label=label_prop,
                                          method_peak='find_peaks',
                                          qorder=qorder,
                                          interp=True,smooth=True,
                                          fix_peak_nb=1,
                                          returnAmp=True,
                                          showfig=True,
                                          Xaxis=x_axis)  
    dfI, dfII, dfIII =  D[0:3]
    hI, hII, hIII  = D[3:6]
    H  = D[3:6]
    
    #%%
    # filter ridges using a minimum length criterium and and filter for a specific range of altitude
    # a =2.25
    # if x_axis=='y':
    #     xf_min = a*x1
    #     xf_max = a*x2
    # else:   
    #     xf_min = a*y1
    #     xf_max = a*y2        

    D_f = dEXP.filter_ridges(dfI,dfII,dfIII,
                                minDepth=500,
                                maxDepth=5000,
                                minlength=8,
                                rmvNaN=True,
                                xmin=1000, xmax=19000,
                                heights=[hI, hII, hIII])
    
    dfI_f, dfII_f, dfIII_f =  D_f[0:3]
    hI_f, hII_f, hIII_f = D_f[3:6]
    df_f = D_f[0:3]

    #%%
    # plot filtered ridges fitted over continuated section
        
    fig = plt.figure()
    ax = plt.gca()
    
    plt, cmap = pEXP.plot_xy(mesh, label=label_prop, Xaxis=x_axis, p1p2=np.array([p1, p2]), ax=ax) #, ldg=)
    plt.colorbar(cmap)
    pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax,label=True)
    df_fit = dEXP.fit_ridges(df_f, rmvOutliers=False) # fit ridges on filtered data
    pEXP.plot_ridges_sources(df_fit, ax=ax, z_max_source=-max_elevation*1.2,
                              ridge_type=[0,1,2],ridge_nb=None)
    if x_axis=='x':
        plt.scatter(easting, upward,marker=(5, 1),c='red')
        plt.annotate(str(masses),[easting, upward])
    else:   
        plt.scatter(northing, upward,marker=(5, 1),c='red')
        plt.annotate(str(masses),[northing, upward])

    #%% 

    qratio = [1,0]
    mesh_dexp, label_dexp = dEXP.dEXP_ratio(xp, yp, zp, gravity, shape_grid, 
                     zmin=0, zmax=max_elevation, nlayers=nlay, 
                     qorders=qratio)
    fig = plt.figure()
    ax = plt.gca()
    
    plt, cmap = pEXP.plot_xy(mesh_dexp, label=label_dexp,
                 markerMax=True,qratio=str(qratio),
                 p1p2=np.array([p1,p2]), ax=ax, Xaxis=x_axis) #, ldg=)
    plt.colorbar(cmap)

    if x_axis=='x':
        plt.scatter(easting, np.abs(upward),marker=(5, 1),c='red')
        plt.annotate(str(masses),[easting, upward])
    else:   
        plt.scatter(northing, upward,marker=(5, 1),c='red')
        plt.annotate(str(masses),[northing, upward])
        
        
    #%% 
    #  ridges analysis: scaling function to determine the SI index
    # RI and RII : zeros of the first horizontal and first vertical derivatives of the potential field
    # RIII :zeros of the potential field itself
    
    z0 = upward # choose an estimate of the depth of the anomaly
    # z0 = -100 # choose an estimate of the depth of the anomaly

    ncol = 0
    for r_type in range(len(df_f)): # loop over ridges type I, II, III
        if df_f[r_type].shape[1]>1:
            ncol = ncol + df_f[r_type].shape[1]-1


    fig, axs = plt.subplots(1,ncol, figsize=(15, 6), facecolor='w', edgecolor='k')
    # fig.subplots_adjust(hspace = .5, wspace=.001)
    axs = axs.ravel()
   
    # for i in range(ncol):
    
    nc = 0
    SI_est = []
    for r_type in range(len(df_f)): # loop over ridges type I, II, III
        for k in enumerate(df_f[r_type].columns[1:]): # loop over ridges of the same familly
            # df_f[r_type].columns[1:]    
            points, fit, SI_est_tmp , EXTnb = dEXP.scalFUN(df_f[r_type],EXTnb=k[1],z0=z0,
                                                           rmvOutliers=True)
            pEXP.plot_scalFUN(points, fit, z0=z0, 
                              label='R'+ str(r_type+1) + k[1], 
                              ax=axs[nc]) # scaling 
            SI_est.append(SI_est_tmp)
            
            nc = nc + 1
        
    SI_mean = np.mean(np.abs(SI_est)) 
    
    #%% 
    #  ridges analysis
    SI_mean = 0
    mesh_dexp, label_dexp = dEXP.dEXP(xp, yp, zp, gravity, shape_grid, 
                     zmin=0, zmax=max_elevation, nlayers=nlay, 
                     qorder=qorder,
                     SI=SI_mean)
    
    fig = plt.figure()
    ax = plt.gca()
    
    plt, cmap = pEXP.plot_xy(mesh_dexp, label=label_dexp,
                 markerMax=True,SI=SI_mean,
                 p1p2=np.array([p1,p2]), ax=ax, Xaxis=x_axis) #, ldg=)
    plt.colorbar(cmap)

    if x_axis=='x':
        plt.scatter(easting, np.abs(upward),marker=(5, 1),c='red')
        plt.annotate(str(masses),[easting, upward])
    else:   
        plt.scatter(northing, upward,marker=(5, 1),c='red')
        plt.annotate(str(masses),[northing, upward])

    #%% 

uEXP.multipage(dataname + '_test_DEXP_SI_punctual.pdf')
