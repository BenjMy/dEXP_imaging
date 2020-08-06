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
easting = [5e3]
northing = [7e3]

# The vertical coordinate is positive upward so negative numbers represent
# depth
# upward = [-0.5e3, -1e3]
upward = [-0.5e3]
points = [easting, northing, upward]
# We're using "negative" masses to represent a "mass deficit" since we assume
# measurements are gravity disturbances, not actual gravity values.
# masses = [3e11, -10e11]
masses = [3e11]

shape_grid = (100,100)
# Define computation points on a grid at 500m above the ground
coordinates = vd.grid_coordinates(
    region=[0, 20e3, 0, 20e3], shape=shape_grid, extra_coords=500
)

# Compute the downward component of the gravitational acceleration
gravity = hm.point_mass_gravity(
    coordinates, points, masses, field="g_z", coordinate_system="cartesian"
)
#%%

max_elevation = np.abs(upward)*3
nlay = 25
qorder= 0
zp = 500
interp = True
smooth = False
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

xp,yp = coordinates[:2]
xp = np.hstack(np.reshape(xp,[1,xp.shape[0]**2]))
yp = np.hstack(np.reshape(yp,[1,yp.shape[0]**2]))
gravity = np.hstack(np.reshape(gravity,[1,gravity.shape[0]**2]))
    
for slicedir in enumerate('x'):
    if slicedir[1]=='y':
        print(slicedir[1])
        x_axis = slicedir[1] # direction of p1p2 profil
        SI = 1.5
        p1 =[0,(northing[0])]
        p2 =[20e3,(northing[0])]
    else:
       p1 =[(easting[0]),0]
       p2 =[(easting[0]),20e3]
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
                                          Xaxis=x_axis)  

    #%%
    # Ridges identification at all levels: plot extremas obtained via find_peaks function (numpy) for all 3 types of extremas familly RI, RII and RIII
    dfI,dfII, dfIII = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,
                                          label=label_prop,
                                          method_peak='find_peaks',
                                          qorder=qorder,
                                          fix_peak_nb=2,
                                          showfig=True,
                                          Xaxis=x_axis)  

    #%%
    # filter ridges using a minimum length criterium and and filter for a specific range of altitude
    # a =2.25
    # if x_axis=='y':
    #     xf_min = a*x1
    #     xf_max = a*x2
    # else:   
    #     xf_min = a*y1
    #     xf_max = a*y2        

    dfI_f,dfII_f, dfIII_f = dEXP.filter_ridges(dfI,dfII,dfIII,
                                                minDepth=100,
                                                maxDepth=3000,
                                                minlength=8,
                                                rmvNaN=True,
                                                xmin=0, xmax=20e3)
    df_f = dfI_f, dfII_f, dfIII_f

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
    if x_axis=='y':
        plt.scatter(easting, upward,marker=(5, 1),c='red')
        plt.annotate(str(masses),[easting, upward])
    else:   
        plt.scatter(northing, upward,marker=(5, 1),c='red')
        plt.annotate(str(masses),[northing, upward])

    #%% 
    #  ridges analysis: scaling function to determine the SI index
    # RI and RII : zeros of the first horizontal and first vertical derivatives of the potential field
    # RIII :zeros of the potential field itself
    
    z0 = -(z1 + z2)/2 # choose an estimate of the depth of the anomaly
    # z0 = -100 # choose an estimate of the depth of the anomaly

    ncol = 0
    for r_type in range(len(df_f)): # loop over ridges type I, II, III
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
        
    SI_mean = np.mean(SI_est)  
    
    #%% 
    #  ridges analysis
    # SI = 1.5
    mesh_dexp, label_dexp = dEXP.dEXP(xp, yp, zp, U, shape, 
                     zmin=0, zmax=max_elevation, nlayers=nlay, 
                     qorder=qorder,
                     SI=SI_mean)
    
    fig = plt.figure()
    ax = plt.gca()
    
    plt, cmap = pEXP.plot_xy(mesh_dexp, label=label_dexp,
                 markerMax=True,SI=SI_mean,
                 p1p2=np.array([p1,p2]), ax=ax, Xaxis=x_axis) #, ldg=)
    plt.colorbar(cmap)

    if x_axis=='y':
        square([x1, x2, z1, z2])
        plt.annotate(dens,[(x1 + x2)/2, -(z1+z2)/2])
    else:   
        square([y1, y2, z1, z2])
        plt.annotate(dens,[(y1 + y2)/2, -(z1+z2)/2])



uEXP.multipage(dataname + '_test_DEXP_SI_punctual.pdf')
