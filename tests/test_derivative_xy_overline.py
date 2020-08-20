# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 17:55:57 2020

@author: Benjamin
"""

import matplotlib.pyplot as plt
import numpy as np
import dEXP as dEXP
from dEXP import _fit
import plot_dEXP as pEXP
import set_parameters as para
from fatiando.gravmag import prism, imaging, transform
from fatiando.vis.mpl import square

# exemples
import examples.magnetic.fwdmag.fwd_mag_sphere as magfwd
import examples.gravimetry.loadgrav.grav_models as grav
import utils_dEXP as uEXP

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 15


#%% 
# load model previously generated using Fatiando a terra package
# asymetric square: za3000_zb3500_l500_ofs3000_dens1200
# asymetric rectangle: mod_asy_rect.pkl or mod_asy_rect_z1.2km or mod_asy_rect_z3.2km
# symetric:  za_1000zb_1500dens_1200 # or sym_square_z3.2km
# a/symetric with a rectangle space around sym_square_z3.2km_rectArea / mod_asy_rect_z3.2km_rectArea
# za3000_zb3500_l500_ofs0_dens1200

dataname = 'mod_asy_rect_z1.2km'
data_struct = grav.load_grav_fatiando(name='../examples/gravimetry/loadgrav/' + dataname)
xp,yp,zp,U = data_struct['xyzg']
shape = data_struct['shape']
model = data_struct['model']
dens  = data_struct['density']
# scaled, SI, zp, qorder, nlay, minAlt_ridge, maxAlt_ridge = para.set_par(shape=shape,max_elevation=max_elevation)

x1, x2, y1, y2, z1, z2 = np.array(model[0].get_bounds())

# p1 =[min(yp),0]
# p2 =[max(yp),0]

max_elevation=z2*1.2
scaled, SI, zp, qorder, nlay, minAlt_ridge, maxAlt_ridge = para.set_par(shape=shape,max_elevation=max_elevation)
interp = True
qorder = 0


for slicedir in enumerate('xy'):
    if slicedir[1]=='x':
        print(slicedir[1])
        p1 =[(x1+x2)/2,min(yp)]
        p2 =[(x1+x2)/2,max(yp)]
        x_axis = slicedir[1] # direction of p1p2 profil
        SI = 1.5
    else:
        p1 =[min(xp),(y1+y2)/2]
        p2 =[max(xp),(y1+y2)/2]
        x_axis = slicedir[1] # direction of p1p2 profil
        SI = 2

    #%%
    # Plot field data over a 2d line crossing the anomalies
    # pEXP.plot_line(xp, yp, U,p1,p2, interp=False, Xaxis='x')
    
    pEXP.plot_line(xp, yp, U,p1,p2, interp=True, Xaxis=x_axis)
    square([x1, x2, y1, y2])

    pEXP.plot_field(xp,yp,U, shape)
    square([x1, x2, y1, y2])

    #%% ------------------------------- Plot the derivatives
    # http://campar.in.tum.de/Chair/HaukeHeibelGaussianDerivatives
    
    xderiv = transform.derivx(xp, yp, U, shape,order=qorder)
    yderiv = transform.derivy(xp, yp, U, shape,order=qorder)
    zderiv = transform.derivz(xp, yp, U, shape,order=qorder)
    
    # # plt.plot(xderiv)
    # # plt.plot(yderiv)
    
    # # interp = True
    pEXP.plot_line(xp, yp, xderiv ,p1,p2,title='xderiv',savefig=False, interp=interp, Xaxis=x_axis)
    pEXP.plot_line(xp, yp, yderiv ,p1,p2,title='yderiv',savefig=False, interp=interp, Xaxis=x_axis)
    
    # # p1_perp,p2_perp = uEXP.perp_p1p2(p1,p2, offset=0)
    # # pEXP.plot_line(xp, yp, yderiv ,p1_perp,p2_perp,title='yderiv',savefig=False, interp=interp, Xaxis=x_axis)    
    pEXP.plot_line(xp, yp, zderiv ,p1,p2,title='zderiv',savefig=False, interp=interp, Xaxis=x_axis)

    #%%
    # Upward continuation of the field data with discretisation in altitude controlled by the number of layers (nlay) and the maximum elevation desired (max_elevation)
    mesh, label_prop = dEXP.upwc(xp, yp, zp, U, shape, 
                     zmin=0, zmax=max_elevation, nlayers=nlay, 
                     qorder=qorder)

    plt, cmap = pEXP.plot_xy(mesh, label=label_prop, Xaxis=x_axis, p1p2=np.array([p1, p2]))
    plt.colorbar(cmap)

# plt, cmap = pEXP.plot_xy(mesh, label=label_prop, Xaxis='y', p1p2=np.array([p1, p2])) #, ldg=)
# plt.colorbar(cmap)

# plt, cmap = pEXP.plot_xy(mesh, label=label_prop, Xaxis='y', p1p2=p) #, ldg=)
# plt.colorbar(cmap)


# plt, cmap = pEXP.slice_mesh(xp, yp, mesh, label_prop, p1, p2, interp=True, Xaxis='x')
# plt.colorbar(cmap)


    #%%
    # Ridges identification: plot all extremas obtained via find_peaks function (numpy) for a given altitude
    dEXP.ridges_minmax_plot(xp, yp, mesh, p1, p2,
                                          label=label_prop,
                                          method_peak='find_peaks',
                                          showfig=True,
                                          Xaxis=x_axis,
                                          qorder=qorder)  

    #%%
    # Ridges identification at all levels: plot extremas obtained via find_peaks function (numpy) for all 3 types of extremas familly RI, RII and RIII
    D = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,
                                          label=label_prop,
                                          method_peak='find_peaks',
                                          fix_peak_nb=3,
                                          returnAmp=True,
                                          showfig=True,
                                          Xaxis=x_axis,
                                          interp=interp,
                                          qorder=qorder)  

    dfI, dfII, dfIII =  D[0:3]
    hI, hII, hIII  = D[3:6]
    heights  = D[3:6]
    #%%
    # plot filtered ridges fitted over continuated section
        
    fig = plt.figure()
    ax = plt.gca()
    
    plt, cmap = pEXP.plot_xy(mesh, label=label_prop, Xaxis=x_axis, p1p2=np.array([p1, p2]), ax=ax) #, ldg=)
    plt.colorbar(cmap)
    pEXP.plot_ridges_harmonic(dfI, dfII, dfIII,ax=ax,label=True)
    # df_fit = dEXP.fit_ridges(D[0:3], rmvOutliers=True) # fit ridges on filtered data
    # pEXP.plot_ridges_sources(df_fit, ax=ax, z_max_source=-max_elevation*1.2,
    #                           ridge_type=[0,1,2],ridge_nb=None)
    if x_axis=='y':
        square([x1, x2, -z1, -z2])
    else:   
        square([y1, y2, -z1, -z2])
    plt.annotate(dens,[(x1 + x2)/2, -(z1+z2)/2])
    
    #%%
    # filter ridges using a minimum length criterium and and filter for a specific range of altitude
    a =2.25
    if x_axis=='y':
        xf_min = a*x1
        xf_max = a*x2
    else:   
        xf_min = -5800
        xf_max = a*x2        

    D_f = dEXP.filter_ridges(dfI,dfII,dfIII,
                            minDepth=200,
                            maxDepth=2000,
                            minlength=8,
                            rmvNaN=True,
                            xmin=xf_min, xmax=xf_max,
                            heights=[hI, hII, hIII])
    
    D_f=D
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
    df_fit = dEXP.fit_ridges(df_f, rmvOutliers=True) # fit ridges on filtered data
    pEXP.plot_ridges_sources(df_fit, ax=ax, z_max_source=-max_elevation*1.2,
                              ridge_type=[0,1,2],ridge_nb=None)
    if x_axis=='y':
        square([x1, x2, -z1, -z2])
    else:   
        square([y1, y2, -z1, -z2])
    plt.annotate(dens,[(x1 + x2)/2, -(z1+z2)/2])

    #%% DEXP ratio

    qratio = [1,0]
    mesh_dexp, label_dexp = dEXP.dEXP_ratio(xp, yp, zp, U, shape, 
                     zmin=0, zmax=max_elevation, nlayers=nlay, 
                     qorders=qratio)
    fig = plt.figure()
    ax = plt.gca()
    
    plt, cmap = pEXP.plot_xy(mesh_dexp, label=label_dexp,
                 markerMax=True,qratio=str(qratio),
                 p1p2=np.array([p1,p2]), ax=ax, Xaxis=x_axis) #, ldg=)
    plt.colorbar(cmap)

    if x_axis=='y':
        square([x1, x2, z1, z2])
        plt.annotate(dens,[(x1 + x2)/2, -(z1+z2)/2])
    else:   
        square([y1, y2, z1, z2])
        plt.annotate(dens,[(y1 + y2)/2, -(z1+z2)/2])
        
    plt.savefig(dataname +  '_fig_DEXP_Ratio_' + x_axis + '.png', r=400)
     
        
    #%% 
    #  ridges analysis: scaling function to determine the SI index

    z0 = -(z1 + z2)/2 + 100 # choose an estimate of the depth of the anomaly
    # z0 = -100 # choose an estimate of the depth of the anomaly
    df_height = D_f[3:6]
    
    ncol = 0
    for r_type in range(len(df_height)): # loop over ridges type I, II, III
        ncol = ncol + df_f[r_type].shape[1]-1



    if ncol<2:
        fig, ax = plt.subplots()
    else:
        fig, axs = plt.subplots(1,ncol, figsize=(15, 6), facecolor='w', edgecolor='k')
        # fig.subplots_adjust(hspace = .5, wspace=.001)
        axs = axs.ravel()

    # for i in range(ncol):
    # df_f[1]['elevation'].iloc[1]
    # dzz = np.log(df_f[1]['elevation'].iloc[1])-np.log(df_f[0]['elevation'].iloc[0])
    # dzz = np.log(df_f[1]['elevation'].iloc[1]-df_f[0]['elevation'].iloc[0])

    nc = 0
    SI_est = []
    for r_type in range(len(df_height)): # loop over ridges type I, II, III
        for k in enumerate(df_height[r_type].columns[1:]): # loop over ridges of the same familly
            # df_f[r_type].columns[1:]    
            points, fit, SI_est_tmp , EXTnb = dEXP.scalFUN(df_height[r_type],EXTnb=k[1],z0=z0,
                                                           rmvOutliers=True)
            if ncol<2:
                pEXP.plot_scalFUN(points, fit, z0=z0, 
                  label='R'+ str(r_type+1) + k[1], 
                  ax=ax) # scaling 
            else:
                pEXP.plot_scalFUN(points, fit, z0=z0, 
                              label='R'+ str(r_type+1) + k[1], 
                              ax=axs[nc]) # scaling 
            SI_est.append(SI_est_tmp)
            
            nc = nc + 1
        
    SI_mean = np.mean(np.abs(SI_est))

    #%% 
    # #  ridges analysis
    # df_f[1]
    # Tau = np.gradient(np.log(up_f_Centralridge)) / np.gradient(np.log(z_r))


    #%% 
    #  ridges analysis
    # SI_mean = 3
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

    plt.savefig(dataname +  '_fig_DEXP_' + x_axis + '.png', r=400)

    #%% 

# uEXP.multipage(dataname + '_test_DEXP_SI.pdf',  dpi=25)
