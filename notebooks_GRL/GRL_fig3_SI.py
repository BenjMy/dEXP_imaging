# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 18:45:42 2020

@author: Benjamin
"""

import pickle

import os

from fatiando.vis import mpl #, myv
from fatiando import gridder, mesher
from fatiando.gravmag import prism, imaging, transform
from fatiando.vis.mpl import square

# my own functions
import lib.dEXP as dEXP
from lib.dEXP import _fit
import lib.plot_dEXP as pEXP
import lib.utils_dEXP as uEXP

# import exemples.fwd_mag_sphere as magfwd
import notebooks_GRL.load_MALM_model as MALMmod
import notebooks_GRL.load_MALM_real as MALMreal

import lib.set_parameters as para

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 15
from mpl_axes_aligner import align

# icsd functions
#from icsd3d.importers.read import load_obs, load_geom


#%% ------------------------------- MALM DATA

interp_size = 300
smooth = 'CubicSmoothingSpline' #'Hanning+Lowpass'
# smooth = 'CubicSmoothingSpline + interp1d' #'Hanning+Lowpass'
interp = False
# path2files = 'E:/Padova/Experiments/GC_2019_Landfill_Porto_MALM/Test/ph_low/'
file = ['NoAno', 'Ano']
path2files = ['./data/phNO/',
              './data/ph_low_top_12.5m/'
              ]

mirror = True
# file = ['Ano']

# path2files = ['E:/Padova/Redaction/Articles/10_InversionUsingGravityPotMethod/data/phNO/',
#               'E:/Padova/Redaction/Articles/10_InversionUsingGravityPotMethod/data/ph_low_top_2m/'
#               'E:/Padova/Redaction/Articles/10_InversionUsingGravityPotMethod/data/ph_low_top_12.5m/'
#               ]

for i, file in enumerate(file):
    
    print(i,file)
    path =  path2files[i] 
    dataset = MALMmod.load_MALM_Landfill_model(path=path, 
                                        filename=file,
                                        shape = (300,300),
                                        field=True,
                                        interp = interp,
                                        radius=60) # length of p1p2 profile

    coord_xyz, coord_xyz_int = dataset[0:2]
    coord_xyz_int
    Uload = dataset[2]
    coords_liner = dataset[3]
    shape, max_elevation = dataset[4:6]

    dict_data = dataset[7]
    dict_data['AnoBool']
    
    xA = (dict_data['HZ'][0][0]+dict_data['HZ'][0][1])/2
    x1 = dict_data['HZ'][0][2]
    y1 = dict_data['HZ'][0][0]
    y2 = dict_data['HZ'][0][1]
    
    z1 = dict_data['HZ'][1]
    z2 = z1 - dict_data['HZ'][2]
    
    p = dataset[6]         # line points                                       
    # set imaging pseudo-inversion parameters                                                                        
    parameters = para.set_par(shape=shape,max_elevation=max_elevation)
    
    scaled = parameters[0]
    SI = parameters[1]
    zp, qorder, nlay = parameters[2:5]
    minAlt_ridge, maxAlt_ridge = parameters[5:7]
    
    max_elevation = 30
    # nlay = 50
    
    # xp, yp, zp = coord_xyz_int
    xp, yp, zp = coord_xyz
    # len(xp)
    Uini = Uload[0] # U_raw, Ucor, U_int, Ucor_int
    p1 , p2 = p

    #%%
    # find point position with respect to line equation defined by p1 and p2
    U_a, p_a, bool_above, U_b, p_b = MALMmod.isabove(xp, yp, Uini, 
                                      np.array(p1),np.array(p2))
    # Compute p1 and p2 line equation ax + by + c = 0
    a, b, c = MALMmod.slope(p1,p2)
    
    # Mirror points with respect to p1p2 line
    Umirror, xy_mirror = MALMmod.mirrorU_alongLine(U_a,p_a,bool_above,a,b,c)
    
    U_a_int = gridder.interp_at(xy_mirror[:,0], xy_mirror[:,1], Umirror, xp, yp, algorithm='cubic', 
                            extrapolate=True)   
    # U_mirror_int = np.copy(U_a_int)
    U_mirror_int = np.copy(Uini)
    U_mirror_int[np.where(bool_above == True)[0]]= U_a_int[np.where(bool_above == True)[0]]

    plt.figure()
    plt.scatter(xp, yp, c=U_mirror_int, cmap='viridis',vmax=0.25)
    plt.colorbar()
    plt.axis('square')
    plt.show()

    
    #%% choose raw or mirrored field

    
    U = np.copy(Uini)
    
    if mirror == True:
        print('mirroring field')
        U = np.copy(U_mirror_int)
        # U = dEXP.smooth2d(xp, yp, U, sigma=2)
    # plt.savefig('smooth2d' + str(file) + '.png', dpi=450)

    # plt.close('all')
    # plt.figure()
    # plt.scatter(xp, yp, c=U, cmap='viridis',vmax=0.25)
    # plt.colorbar()
    # plt.axis('square')
    # plt.show()
    
    #%%
    MainPath= './data/'
    # os.chdir(MainPath)
    ## --------- read MALM measured data file (with electrode positions) --------- ##
    # RealData = np.loadtxt("./1_Data_2_plot/to_plot.dat",skiprows=0,delimiter='\t') 
    out = MALMreal.load_MALM_Porto_real(MainPath,
                              MainPath + 'XYObs_real_f_m3.txt',
                              shape=(200,200),
                              radius=200,
                              rcor=10,
                              rot=60,
                              showfig=False)
    
    # coords_liner = out[3]
    p = out[6]         # line points  
    p1 , p2 = p
    xf, yf, zf = coord_xyz
    # len(xp)
    Uf = Uload[0] # U_raw, Ucor, U_int, Ucor_int

    
    # %% rotate and rescale all
    
    rot = 60
    origin=(max(xp), min(yp))
    point_torotate = np.array([xp, yp])
    xp_r, yp_r = MALMreal.rotate_60(origin, point_torotate, angle_deg=rot, showfig=False)
    Xs = xp_r-min(xp_r)
    Ys = yp_r-min(yp_r) 
    
    
    point_torotate = np.array([[dict_data['HZ'][0][0],dict_data['HZ'][0][1]],
                              [dict_data['HZ'][0][2],dict_data['HZ'][0][2]]])
    xA_r, yA_r =  MALMreal.rotate_60(origin,point_torotate,angle_deg=rot, showfig=False) #58.511
    xA_r = xA_r-min(xp_r)
    yA_r = yA_r-min(yp_r)
    
    # point_torotate = np.array([[p1[0],p2[0]],[p1[1],p2[1]]])
    # px, py = MALM_pm.rotate_60(origin,point_torotate,angle_deg=rot, showfig=False) #58.511
    # px = px-min(xp_r)
    # py = py-min(yp_r)
    
    # px_fix = 284622.86 - min(xp_r)
    # p1_r = [px_fix,py[0]]
    # p2_r = [px_fix,py[1]]
    
    
    # p1_s = [px[0],py[0]]
    # p2_s = [px[0],py[1]]

    point_torotate = np.array(coords_liner).T
    coords_linerx, coords_linery =  MALMreal.rotate_60(origin,point_torotate,angle_deg=rot, showfig=False)
    coords_linerx = coords_linerx-min(xp_r)
    coords_linery = coords_linery-min(yp_r)
    coords_liner_s = np.array([coords_linerx, coords_linery]).T
    
    ax, plt = pEXP.plot_field(Xs,Ys,U, shape,
                    Vminmax=[0,0.35])
    ax.plot(coords_liner_s[2:5,0],coords_liner_s[2:5,1],'k')
    # ax.plot(p1_s,p2_s,'r')
    # # ax.plot(coords_liner[2:5,0],coords_liner[2:5,1],'k')
    # plt.axis('square')
    # # plt.xlim(min(xint_scipy),max(xint_scipy))
    # # plt.ylim(min(yint_scipy),max(yint_scipy))
    plt.xlim(300,500)
    plt.ylim(300,500)

    if i==1:
        plt.savefig('fig2b.pdf', dpi=450)
        plt.savefig('fig2b.svg', dpi=450)
        plt.savefig('fig2b.png', dpi=450)
    else:
        plt.savefig('fig2a.pdf', dpi=450)
        plt.savefig('fig2a.svg', dpi=450)
        plt.savefig('fig2a.png', dpi=450)
        
        
    #%%
    prl = 60
    # shape = shape  (max(xp)-min(xp))/
    shape = (150,150)
    xint_scipy, yint_scipy = gridder.regular((min(Xs)-prl, max(Xs)+prl, 
                              min(Ys)-prl, max(Ys)+prl),shape=shape)
    
    #%% Solution 1
    # extrapolate False and fill with 0 before derivative - mask them later on 
    U_int_scipy = gridder.interp_at(Xs,Ys,U, xint_scipy, yint_scipy, algorithm='nearest', extrapolate=False)
    InterpData = np.array([xint_scipy, yint_scipy, U_int_scipy]).T
    where_are_NaNs = np.isnan(InterpData)
    InterpData[where_are_NaNs] = 0.0074
    xint_scipy, yint_scipy, U_int_scipy = InterpData.T
    
    #%% Solution 2
    # Extrapolate = True
    # U_int_scipy = gridder.interp_at(Xs,Ys, U, xint_scipy, yint_scipy, algorithm='nearest', extrapolate=True)
    
    
    #%%
    Xs, Ys = xint_scipy,yint_scipy
    U = U_int_scipy

    #%%
    # ax, plt = pEXP.plot_field(xp,yp,U, shape,
    #                 Vminmax=[0,0.35])
    # ax.plot(coords_liner[2:5,0],coords_liner[2:5,1],'k')
    # p1_s = np.array([p1[0] -min(xp_r),p2[0] -min(xp_r) ])
    # p2_s = np.array([p1[1] -min(yp_r),p2[1] -min(yp_r) ])    
    
    offset = 255#255 #260
    p1_s = np.array([p1[0] -min(xp_r)+offset,p2[0] -min(xp_r)+offset])
    p2_s = np.array([p1[1] -min(yp_r),p2[1] -min(yp_r) ])
    
    #%%
    
    ax, plt = pEXP.plot_field(Xs,Ys,U, shape,
                    Vminmax=[0,0.35])
    ax.plot(coords_liner_s[2:5,0],coords_liner_s[2:5,1],'k')
    ax.plot(p1_s,p2_s,'r')
    # ax.plot(coords_liner[2:5,0],coords_liner[2:5,1],'k')
    plt.axis('square')
    # plt.xlim(min(xint_scipy),max(xint_scipy))
    # plt.ylim(min(yint_scipy),max(yint_scipy))
    plt.xlim(300,500)
    plt.ylim(300,500)
    


    x_axis = 'x'
    # xx, yy, distance, profile, ax, plt = pEXP.plot_line(Xs, Ys, U ,p1_s,p2_s, 
    #                                                     interp=True, smooth=True, 
    #                                                     xaxis = x_axis)
    p1_s = np.array([p1[0] -min(xp_r)+offset,p1[1] -min(yp_r)])
    p2_s = np.array([p2[0] -min(xp_r)+offset,p2[1] -min(yp_r) +
                      abs(p1[1] -min(yp_r) - (yA_r[0]+yA_r[1])/2)
                      -abs(p2[1] -min(yp_r) - (yA_r[0]+yA_r[1])/2)])
    # p1_s = np.array([p1[0] -min(xp_r)+260,0])
    # p2_s = np.array([p2[0] -min(xp_r)+260,800])
    xx, yy, distance, profile, ax,plt = pEXP.plot_line(Xs, Ys, U ,p1_s,p2_s, 
                                                interp=False,
                                                x_resolution = interp_size,
                                                smooth=smooth, 
                                                xaxis = x_axis,
                                                Vminmax=[0,0.35],
                                                limx=[100,650],
                                                limy=[100,650],
                                                showfig=True)
    # plt.savefig('profile' + str(file) + '.png', dpi=450)

    xA_r_new = [p1_s[0]+xA_r[0]-xA_r[1], p1_s[0]-xA_r[0]+xA_r[1]] 

    # %% SAVE DATA AND PARAMETERS
    
    dict_model_data = { "XYU" : [Xs,Ys,U], 
                "prl" : 60,
                "shape": (150,150),
                "xAyAzA1zA2": [xA_r_new,yA_r,z1,z2],
                "coords_liner": coords_liner_s,
                "p12": [p1_s,p2_s]}
    
    if i==0:
        file2save = './data/NoAno_synth_landflill_data'
    else: 
        file2save = './data/Ano_synth_landflill_data'

    afile = open(file2save + '.pkl', 'wb')
    pickle.dump(dict_model_data, afile)
    afile.close()
    

    # %% ------------------------------- plot publi mirror

    ax, plt = pEXP.plot_field(Xs,Ys,U, shape,Vminmax=[0,0.35])
    ax.plot(coords_liner_s[2:5,0],coords_liner_s[2:5,1],'k')
    plt.axis('square')
    plt.xlim(300,500)
    plt.ylim(300,500)
    # plt.savefig('publi_mirror' + file + '.png', dpi=450)
    
    # %% ------------------------------- Pad the edges of grids
    
    # xp,yp,U, shape = dEXP.pad_edges(xp,yp,U,shape,pad_type=0) # reflexion=5
    # pEXP.plot_line(xp, yp,U,p1,p2, interp=interp)
    
    # %% ------------------------------- Plot the derivatives
    
    xderiv = transform.derivx(Xs, Ys, U, shape,order=0)
    yderiv = transform.derivy(Xs, Ys, U, shape,order=0)
    zderiv = transform.derivz(Xs, Ys, U, shape,order=0)
    
    # interp = True
    pEXP.plot_line(Xs, Ys, xderiv ,p1_s,p2_s,title='xderiv',x_resolution= interp_size,
                    savefig=False, interp=interp, smooth=smooth,  Xaxis=x_axis)

    # plt.savefig('xderiv' + str(file) + '.png', dpi=450)

    pEXP.plot_line(Xs, Ys, yderiv ,p1_s,p2_s,title='yderiv',x_resolution= interp_size,
                    savefig=False, interp=interp, smooth=smooth,  Xaxis=x_axis)

    # plt.savefig('yderiv' + str(file) + '.png', dpi=450)

    pEXP.plot_line(Xs, Ys, zderiv ,p1_s,p2_s,title='zderiv',x_resolution= interp_size,
                    savefig=False, interp=interp, smooth=smooth,  Xaxis=x_axis)
    
    # plt.savefig('zderiv' + str(file) + '.png', dpi=450)

    #%% ------- upward continuation of the field data
    p = [p1_s,p2_s]

    mesh, label_prop = dEXP.upwc(Xs, Ys, zp, U, shape, 
                      zmin=0, zmax=max_elevation, nlayers=nlay, 
                      qorder=qorder)
    
    plt, cmap = pEXP.plot_xy(mesh, label=label_prop,Xaxis=x_axis,
                              Vminmax=[0,0.35], p1p2=p)
    cbar = plt.colorbar(cmap) #,shrink=0.25, pad=0.04)
    cbar.set_label('upwc voltage (V)')
    plt.tight_layout()
    # plt.savefig('upwc voltage' + str(file) + '.png', dpi=450)
    
    #%% DEXP ratio
    # x_axis = 'x'
    qratio = [1,0]
    mesh_dexp, label_dexp = dEXP.dEXP_ratio(Xs, Ys, zp, U, shape, 
                      zmin=0, zmax=max_elevation, nlayers=nlay, 
                      qorders=qratio)
    fig, ax = plt.subplots(figsize=(15,3))
    plt, cmap = pEXP.plot_xy(mesh_dexp, label=label_dexp,
                  markerMax=True,qratio=str(qratio),Vminmax=[0,0.075],
                  p1p2=np.array([p1_s,p2_s]), ax=ax, Xaxis=x_axis) #, ldg=)
    # plt, cmap = pEXP.plot_xy(mesh_dexp, label=label_dexp,
    #              markerMax=True,qratio=str(qratio)
    #              ax=ax, Xaxis=x_axis) #, ldg=)
    cbar = plt.colorbar(cmap,shrink=0.25, pad=0.04)
    cbar.set_label('ratio voltage (V)')
    ax.set_aspect("equal")

    if x_axis=='y':
        square([xA_r_new[0], xA_r_new[1], -z1, -z2])
    else:   
        square([yA_r[0], yA_r[1], -z1, -z2])
    plt.xlim([200,600])
    
    if i==1:
        plt.savefig('fig3c_SI.pdf', dpi=450)
        plt.savefig('fig3c_SI.svg', dpi=450)
        plt.savefig('fig3c_SI.png', dpi=450)
    else:
        plt.savefig('fig3d_SI.pdf', dpi=450)
        plt.savefig('fig3d_SI.svg', dpi=450)
        plt.savefig('fig3d_SI.png', dpi=450)


    # %% ridges identification

    dEXP.ridges_minmax_plot(Xs, Ys, mesh, p1_s, p2_s,
                                          label=label_prop,
                                          interp=interp,x_resolution= interp_size,
                                          smooth=smooth,
                                          fix_peak_nb=2,
                                          method_peak='find_peaks',
                                          showfig=False,
                                          Xaxis=x_axis)  
    #%%
    # or  find_peaks or peakdet or spline_roots
    # dfI,dfII, dfIII = dEXP.ridges_minmax(Xs, Ys, mesh, p1_s, p2_s,interp=interp,x_resolution= interp_size,
    #                                       label=label_prop,fix_peak_nb=2,
    #                                       smooth=smooth, # true = low pass, otherwise specify the filter to apply
    #                                       method_peak='find_peaks',
    #                                       showfig=True,
    #                                       Xaxis=x_axis) 

    D = dEXP.ridges_minmax(Xs, Ys, mesh, p1_s, p2_s,
                                          label=label_prop,
                                          method_peak='find_peaks',
                                          fix_peak_nb=2,
                                          returnAmp=True,
                                          showfig=True,
                                          Xaxis=x_axis,
                                          interp=interp,
                                          x_resolution= interp_size,
                                          smooth = smooth,
                                          qorder=qorder)  
    
    dfI, dfII, dfIII =  D[0:3]
    hI, hII, hIII  = D[3:6]
    heights  = D[3:6]

    #%% ------------------------------- plot ridges over continuated section
        
    # fig = plt.figure()
    # ax = plt.gca()
    # plt, cmap = pEXP.plot_xy(mesh, label=label_prop, ax=ax, Xaxis=x_axis,
    #             Vminmax=[0,0.35], p1p2=p)
    # cbar = plt.colorbar(cmap,shrink=0.25, pad=0.04)
    # cbar.set_label('upwc voltage (V)')
    # plt.tight_layout()
    # pEXP.plot_ridges_harmonic(dfI,dfII,dfIII,ax=ax)
    # plt.xlim([200,600])
    # plt.savefig('ridges_raw_' + str(file) + '.png', dpi=450)

    #%% ------------------------------- filter ridges regionally constrainsted)

    
    dfI_f,dfII_f, dfIII_f = dEXP.filter_ridges(dfI,dfII,dfIII,
                                                minDepth=minAlt_ridge, 
                                                maxDepth=maxAlt_ridge,
                                                minlength=7,rmvNaN=True,
                                                xmin=250, xmax=500,
                                                Xaxis=x_axis)
    
    df_f = dfI_f, dfII_f, dfIII_f
    
    #%% ------------------------------- plot ridges fitted over continuated section
    
    fig, ax1 = plt.subplots(figsize=(15,3))

    plt, cmap = pEXP.plot_xy(mesh, label=label_prop, ax=ax1, Xaxis=x_axis,
              Vminmax=[0,0.35], p1p2=p)
    pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax1,label=False)
    
    df_fit = dEXP.fit_ridges(df_f, rmvOutliers=True) # fit ridges on filtered data
    
    ax2 = pEXP.plot_ridges_sources(df_fit, ax=ax1, z_max_source=-max_elevation,
                          ridge_type=[0,1,2],ridge_nb=None, xmin=250,xmax=500)
    
    cbar = plt.colorbar(cmap) #,shrink=0.25, pad=0.04)
    cbar.set_label('upwc voltage (V)')
    ax1.set_xlabel('y (m)')
    if i==1:
        if x_axis=='y':
            square([xA_r_new[0], xA_r_new[1], z1, z2])
        else:   
            square([yA_r[0], yA_r[1], z1, z2])
        
        fig.savefig('fig3a_SI.pdf', dpi=450)
        fig.savefig('fig3a_SI.svg', dpi=450)
        fig.savefig('fig3a_SI.png', dpi=450)
        
    else:
        
        fig.savefig('fig3b_SI.pdf', dpi=450)
        fig.savefig('fig3b_SI.svg', dpi=450)
        fig.savefig('fig3b_SI.png', dpi=450)
        
