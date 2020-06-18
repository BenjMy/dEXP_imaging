import os

from fatiando.vis import mpl #, myv
from fatiando import gridder, mesher, utils
from fatiando.gravmag import prism, imaging, transform
from fatiando.vis.mpl import square

# my own functions
import dEXP as dEXP
from dEXP import _fit
import plot_dEXP as pEXP
import utils_dEXP as uEXP

# exemples
# import exemples.fwd_mag_sphere as magfwd
import examples_in_prep.load_MALM_model as MALM

import set_parameters as para

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 15

# icsd functions
from icsd3d.importers.read import load_obs, load_geom


#%% ------------------------------- MALM DATA synth sensibility

xp, yp, zp, U, max_elevation, shape = MALM.load_MALM_synthetic(ZZ=-3.75,shape=(100,100), field=False)
parameters = para.set_par(shape=shape,max_elevation=max_elevation)
p1, p2 = [min(xp),0], [max(xp), 0]
interp = True
scaled = parameters[0]
SI = parameters[1]
zp, qorder, nlay = parameters[2:5]
minAlt_ridge, maxAlt_ridge = parameters[5:7]

# ----------- ridges analysis
nlay = 25
max_elevation = 20
minAlt_ridge = max_elevation*0.05
maxAlt_ridge = max_elevation*0.65

interp = True
smooth = True 

#%% ------------------------------- MALM DATA
# path2files="example_2add_later/Landfill_3d/Ano_0_EA/"  # Real A position without anomaly
# path2files="example_2add_later/Landfill_3d/Ano_1_BH_EA/"  # Real A position with big hole anomaly
# path2files="example_2add_later/Landfill_3d/Ano_1_SH_EA/"  # Real A position  with small hole anomaly
# Main = 'E:/Padova/Software/SourceInversion/icsd_dev/example_2add_later/Landfill_3d/Ano_1_BH_EA/'
# file = 'OAno_synt'

# Main = 'E:/Padova/Experiments/GC_2019_Landfill_Porto_MALM/Test/ph/'
# file = 'Ano'
# interp = True
# smooth = False 

# dataset = MALM.load_MALM_LandfillPorto(path=Main, 
#                                         filename=file,
#                                         shape = (300,300),
#                                         field=True,
#                                         interp = interp,
#                                         radius=300)
# coord_xyz, coord_xyz_int = dataset[0:2]
# coord_xyz_int
# Uload = dataset[2]
# coords_liner = dataset[3]
# shape, max_elevation = dataset[4:6]



# p = dataset[6]         # line points                                       
# # set imaging pseudo-inversion parameters                                                                        
# parameters = para.set_par(shape=shape,max_elevation=max_elevation)

# scaled = parameters[0]
# SI = parameters[1]
# zp, qorder, nlay = parameters[2:5]
# minAlt_ridge, maxAlt_ridge = parameters[5:7]

# max_elevation = 50
# # nlay = 50

# # xp, yp, zp = coord_xyz_int
# xp, yp, zp = coord_xyz
# # len(xp)
# Uini = Uload[0] # U_raw, Ucor, U_int, Ucor_int
# p1 , p2 = p

# # len(xp)
# # MALM.definep1p2(path=Main, radius=130)
# # MALM.squaremat(r=130)

# #%%
# # find point position with respect to line equation defined by p1 and p2
# U_a, p_a, bool_above, U_b, p_b = MALM.isabove(xp, yp, Uini, 
#                                  np.array(p1),np.array(p2))
# # Compute p1 and p2 line equation ax + by + c = 0
# a, b, c = MALM.slope(p1,p2)

# # Mirror points with respect to p1p2 line
# Umirror, xy_mirror = MALM.mirrorU_alongLine(U_a,p_a,bool_above,a,b,c)

# U_a_int = gridder.interp_at(xy_mirror[:,0], xy_mirror[:,1], Umirror, xp, yp, algorithm='nearest', 
#                         extrapolate=True)   

# plt.figure()
# plt.scatter(xp,yp,c=U_a_int, cmap='viridis',vmax=0.25)
# plt.colorbar()
# plt.axis('square')
   
   
# U = np.copy(Uini)
# U[np.where(bool_above == True)[0]]= U_a_int[np.where(bool_above == True)[0]]
# plt.figure()
# plt.scatter(xp, yp, c=U, cmap='viridis',vmax=0.25)
# plt.colorbar()
# plt.axis('square')
# plt.show()

#%%

# def mirrorImage( a, b, c, x1, y1): 
#  	temp = -2 * (a * x1 + b * y1 + c) /(a * a + b * b) 
#  	x = temp * a + x1 
#  	y = temp * b + y1 
#  	return (x, y) 

# plt.figure()
 
# # Umirror = np.copy(U).tolist()
#   # pmirror = np.copy(p).tolist()
# # pmirror = np.ones(p.shape).tolist()
# Umirror= []
# pmirror= []
# for i, bool_pi in enumerate(zip(c_above,p_a)):
#     # if bool_pi[0] == False:
#     xmir, ymir = mirrorImage(a, b, c, bool_pi[1][0], bool_pi[1][1]); 
#     # plt.scatter(xmir, ymir, c=U_a[i], cmap='viridis')
#     # plt.annotate(str(i)+ '_m', [xmir, ymir])
#     # plt.scatter(bool_pi[1][0],  bool_pi[1][1], c='black', cmap='viridis')
#     # plt.annotate(str(i), [bool_pi[1][0],  bool_pi[1][1]])
#     # plt.annotate(str(i)+ '_m', [xmir, ymir])
#     Umirror.append(U_a[i])
#     pmirror.append([xmir, ymir])
    
    
# p12x=[p1[0],p2[0]]
# p12y=[p1[1],p2[1]]

# plt.plot(p12x,p12y)
# # plt.axis('square')


# pmirror = np.vstack(pmirror)
# Umirror = np.array(Umirror)
# plt.scatter(pmirror[:,0],pmirror[:,1],c=Umirror, cmap='viridis',vmax=0.5)
# plt.axis('square')
    
   
# U = np.copy(Uini)
# U[np.where(bool_above == True)[0]]= Ua[np.where(bool_above == True)[0]]
# plt.figure()
# plt.scatter(xp, yp, c=U, cmap='viridis',vmax=0.25)
# plt.colorbar()
# plt.axis('square')
# plt.show()


# %% change p1p2 axis 

# p1,p2 = uEXP.perp_p1p2(p1,p2, offset=0)

# U = MALM.zeros_sym(Uini,c_above,value=0)
# shape = (300,300)

#%% select part of the grid
# uEXP.mirror(p[:,0],p[:,1],data,a,b)
# uEXP.mirror(xp,yp,data,p1,p2)

#%% ------------------------------- smooth the data 

# U = dEXP.smooth2d(xp, yp, U, sigma=0.5)

#%% ------------------------------- Plot the data 

# _, p1, p2, _ = MALM.definep1p2(path=Main,radius=300)
xx, yy, distance, profile = pEXP.plot_line(xp, yp, U ,p1,p2, interp=interp)
# xx, yy, distance, profile = pEXP.plot_line(xp, yp, U_f ,p1,p2, interp=interp)

#%% ------------------------------- Pad the edges of grids

# xp,yp,U, shape = dEXP.pad_edges(xp,yp,U,shape,pad_type=0) # reflexion=5
# pEXP.plot_line(xp, yp,U,p1,p2, interp=interp)

#%% ------------------------------- Plot the derivatives

xderiv = transform.derivx(xp, yp, U, shape,order=1)
yderiv = transform.derivy(xp, yp, U, shape,order=1)
zderiv = transform.derivz(xp, yp, U, shape,order=1)

# interp = True
pEXP.plot_line(xp, yp, xderiv ,p1,p2,title='xderiv',savefig=False, interp=interp)
pEXP.plot_line(xp, yp, yderiv ,p1,p2,title='yderiv',savefig=False, interp=interp)
pEXP.plot_line(xp, yp, zderiv ,p1,p2,title='zderiv',savefig=False, interp=interp)

#%% ------- upward continuation of the field data

mesh, label_prop = dEXP.upwc(xp, yp, zp, U, shape, 
                  zmin=0, zmax=max_elevation, nlayers=nlay, 
                  qorder=qorder)

plt, cmap = pEXP.plot_xy(mesh, label=label_prop)
plt.colorbar(cmap)
        

#%% ------------------------------- ridges identification
# import Find_peaks_Du_et_al_2006 as fpeak
# from Find_peaks_Du_et_al_2006 import _boolrelextrema, _identify_ridge_lines, _filter_ridge_lines

# # exemple
# xs = np.arange(0, 6, 0.05)
# data = np.sin(xs)
# peakind = fpeak.find_peaks_cwt(data, np.arange(1,10))

# # peakind, xs[peakind],data[peakind]
# # # len(xs)q
# plt.plot(xs,data)
# plt.scatter(xs[peakind],data[peakind])
   

# %% ridges identification

dEXP.ridges_minmax_plot(xp, yp, mesh, p1, p2,
                                      label=label_prop,
                                      fix_peak_nb=4,
                                      interp=interp,smooth=smooth,
                                      method_peak='find_peaks')  

# or  find_peaks or peakdet or spline_roots
dfI,dfII, dfIII = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,interp=interp,
                                      label=label_prop,fix_peak_nb=4,
                                      smooth=smooth,
                                      method_peak='find_peaks')  

# dfI, dfII, dfIII = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,
#                                       label=label_prop,
#                                       # minAlt_ridge=minAlt_ridge,
#                                       maxAlt_ridge=maxAlt_ridge,
#                                       fix_peak_nb=None) 
 
#%% ------------------------------- plot ridges over continuated section
    
fig = plt.figure()
ax = plt.gca()
pEXP.plot_xy(mesh, label=label_prop, ax=ax) #, ldg=)
pEXP.plot_ridges_harmonic(dfI,dfII,dfIII,ax=ax)

#%% ------------------------------- filter ridges regionally constrainsted)
   

# dfI_f,dfII_f, dfIII_f = dEXP.filter_ridges(dfI,dfII,dfIII,
#                                            minAlt_ridge,maxAlt_ridge,
#                                            minlength=5,rmvNaN=True)

dfI_f,dfII_f, dfIII_f = dEXP.filter_ridges(dfI,dfII,dfIII,
                                            1,maxAlt_ridge,
                                            minlength=8,rmvNaN=True)

# dfI_f,dfII_f, dfIII_f = dEXP.filter_ridges(dfI,dfII,dfIII,
#                                            minAlt_ridge,maxAlt_ridge,
#                                            minlength=5,rmvNaN=True,
#                                            xmin=284200)

df_f = dfI_f, dfII_f, dfIII_f

#%% ------------------------------- plot ridges fitted over continuated section

fig = plt.figure()
ax = plt.gca()

pEXP.plot_xy(mesh, label=label_prop, ax=ax) #, ldg=)
pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax,label=True)

df_fit = dEXP.fit_ridges(df_f, rmvOutliers=True) # fit ridges on filtered data

pEXP.plot_ridges_sources(df_fit, ax=ax, z_max_source=-max_elevation*2,
                          ridge_type=[0,1,2],ridge_nb=None)

# pEXP.plot_ridges_sources(df_fit, ax=ax, z_max_source=-max_elevation*2,
#                           ridge_type=[1],ridge_nb=None)
# pEXP.plot_ridges_sources(df_fit, ax=ax, z_max_source=-30, 
#                           ridge_type=[0,1], ridge_nb = [[0,1]])

#%% ------------------------------- save intersection
# TO IMPLEMENT !
# dEXP.ridges_intersection_Z0(df_fit,ax=ax,ridge_nb=[3,4]) # 
# f= fit
# idx = np.argwhere(np.diff(np.sign(f - g))).flatten()
# plt.plot(x[idx], f[idx], 'ro')
# import itertools  
# l = ['Geeks', 'for', 'Geeks']  
# import itertools
# # defining iterator  
# iterators = itertools(l)  
 # # for in loop  
# for i in range(6):  
        
#     # Using next function  
#     print(next(iterators), end = " ")  

# np.diff(df[0]['EX_xpos1'])
# np.diff(df[0]['EX_xpos2']).any()
# print(np.diff(df[i][k[1]]))

#%% ------------------------------- ridges analysis

fig = plt.figure()
ax = plt.gca()

pEXP.plot_xy(mesh, label=label_prop, ax=ax) #, ldg=)
pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax,label=True)

df_fit = dEXP.fit_ridges(df_f) # fit ridges on filtered data

pEXP.plot_ridges_sources(df_fit, ax=ax, z_max_source=-max_elevation*2,
                          ridge_type=[0],ridge_nb=[[0,1,2]])


# --------------------
z0 = -6000
points, fit, SI, EXTnb = dEXP.scalFUN(dfI_f,EXTnb=[3],z0=z0)
pEXP.plot_scalFUN(points, fit, z0=z0)

points.shape
max(points[0][:,0])


fit = np.array([fit])
fit[0][:,1]


# points, fit = dEXP.scalFUN(dfI_f,EXTnb=[1],z0=-10000)
z0est = [1000,2000,3000]
P = []
F = []
SI = []

# Test scalfun for different source depth estimates
for zi in z0est:
    points, fit, si = dEXP.scalFUN(dfI_f,EXTnb=[1],z0=zi)
    P.append(points)
    F.append(fit)
    SI.append(si)
    
fig = plt.figure()
ax = plt.gca()
ax1 = plt.subplot(3,1,1)
pEXP.plot_scalFUN(P, F, ax=ax1, z0=z0est)



# z0est = [0,30,60]
# P = []
# F = []
# for zi in z0est:
#     points, fit = dEXP.scalEuler(dfII,EXTnb=[1],z0=zi)
#     P.append(points)
#     F.append(fit)
    



#%% ------- dEXP continuation of the field data
mesh, label_prop = dEXP.dEXP(xp, yp, zp, U, shape, 
                          zmin=0, zmax=max_elevation, nlayers=nlay, 
                          qorder=qorder, SI=1)
# pEXP.plot_z(mesh)
pEXP.plot_xy(mesh, label=label_prop,markerMax=True)      



#%% sum up plot

# fig = plt.figure()
# ax = plt.gca()

# for orderi in range(3):
#     i = orderi +1
#     axupwc = plt.subplot(2, 3, i)
#     mesh, label_prop = dEXP.upwc(xp, yp, zp, U, shape, 
#                       zmin=0, zmax=max_elevation, nlayers=nlay, 
#                       qorder=orderi) 
#     plt, cmap = pEXP.plot_xy(mesh, label=label_prop, ax=axupwc)
#     plt.colorbar(cmap)
#     axdexp = plt.subplot(2, 3, i+3)
#     mesh, label_prop = dEXP.dEXP(xp, yp, zp, U, shape, 
#                                   zmin=0, zmax=max_elevation, nlayers=nlay, 
#                                   qorder=orderi, SI=1)
#     pEXP.plot_xy(mesh, label=label_prop,markerMax=True, ax=axdexp) 
#     plt.colorbar(cmap)
    