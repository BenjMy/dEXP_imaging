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

# import exemples.fwd_mag_sphere as magfwd
import examples_in_prep.load_MALM_model as MALM
import examples_in_prep.load_MALM_PM as MALM_pm

import set_parameters as para

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 15

# icsd functions
from icsd3d.importers.read import load_obs, load_geom


#%% ------------------------------- MALM DATA

# Main = 'E:/Padova/Experiments/GC_2019_Landfill_Porto_MALM/Test/phNO/'

# Potential electrode grid
# ------------------------
# Main = 'E:/Padova/Software/SourceInversion/icsd_dev/'
# os.chdir(Main)
# path2files="example_2add_later/Landfill_3d/Ano_0_EA/"  # Real A position without anomaly
# # path2files="example_2add_later/Landfill_3d/Ano_1_BH_EA/"  # Real A position with big hole anomaly
# # path2files="example_2add_later/Landfill_3d/Ano_1_SH_EA/"  # Real A position  with small hole anomaly
# # Main = 'E:/Padova/Software/SourceInversion/icsd_dev/example_2add_later/Landfill_3d/Ano_1_BH_EA/'
# # file = 'OAno_synt'
# file = 'ONoAno_synt'


# Mesh grid
# ---------
# path2files = 'E:/Padova/Experiments/GC_2019_Landfill_Porto_MALM/Test/ph_low/'
# file = 'Ano'
path2files = 'E:/Padova/Experiments/GC_2019_Landfill_Porto_MALM/Test/phNO/'
file = 'NoAno'
interp = True
smooth = False
       
dataset = MALM.load_MALM_LandfillPorto(path=path2files, 
                                        filename=file,
                                        shape = (300,300),
                                        field=True,
                                        interp = interp,
                                        radius=30) # length of p1p2 profile
coord_xyz, coord_xyz_int = dataset[0:2]
coord_xyz_int
Uload = dataset[2]
coords_liner = dataset[3]
shape, max_elevation = dataset[4:6]

dict_data = dataset[7]

# HZ = [Ha,zA,thickness]

# dstruct = { "SoilR" : 100, "AnoR" : 1e10, 
#            'AnoBool': Ano==1, 
#            "HZ" : [Ha,zA,thickness], 
#            "XYU" : uz0_grid}

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
xp_r, yp_r, zp_r = coord_xyz
# len(xp)
Uini = Uload[0] # U_raw, Ucor, U_int, Ucor_int
p1 , p2 = p




#%%
MainPath= r'E:\Padova\Software\SourceInversion\Potential_field_imaging\dEXP_imaging\examples_in_prep\\'
# os.chdir(MainPath)
## --------- read MALM measured data file (with electrode positions) --------- ##
# RealData = np.loadtxt("./1_Data_2_plot/to_plot.dat",skiprows=0,delimiter='\t') 
out = MALM_pm.load_MALM_Porto_real(MainPath + '/malm_models/',
                          MainPath + './malm_models/XYObs_real_f_m3.txt',
                          shape=(100,100),
                          radius=200,
                          rcor=10,
                          rot=60,
                          showfig=False)

coords_liner = out[3]
p = out[6]         # line points  
p1 , p2 = p


#%%
# find point position with respect to line equation defined by p1 and p2
U_a, p_a, bool_above, U_b, p_b = MALM.isabove(xp_r, yp_r, Uini, 
                                  np.array(p1),np.array(p2))
# Compute p1 and p2 line equation ax + by + c = 0
a, b, c = MALM.slope(p1,p2)

# Mirror points with respect to p1p2 line
Umirror, xy_mirror = MALM.mirrorU_alongLine(U_a,p_a,bool_above,a,b,c)

U_a_int = gridder.interp_at(xy_mirror[:,0], xy_mirror[:,1], Umirror, xp_r, yp_r, algorithm='nearest', 
                        extrapolate=True)   
U = np.copy(Uini)
U[np.where(bool_above == True)[0]]= U_a_int[np.where(bool_above == True)[0]]
plt.figure()
plt.scatter(xp_r, yp_r, c=U, cmap='viridis',vmax=0.25)
plt.colorbar()
plt.axis('square')
plt.show()

plt.figure()
plt.scatter(xp_r,yp_r,c=U_a_int, cmap='viridis',vmax=0.25)
plt.colorbar()
plt.axis('square')

# plt.figure()
# plt.scatter(xp,yp,c=U_a_int, cmap='viridis',vmax=0.25)
# plt.colorbar()
# plt.axis('square')

#%% Prepare grid

# pEXP.plot_field(xp, yp, U_a_int, shape, Vminmax=[1,2])

# prl = 10
# # shape = shape  (max(xp)-min(xp))/
# shape = (115,115)
# xint_scipy, yint_scipy = gridder.regular((min(xp)-prl, max(xp)+prl, 
#                           min(yp)-prl, max(yp)+prl),shape=shape)

# #%% Solution 1
# # extrapolate False and fill with 0 before derivative - mask them later on 
# # U_int_scipy = gridder.interp_at(xp,yp,U, xint_scipy, yint_scipy, algorithm='nearest', extrapolate=False)
# # InterpData = np.array([xint_scipy, yint_scipy, U_int_scipy]).T
# # where_are_NaNs = np.isnan(InterpData)
# # InterpData[where_are_NaNs] = 0.0074
# # xint_scipy, yint_scipy, U_int_scipy = InterpData.T

# #%% Solution 2
# # Extrapolate = True
# # U_int_scipy = gridder.interp_at(xp,yp,U, xint_scipy, yint_scipy, algorithm='nearest', extrapolate=True)
# U_int_scipy = gridder.interp_at(xp,yp,Uini, xint_scipy, yint_scipy, algorithm='nearest', extrapolate=True)


# xp, yp = xint_scipy,yint_scipy
# U = U_int_scipy


# %% rotate and rescale all

rot = 58.511
origin=(max(xp_r), min(yp_r))
    # point_torotate = np.array([X_raw, Y_raw])

# point_torotate = np.array([xp, yp])
# xp_r, yp_r = MALM_pm.rotate_60(origin, point_torotate, angle_deg=rot, showfig=False)
Xs = xp_r-min(xp_r)
Ys = yp_r-min(yp_r)

point_torotate = np.array([[p1[0],p2[0]],[p1[1],p2[1]]])
px, py = MALM_pm.rotate_60(origin,point_torotate,angle_deg=rot, showfig=False) #58.511
# p1_r = [px[0],py[0]]
# p2_r = [px[1],py[1]]
px = px-min(xp_r)
py = py-min(yp_r)

# px_fix = 284622.86 - min(xp_r)
# p1_r = [px_fix,py[0]]
# p2_r = [px_fix,py[1]]

p1_r = [px[0],py[0]]
p2_r = [px[0],py[1]]


# pEXP.plot_field(xp_r,yp_r,U, shape)
# x_axis = 'y'
# xx, yy, distance, profile, _,_ = pEXP.plot_line(Xs, Ys, U ,p1_r,p2_r, 
#                                             interp=False, smooth=False, 
#                                             xaxis = x_axis,
#                                             Vminmax=[0,1],
#                                             limx=[300,450],
#                                             limy=[300,450]
#                                             )

point_torotate = np.array([[dict_data['HZ'][0][0],dict_data['HZ'][0][1]],
                          [dict_data['HZ'][0][2],dict_data['HZ'][0][2]]])
xA_r, yA_r =  MALM_pm.rotate_60(origin,point_torotate,angle_deg=rot, showfig=False) #58.511
xA_r = xA_r-min(xp_r)
yA_r = yA_r-min(yp_r)

# x1_r = dict_data['HZ'][0][2]
# y1_r = dict_data['HZ'][0][0]
# y2_r = dict_data['HZ'][0][1]


# p1,p2 = uEXP.perp_p1p2([px[0],py[0]],[px[1],py[1]], offset=0)

# p1_r = [p1[0],p2[0]]
# p2_r = [p1[0],p2[1]]

# x_axis = 'x'
# xx, yy, distance, profile = pEXP.plot_line(xp_r, yp_r, U ,p1,p2, 
#                                            interp=True, smooth=True, xaxis = x_axis)
# x_axis = 'y'
# xx, yy, distance, profile = pEXP.plot_line(xp_r, yp_r, U ,p1_r,p2_r, 
#                                            interp=False, smooth=True, xaxis = x_axis)

# rotate_and_rescale_all(X_raw,Y_raw,coordE,p1,p2,coords_liner)


# point_torotate = np.array(coords_liner).T
# coords_linerx, coords_linery =  MALM_pm.rotate_60(origin,point_torotate,angle_deg=rot, showfig=False)
# coords_linerx = coords_linerx-min(xp_r)
# coords_linery = coords_linery-min(yp_r)
# coords_liner = np.array([coords_linerx, coords_linery]).T

#%%
ax, plt = pEXP.plot_field(Xs,Ys,U, shape,
                Vminmax=[0,0.15])
coords_liner_s = np.c_[coords_liner[:,0] -min(xp_r),
                       coords_liner[:,1] -min(yp_r) ]
ax.plot(coords_liner_s[2:5,0],coords_liner_s[2:5,1],'k')
p1_s = np.array([p1[0] -min(xp_r),p2[0] -min(xp_r) ])
p2_s = np.array([p1[1] -min(yp_r),p2[1] -min(yp_r) ])
ax.plot(p1_s,p2_s,'r*')
# ax.plot(p1_r,p2_r,'k')


# %% zeros sym
# U  = MALM.zeros_sym(Uini,bool_above,value=0)
# # shape = (300,300)
# plt.figure()
# plt.scatter(xp, yp, c=U_test, cmap='viridis',vmax=0.25)
# plt.colorbar()
# plt.axis('square')
# plt.show()

#%% select part of the grid
# uEXP.mirror(p[:,0],p[:,1],data,a,b)
# uEXP.mirror(xp,yp,data,p1,p2)

#%% ------------------------------- smooth the data 

Us = dEXP.smooth2d(Xs, Ys, U, sigma=10)
ax, plt = pEXP.plot_field(Xs,Ys,Us, shape, Vminmax=[0,0.15])
ax.plot(coords_liner_s[2:5,0],coords_liner_s[2:5,1],'k')
ax.plot(p1_s,p2_s,'r*')
# ax.plot(coords_liner[2:5,0],coords_liner[2:5,1],'k')

# ax.plot(p1_r,p2_r,'k')
# plt.xlim(300,650)
# plt.ylim(300,650)

#%% ------------------------------- Plot the data 
U = np.copy(Us)
# _, p1, p2, _ = MALM.definep1p2(path=Main,radius=300)
# xx, yy, distance, profile = pEXP.plot_line(xp, yp, U ,p1,p2, interp=interp, smooth=False, xaxis = x_axis)
# xx, yy, distance, profile = pEXP.plot_line(xp, yp, U_f ,p1,p2, interp=interp)
x_axis = 'x'
p1_s = np.array([p1[0] -min(xp_r),p1[1] -min(yp_r) ])
p2_s = np.array([p2[0] -min(xp_r),p2[1] -min(yp_r) ])
xx, yy, distance, profile, ax, plt = pEXP.plot_line(Xs, Ys, U ,p1_s,p2_s, interp=True, smooth=True, xaxis = x_axis)
# ax.set_xlim()

#%% ------------------------------- Pad the edges of grids

# xp,yp,U, shape = dEXP.pad_edges(xp,yp,U,shape,pad_type=0) # reflexion=5
# pEXP.plot_line(xp, yp,U,p1,p2, interp=interp)

#%% ------------------------------- Plot the derivatives

xderiv = transform.derivx(Xs, Ys, U, shape,order=0)
yderiv = transform.derivy(Xs, Ys, U, shape,order=0)
zderiv = transform.derivz(Xs, Ys, U, shape,order=0)

# interp = True
pEXP.plot_line(Xs, Ys, xderiv ,p1_s,p2_s,title='xderiv',
               savefig=False, interp=interp, smooth=smooth,  Xaxis=x_axis)
pEXP.plot_line(Xs, Ys, yderiv ,p1_s,p2_s,title='yderiv',
               savefig=False, interp=interp, smooth=smooth,  Xaxis=x_axis)
pEXP.plot_line(Xs, Ys, zderiv ,p1_s,p2_s,title='zderiv',
               savefig=False, interp=interp, smooth=smooth,  Xaxis=x_axis)


#%% ------- upward continuation of the field data

mesh, label_prop = dEXP.upwc(Xs, Ys, zp, U, shape, 
                  zmin=0, zmax=max_elevation, nlayers=nlay, 
                  qorder=qorder)

plt, cmap = pEXP.plot_xy(mesh, label=label_prop,Xaxis=x_axis)
plt.colorbar(cmap)
        

#%% DEXP ratio
# x_axis = 'x'
qratio = [1,0]
mesh_dexp, label_dexp = dEXP.dEXP_ratio(Xs, Ys, zp, U, shape, 
                 zmin=0, zmax=max_elevation, nlayers=nlay, 
                 qorders=qratio)
fig = plt.figure()
ax = plt.gca()
plt, cmap = pEXP.plot_xy(mesh_dexp, scaled=0, label=label_dexp,
              markerMax=True,qratio=str(qratio),
              p1p2=np.array([p1_s,p2_s]), ax=ax, Xaxis=x_axis) #, ldg=)
# plt, cmap = pEXP.plot_xy(mesh_dexp, label=label_dexp,
#              markerMax=True,qratio=str(qratio)
#              ax=ax, Xaxis=x_axis) #, ldg=)
plt.colorbar(cmap)

if x_axis=='y':
    square([xA_r[0], xA_r[1], -z1, -z2])
    # plt.annotate('\Omega',[x1, -(z1+z2)/2])
else:   
    square([yA_r[0], yA_r[1], -z1, -z2])
    # plt.annotate('\Omega',[(y1 + y2)/2, -(z1+z2)/2])
        
# p1p2=np.array([p1_r,p2_r])
# p_xaxis = []
# for i in p1p2[0]:
#     if(i in p1p2[1]):
#         p_xaxis.append(i)
                
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

# dEXP.ridges_minmax_plot(xp, yp, mesh, p1, p2,
#                                       label=label_prop,
#                                       fix_peak_nb=4,
#                                       interp=interp,smooth=smooth,
#                                       method_peak='find_peaks')  

# or  find_peaks or peakdet or spline_roots
dfI,dfII, dfIII = dEXP.ridges_minmax(Xs, Ys, mesh, p1_s, p2_s,interp=interp,
                                      label=label_prop,fix_peak_nb=2,
                                      smooth=smooth,
                                      method_peak='find_peaks',
                                      showfig=True,
                                      Xaxis=x_axis) 

# dfI, dfII, dfIII = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,
#                                       label=label_prop,
#                                       # minAlt_ridge=minAlt_ridge,
#                                       maxAlt_ridge=maxAlt_ridge,
#                                       fix_peak_nb=None) 
 
#%% ------------------------------- plot ridges over continuated section
    
fig = plt.figure()
ax = plt.gca()
pEXP.plot_xy(mesh, label=label_prop, ax=ax, Xaxis=x_axis)
pEXP.plot_ridges_harmonic(dfI,dfII,dfIII,ax=ax)

#%% ------------------------------- filter ridges regionally constrainsted)
   

# dfI_f,dfII_f, dfIII_f = dEXP.filter_ridges(dfI,dfII,dfIII,
#                                            minAlt_ridge,maxAlt_ridge,
#                                            minlength=5,rmvNaN=True)

dfI_f,dfII_f, dfIII_f = dEXP.filter_ridges(dfI,dfII,dfIII,
                                            minAlt_ridge,maxAlt_ridge,
                                            minlength=8,rmvNaN=True,
                                            xmin=150, xmax=450,
                                            Xaxis=x_axis)

# dfI_f,dfII_f, dfIII_f = dEXP.filter_ridges(dfI,dfII,dfIII,
#                                            minAlt_ridge,maxAlt_ridge,
#                                            minlength=5,rmvNaN=True,
#                                            xmin=284200)

df_f = dfI_f, dfII_f, dfIII_f




#%% ------------------------------- plot ridges fitted over continuated section

fig = plt.figure()
ax = plt.gca()

pEXP.plot_xy(mesh, label=label_prop, ax=ax, Xaxis=x_axis)#, ldg=)
pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax,label=True)

df_fit = dEXP.fit_ridges(df_f, rmvOutliers=True) # fit ridges on filtered data

pEXP.plot_ridges_sources(df_fit, ax=ax, z_max_source=-max_elevation*2,
                          ridge_type=[0,1,2],ridge_nb=None)

# x1, x2, z1, z2 = XXZZ[i]
# square([x1, x1+0.01, z1, z2])
square([y1, y2, z1, z2])

# uEXP.multipage(file+ '.pdf')

# plt.annotate(CTm[i],[(x1 + x2)/2, (z1+z2)/2])


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

# fig = plt.figure()
# ax = plt.gca()

# pEXP.plot_xy(mesh, label=label_prop, ax=ax) #, ldg=)
# pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax,label=True)

# df_fit = dEXP.fit_ridges(df_f) # fit ridges on filtered data

# pEXP.plot_ridges_sources(df_fit, ax=ax, z_max_source=-max_elevation*2,
#                           ridge_type=[0],ridge_nb=[[0,1,2]])


# # --------------------
# z0 = -6000
# points, fit, SI, EXTnb = dEXP.scalFUN(dfI_f,EXTnb=[3],z0=z0)
# pEXP.plot_scalFUN(points, fit, z0=z0)

# points.shape
# max(points[0][:,0])


# fit = np.array([fit])
# fit[0][:,1]


# # points, fit = dEXP.scalFUN(dfI_f,EXTnb=[1],z0=-10000)
# z0est = [1000,2000,3000]
# P = []
# F = []
# SI = []

# # Test scalfun for different source depth estimates
# for zi in z0est:
#     points, fit, si = dEXP.scalFUN(dfI_f,EXTnb=[1],z0=zi)
#     P.append(points)
#     F.append(fit)
#     SI.append(si)
    
# fig = plt.figure()
# ax = plt.gca()
# ax1 = plt.subplot(3,1,1)
# pEXP.plot_scalFUN(P, F, ax=ax1, z0=z0est)



# # z0est = [0,30,60]
# # P = []
# # F = []
# # for zi in z0est:
# #     points, fit = dEXP.scalEuler(dfII,EXTnb=[1],z0=zi)
# #     P.append(points)
# #     F.append(fit)
    



# #%% ------- dEXP continuation of the field data
# mesh, label_prop = dEXP.dEXP(xp, yp, zp, U, shape, 
#                           zmin=0, zmax=max_elevation, nlayers=nlay, 
#                           qorder=qorder, SI=1)
# # pEXP.plot_z(mesh)
# pEXP.plot_xy(mesh, label=label_prop,markerMax=True)      



# #%% sum up plot

# # fig = plt.figure()
# # ax = plt.gca()

# # for orderi in range(3):
# #     i = orderi +1
# #     axupwc = plt.subplot(2, 3, i)
# #     mesh, label_prop = dEXP.upwc(xp, yp, zp, U, shape, 
# #                       zmin=0, zmax=max_elevation, nlayers=nlay, 
# #                       qorder=orderi) 
# #     plt, cmap = pEXP.plot_xy(mesh, label=label_prop, ax=axupwc)
# #     plt.colorbar(cmap)
# #     axdexp = plt.subplot(2, 3, i+3)
# #     mesh, label_prop = dEXP.dEXP(xp, yp, zp, U, shape, 
# #                                   zmin=0, zmax=max_elevation, nlayers=nlay, 
# #                                   qorder=orderi, SI=1)
# #     pEXP.plot_xy(mesh, label=label_prop,markerMax=True, ax=axdexp) 
# #     plt.colorbar(cmap)
    