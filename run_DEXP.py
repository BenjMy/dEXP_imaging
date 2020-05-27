import os

from fatiando.vis import mpl #, myv
from fatiando import gridder, mesher, utils
from fatiando.gravmag import prism, imaging, transform
from fatiando.vis.mpl import square

import dEXP as dEXP
import plot_dEXP as pEXP

import numpy as np
import matplotlib.pyplot as plt

# -------------------------------  Graphical parameters
scaled=0 

# -------------------------------  Imaging parameters
# shape = (50, 50) # data interpolation
SI = 2 # structural index
zp=0  # initial depth (conditionned upward or downward)
qorder = 0 # derivative order of the continuated field
# ---- z-discretisation - Upward continuation parameters parameters
maxdepth=20
nlay = 25
#square([y1, y2, x1, x2])
minDepth = 5
maxDepth = 15

# # ------------------------------  Model parametes
# ZZ = -3.75 # depth of the synthetic anomaly
# x1, x2, y1, y2, z1, z2 = -5,5,-5,5,ZZ-2.5/2,ZZ+2.5/2
# Rsoil = 1000

# # ------------------------------- Load data
# filename = 'MSoilR' + str(Rsoil) + 'AnoR1Z'+ str(ZZ) +'L5h2.5'
# MainPath='E:/Padova/Simulation/MALM_SensitivityAnalysis/Sens3d/' + filename + '/Data/'
# #MainPath='E:/Padova/Simulation/MALM_SensitivityAnalysis/Sens3d/' + filename + '/Data/'
# os.chdir(MainPath)
# x,y,z,gz=np.loadtxt('3d_SENS_FAT.txt', unpack=True)

# # ------------------------------- remove B return electrode effects 
# B = [-104.16666666666667, -104.16666666666667, 0]
# U = dEXP.cor_field_B(x,y,z,gz,B,rho=100)

# # U_cor = U
# xp,yp,U = gridder.interp(x,y,U,shape)
# # xp,yp,gz_cor= gridder.interp(x,y,gz_cor,shape)

#%% ------------------------------- GRAVITY DATA
# -------------------------------  Model
name = '3000_zbot5000_data'
xp, yp, zp, U= np.loadtxt('./grav_models/' + name + '.txt', unpack=True)
shape = (25,25) # data interpolation
maxdepth=10000 # max depth for upward continuation
minAlt_ridge = 1000
maxAlt_ridge = 5000

#%% ------------------------------- Plot the data
p1, p2 = [min(xp), 0], [max(xp), 0]
pEXP.plot_line(xp, yp, U,p1,p2)
# #plt.title(strname +'ztop' + str(za) +'_zbot'+ str(zb) + '_data', fontsize=20)
# #plt.savefig(pathFig+strname + '_ExempleFig_z' + str(za) + str(zb) + '_data' + '.png')
# #x1, x2, y1, y2, z1, z2 = np.array(model[0].get_bounds())q
# square([y1, y2, x1, x2])

#%% ------------------------------- Pad the edges of grids

# padtypes = ['0', 'mean', 'edge', 'lintaper', 'reflection', 'oddreflection',
#             'oddreflectiontaper']
# fig = plt.figure()
# ax = plt.gca()

# xs = xp.reshape(shape)
# ys = yp.reshape(shape)
# data = U.reshape(shape)

# padtype = padtypes[3]
# padded_data, nps = gridder.pad_array(data, padtype=padtype)
# # Get coordinate vectors
# pad_x, pad_y = gridder.pad_coords([xs, ys], shape, nps)
# padshape = padded_data.shape
# ax.set_title(padtype)
# ax.pcolormesh(pad_y.reshape(padshape), pad_x.reshape(padshape),
#               padded_data, cmap='RdBu_r')
# ax.set_xlim(pad_y.min(), pad_y.max())
# ax.set_ylim(pad_x.min(), pad_x.max())
    
# shape = padded_data.shape
# U = padded_data.reshape(shape[0]*shape[1])
# xp = pad_x
# yp = pad_y

# p1, p2 = [min(xp), 0], [max(xp), 0]


#%% ------------------------------- Plot the derivatives
xderiv = transform.derivx(xp, yp, U, shape,order=1)
yderiv = transform.derivy(xp, yp, U, shape,order=1)
zderiv = transform.derivz(xp, yp, U, shape,order=1)

ax1 = pEXP.plot_line(xp, yp, zderiv ,p1,p2,title='zderiv',savefig=True)

#%% ------- upward continuation of the field data

mesh, label_prop = dEXP.upwc(xp, yp, zp, U, shape, 
                 zmin=0, zmax=maxdepth, nlayers=nlay, 
                 qorder=qorder)

plt, cmap = pEXP.plot_xy(mesh, label=label_prop)
plt.colorbar(cmap)
        
#%% ------------------------------- ridges identification
dfI,dfII, dfIII = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,
                                      label=label_prop) 
dfI.head(5)
dfII.head(5)
dfIII.head(5)

#%% ------------------------------- plot ridges over continuated section

fig = plt.figure()
ax = plt.gca()
pEXP.plot_xy(mesh, label=label_prop, ax=ax) #, ldg=)
pEXP.plot_ridges_harmonic(dfI,dfII,dfIII,ax=ax)
# pEXP.plot_ridges_source(R_fit,ax=ax)

#%% ------------------------------- filter ridges regionally constrainsted)

minAlt_ridge = 500
maxAlt_ridge = 3000
dfI_f,dfII_f, dfIII_f = dEXP.filter_Ridges(dfI,dfII,dfIII,minAlt_ridge,maxAlt_ridge)
df = [dfI_f, dfII_f, dfIII_f]

fig = plt.figure()
ax = plt.gca()
pEXP.plot_xy(mesh, label=label_prop, ax=ax) #, ldg=)
pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax)
points, fit = dEXP.fit_ridges(df) # fit ridges
pEXP.plot_ridges_sources(points, fit, ax=ax, z_max_source=-3000)

# from scipy.optimize import curve_fit

# def f(x, A, B): # this is your 'straight line' y=f(x)
#     return A*x + B

# df = [dfI_f,dfII_f]

# plt.figure()


   
#         plt.plot(df[i][k[1]],df[i]['depth'] , 'b.',
#                   label='points')
#         plt.plot(x_fit,y_fit, 'k--',
#                   label='fit')
#     # plt.scatter(df[i]['EX_xpos1'],df[i]['depth'],marker='*')

# # plt.ylim([zlim[0],zlim[1]])
    

#%% ------------------------------- ridges analysis

# points, fit = dEXP.scalFUN(dfI_f,EXTnb=[1],z0=-10000)
z0est = [0,3,6]
P = []
F = []
for zi in z0est:
    points, fit = dEXP.scalFUN(dfII_f,EXTnb=[1],z0=zi)
    P.append(points)
    F.append(fit)
    
points2, fit = dEXP.scalFUN(dfII_f,EXTnb=[1],z0=3)

fig = plt.figure()
ax = plt.gca()
ax1 = plt.subplot(3,1,1)
pEXP.plot_scalFUN(P, F, ax=ax1, z0=z0est)


# plt.subplot(3,1,2)
# plot_scalFUN(points, fit, ax=ax, z0=z0est)
# plt.subplot(3,1,2)
# plot_scalFUN(points, fit, ax=ax, z0=z0est)


# dfI_f.head(5)
# EXTnb=[1,2]




#%% ------- dEXP continuation of the field data
mesh, label_prop = dEXP.dEXP(xp, yp, zp, U, shape, 
                          zmin=0, zmax=maxdepth, nlayers=nlay, 
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
#                       zmin=0, zmax=maxdepth, nlayers=nlay, 
#                       qorder=orderi) 
#     plt, cmap = pEXP.plot_xy(mesh, label=label_prop, ax=axupwc)
#     plt.colorbar(cmap)
#     axdexp = plt.subplot(2, 3, i+3)
#     mesh, label_prop = dEXP.dEXP(xp, yp, zp, U, shape, 
#                                   zmin=0, zmax=maxdepth, nlayers=nlay, 
#                                   qorder=orderi, SI=1)
#     pEXP.plot_xy(mesh, label=label_prop,markerMax=True, ax=axdexp) 
#     plt.colorbar(cmap)
    