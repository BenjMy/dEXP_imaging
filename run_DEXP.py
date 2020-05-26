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
shape = (50, 50) # data interpolation
SI = 2 # structural index
zp=0  # initial depth (conditionned upward or downward)
qorder = 1 # derivative order of the continuated field
# ---- z-discretisation - Upward continuation parameters parameters
maxdepth=20
nlay = 25
#square([y1, y2, x1, x2])


# ------------------------------  Model parametes
ZZ = -13.75 # depth of the synthetic anomaly
x1, x2, y1, y2, z1, z2 = -5,5,-5,5,ZZ-2.5/2,ZZ+2.5/2
Rsoil = 1000

# ------------------------------- Load data
filename = 'MSoilR' + str(Rsoil) + 'AnoR1Z'+ str(ZZ) +'L5h2.5'
MainPath='E:/Padova/Simulation/MALM_SensitivityAnalysis/Sens3d/' + filename + '/Data/'
#MainPath='E:/Padova/Simulation/MALM_SensitivityAnalysis/Sens3d/' + filename + '/Data/'
os.chdir(MainPath)
x,y,z,gz=np.loadtxt('3d_SENS_FAT.txt', unpack=True)

# ------------------------------- remove B return electrode effects 
B = [-104.16666666666667, -104.16666666666667, 0]
U = dEXP.cor_field_B(x,y,z,gz,B,rho=100)

# U_cor = U
xp,yp,U = gridder.interp(x,y,U,shape)
# xp,yp,gz_cor= gridder.interp(x,y,gz_cor,shape)

#%% ------------------------------- Plot the data
p1, p2 = [min(xp), 0], [max(xp), 0]
pEXP.plot_line(xp, yp, U,p1,p2)
# #plt.title(strname +'ztop' + str(za) +'_zbot'+ str(zb) + '_data', fontsize=20)
# #plt.savefig(pathFig+strname + '_ExempleFig_z' + str(za) + str(zb) + '_data' + '.png')
# #x1, x2, y1, y2, z1, z2 = np.array(model[0].get_bounds())q
# square([y1, y2, x1, x2])

#%% ------------------------------- Pad the edges of grids

# not implemented yet

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

#%% sum up plot

fig = plt.figure()
ax = plt.gca()

for orderi in range(3):
    i = orderi +1
    axupwc = plt.subplot(2, 3, i)
    mesh, label_prop = dEXP.upwc(xp, yp, zp, U, shape, 
                     zmin=0, zmax=maxdepth, nlayers=nlay, 
                     qorder=orderi) 
    plt, cmap = pEXP.plot_xy(mesh, label=label_prop, ax=axupwc)
    plt.colorbar(cmap)
    axdexp = plt.subplot(2, 3, i+3)
    mesh, label_prop = dEXP.dEXP(xp, yp, zp, U, shape, 
                                  zmin=0, zmax=maxdepth, nlayers=nlay, 
                                  qorder=orderi, SI=1)
    pEXP.plot_xy(mesh, label=label_prop,markerMax=True, ax=axdexp) 
    plt.colorbar(cmap)
        
#%% ------------------------------- ridges identification
dfI,dfII, dfIII = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,
                                     label=label_prop)

#%% ------------------------------- plot ridges over continuated section


fig = plt.figure()
ax = plt.gca()
pEXP.plot_xy(mesh, label=label_prop, ax=ax) #, ldg=)
pEXP.plot_ridges_harmonic(dfI,dfII,dfIII,ax=ax)
# pEXP.plot_ridges_source(R_fit,ax=ax)

#%% ------------------------------- filter ridges regionally constrainsted)

minDepth = 5
maxDepth = 15
dfI_f,dfII_f, dfIII_f = dEXP.filter_Ridges(dfI,dfII,dfIII,minDepth,maxDepth)
df = [dfI_f, dfII_f, dfIII_f]


fig = plt.figure()
ax = plt.gca()
pEXP.plot_xy(mesh, label=label_prop, ax=ax) #, ldg=)
pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax)
pEXP.plot_ridges_sources(df,ax=ax,zlim=[-45,20])

#%% ------------------------------- ridges analysis




#%% ------- dEXP continuation of the field data
mesh, label_prop = dEXP.dEXP(xp, yp, zp, U, shape, 
                          zmin=0, zmax=maxdepth, nlayers=nlay, 
                          qorder=qorder, SI=1)
# pEXP.plot_z(mesh)
pEXP.plot_xy(mesh, label=label_prop,markerMax=True)      
        