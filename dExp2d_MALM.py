"""
GravMag applied to MALM: 3D imaging on synthetic MALM data
(simple example)
"""

from fatiando import gridder, mesher, utils
from fatiando.gravmag import prism, imaging
from fatiando.vis import mpl #, myv
from fatiando.vis.mpl import square
import numpy as np
import matplotlib.pyplot as plt
from fatiando.gravmag import transform

import os

import math 

import dEXP as dEXP

# %matplotlib qt

SI = 2
scaled=0
# -------------------------------  Model
ZZ = -3.75
maxdepth=20
Rsoil = 1000
x1, x2, y1, y2, z1, z2 = -5,5,-5,5,ZZ-2.5/2,ZZ+2.5/2
#square([y1, y2, x1, x2])

filename = 'MSoilR' + str(Rsoil) + 'AnoR1Z'+ str(ZZ) +'L5h2.5'
MainPath='E:/Padova/Simulation/MALM_SensitivityAnalysis/Sens3d/' + filename + '/Data/'
#MainPath='E:/Padova/Simulation/MALM_SensitivityAnalysis/Sens3d/' + filename + '/Data/'

#FigPath = '../../'
os.chdir(MainPath)

#if os.path.exists('../' + 'FiguresPotImg/') == False:
#   os.mkdir('../' + 'FiguresPotImg/')
#   
#pathFig='../' + 'FiguresPotImg/'+ strname + '/'
#pathData=pathFig + '/data/' 
#
##pathFig = os.getcwd() + '/' + strname
#if os.path.exists(pathFig) == False:
#   os.mkdir(pathFig)
#
#if os.path.exists(pathData) == False:
#   os.mkdir(pathData)


# ------------------------------- Load data
import exemples.fwd_mag_sphere as magfwd
xp, yp, zp, gz, shape, p1, p2, coord= magfwd.load_mag_synthetic()


# shape = (50, 50)
# x,y,z,gz=np.loadtxt('3d_SENS_FAT.txt', unpack=True)

# # remove B return electrode effects 
# B = [-104.16666666666667, -104.16666666666667, 0]
# gz_cor = dEXP.cor_field_B(x,y,z,gz,B,rho=100)

# gz_cor = gz
# xp,yp,gz= gridder.interp(x,y,gz,shape)
# xp,yp,gz_cor= gridder.interp(x,y,gz_cor,shape)
# zp=0

# # ------------------------------- Export for CWT analysis

# # Extract a profile between points 1 and 2
# #p1, p2 = [min(xp), 0], [max(xp), 0]
# p1, p2 = [0, min(xp)], [0, max(xp)]
xx, yy, distance, profile = gridder.profile(xp, yp, gz, p1, p2, 1000)
# xx, yy, distance, profile_cor = gridder.profile(xp, yp, gz_cor, p1, p2, 1000)


# ------------------------------- Export for CWT analysis

## Extract a profile between points 1 and 2
#p1, p2 = [0.1, 0.1], [12000, 0.1]
#xx, yy, distance, profile = gridder.profile(xp, yp, gz, p1, p2, 1000)

# Plot the profile and the original map data
plt.figure()
plt.subplot(2, 1, 1)
#plt.title(strname + '_ztop' + str(za) +'_zbot'+ str(zb), fontsize=15)
plt.plot(distance, profile, '.k')
# plt.plot(distance, profile_cor, '.r')
plt.xlim(distance.min(), distance.max())
plt.grid()
plt.subplot(2, 1, 2)
plt.title("Original data")
plt.plot(xx, yy, '-k', label='Profile', linewidth=2)
scale = np.abs([gz.min(), gz.max()]).max()
plt.tricontourf(xp, yp, gz, 50, cmap='RdBu_r', vmin=-scale, vmax=scale)
plt.colorbar(orientation='horizontal', aspect=50)
plt.legend(loc='lower right')
plt.tight_layout()
plt.show()
#plt.suptitle(strname + '_ztop' + str(za) +'_zbot'+ str(zb), fontsize=15)
#plt.savefig(pathFig+ strname + '_ExempleFig_prof' + str(za) + str(zb) + '.png')

# ------------------------------- Plot the data
plt.figure()
plt.axis('scaled')
mpl.contourf(yp, xp, gz, shape, 30)
plt.colorbar()
plt.xlabel('East (km)')
plt.ylabel('North (km)')
mpl.m2km()
plt.show()
#plt.title(strname +'ztop' + str(za) +'_zbot'+ str(zb) + '_data', fontsize=20)
#plt.savefig(pathFig+strname + '_ExempleFig_z' + str(za) + str(zb) + '_data' + '.png')

#x1, x2, y1, y2, z1, z2 = np.array(model[0].get_bounds())
square([y1, y2, x1, x2])

# %% ------------------------------- derivatives

xderiv = transform.derivx(xp, yp, gz, shape,order=1)
yderiv = transform.derivy(xp, yp, gz, shape)
zderiv = transform.derivz(xp, yp, gz, shape,order=1)

plt.subplots_adjust()
levels = 30
cbargs = dict(orientation='horizontal')

ax = plt.subplot(1, 3, 1)
plt.axis('scaled')
plt.title('x derivative')
mpl.contourf(yp, xp, xderiv, shape, levels, cmap=plt.cm.RdBu_r)
plt.colorbar(**cbargs).set_label('u/m')
mpl.m2km()
plt.xlabel('x (m)')
plt.ylabel('y (m)')

ax = plt.subplot(1, 3, 2)
plt.axis('scaled')
plt.title('y derivative')
mpl.contourf(yp, xp, yderiv, shape, levels, cmap=plt.cm.RdBu_r)
plt.colorbar(**cbargs).set_label('u/m')
mpl.m2km()
plt.xlabel('x (m)')
plt.ylabel('y (m)')

ax = plt.subplot(1, 3, 3)
plt.axis('scaled')
plt.title('z derivative')
mpl.contourf(yp, xp, zderiv, shape, levels, cmap=plt.cm.RdBu_r)
cb = plt.colorbar(**cbargs).set_label('u/m')
mpl.m2km()
plt.xlabel('x (m)')
plt.ylabel('y (m)')

#%% ------------------------------- upward continuation of the field data

#
from fatiando.mesher import PrismMesh
def _makemesh(x, y, shape, zmin, zmax, nlayers):
    """
    Make a prism mesh bounded by the data.
    """
    ny, nx = shape
    bounds = [x.min(), x.max(), y.min(), y.max(), zmin, zmax]
    mesh = PrismMesh(bounds, (nlayers, ny, nx))
    return mesh

nlayers=25

mesh = _makemesh(xp, yp, shape, zmin=0, zmax=maxdepth, nlayers=nlayers)
# This way, if z is not an array, it is now
zp = zp * np.ones_like(xp)
  
depths = mesh.get_zs()[:-1]
density = []

up_f_sec = []
up_f_d1_sec = []
up_f_d2_sec = []

for depth in depths - zp[0]:
    
    # continued field calculation
    up_f= transform.upcontinue(xp, yp, gz, shape, depth)
    xx, yy, distance, p_up_f = gridder.profile(xp, yp, up_f, p1, p2, 1000)
    up_f_sec = np.concatenate([up_f_sec, p_up_f])

    # 1st vertical derivate of the continued field
    up_f_d1 = transform.derivz(xp, yp, up_f, shape,order=1)
    xx, yy, distance, p_up_f_d1 = gridder.profile(xp, yp, up_f_d1, p1, p2, 1000)
    up_f_d1_sec = np.concatenate([up_f_d1_sec, p_up_f_d1])

 
    # 2nd vertical derivate of the continued field
    up_f_d2 = transform.derivz(xp, yp, up_f, shape,order=2)
    xx, yy, distance, p_up_f_d2 = gridder.profile(xp, yp, up_f_d2, p1, p2, 1000)
    up_f_d2_sec = np.concatenate([up_f_d2_sec, p_up_f_d2])
    

up_f_w_sec = []
up_f_d1_w_sec = []
up_f_d2_w_sec = []
density = []


for depth in (depths - zp[0]):

    # the continued field weigted (=DEXP)
    up_f_w= up_f*depth**(SI/2)
    xx, yy, distance, p_up_f_w = gridder.profile(xp, yp, up_f_w, p1, p2, 1000)
    up_f_w_sec = np.concatenate([up_f_w_sec, p_up_f_w])

    # the derived continued field weigted (=DEXP)
    up_f_d1_w=  up_f_d1*depth**((SI+1)/2)
    xx, yy, distance, p_up_f_d1_w = gridder.profile(xp, yp, up_f_d1_w, p1, p2, 1000)
    up_f_d1_w_sec = np.concatenate([up_f_d1_w_sec, p_up_f_d1_w])


    # the 2nd derived continued field weigted (=DEXP)
    up_f_d2_w= up_f_d2*depth**((SI+2)/2)
    xx, yy, distance, p_up_f_d2_w = gridder.profile(xp, yp, up_f_d2_w, p1, p2, 1000)
    up_f_d2_w_sec = np.concatenate([up_f_d2_w_sec, p_up_f_d2_w])

    density.extend(up_f_w)
    mesh.addprop('density', np.array(density))

up_f_d2_w_sec_lin = up_f_d2_w_sec
len(up_f_d2_w_sec_lin)

up_f_sec = np.reshape(up_f_sec,(nlayers,1000))
up_f_d1_sec = np.reshape(up_f_d1_sec,(nlayers,1000))
up_f_d2_sec = np.reshape(up_f_d2_sec,(nlayers,1000))
up_f_w_sec = np.reshape(up_f_w_sec,(nlayers,1000))
up_f_d1_w_sec = np.reshape(up_f_d1_w_sec,(nlayers,1000))
up_f_d2_w_sec = np.reshape(up_f_d2_w_sec,(nlayers,1000))

up = [up_f_sec, up_f_d1_sec, up_f_d2_sec, up_f_w_sec, up_f_d1_w_sec, up_f_d2_w_sec]


list_indmax = []
uu_rs = []
for uu in enumerate(up):
    field = uu[1]
    uu_rs.append(np.reshape(field,(nlayers,1000)))
    idmax = np.argmax(field)
    indmax = np.unravel_index(idmax,field.shape)
    list_indmax.append(indmax)



X, Y = np.meshgrid(yy, depths)

#%% Plot the profile and the original map data

plt.figure(figsize=[20, 20])
plt.subplot(3, 3, 1)
plt.plot(xx, profile, '.k')
plt.xlim(xx.min(), xx.max())
plt.xlabel('position (m)',size=20)
plt.ylabel('Field u (V)',size=20)
plt.grid()
#
plt.subplot(3, 3, 2)
d1 = transform.derivz(xp, yp, gz, shape,order=1)
xx, yy, distance, p_d1 = gridder.profile(xp, yp, d1, p1, p2, 1000)
plt.plot(xx, p_d1, '.k')
plt.xlim(xx.min(), xx.max())
plt.xlabel('position (m)',size=20)
plt.ylabel( r'$\frac{\partial u}{\partial z}$',size=20)
plt.title( '1st derivative',size=20)
plt.grid()

plt.subplot(3, 3, 3)
d2 = transform.derivz(xp, yp, gz, shape,order=2)
xx, yy, distance, p_d2 = gridder.profile(xp, yp, d2, p1, p2, 1000)
plt.plot(xx, p_d2, '.k')
plt.xlim(xx.min(), xx.max())
plt.xlabel('position (m)',size=20)
plt.ylabel( r'$\frac{\partial^{2} u}{\partial^{2} z}$',size=20)
plt.title( '2nd derivative',size=20)
plt.grid()

plt.subplot(3, 3, 4)
plt.contourf(X, Y, up_f_sec)
plt.colorbar()
plt.xlabel('position (m)',size=20)
plt.ylabel('Field u continuated (altitude)',size=20)
plt.grid()

id_0_L1 = (np.abs(up_f_d1_sec[0,:] - 0)).argmin()
id_0_Lend = (np.abs(up_f_d1_sec[-1,:] - 0)).argmin()
ind_0_L1 = [(0,id_0_L1),(0,id_0_Lend)]
ind_0_Lend = (24,id_0_Lend)

ax = plt.subplot(3, 3, 5)
plt.contourf(X, Y, up_f_d1_sec)
plt.colorbar()
plt.xlabel('position (m)',size=20)
plt.grid()

plt.subplot(3, 3, 6)
plt.contourf(X, Y, up_f_d2_sec)
plt.colorbar()
plt.xlabel('position (m)',size=20)
plt.grid()

plt.subplot(3, 3, 7)
plt.contourf(X, Y, up_f_w_sec)
plt.colorbar()
plt.scatter(X[list_indmax[3]],Y[list_indmax[3]], s=70, c='w', marker='v')
plt.xlabel('position (m)',size=20)
plt.ylabel(r'$\Omega^{(q=0)}$',size=20)
plt.grid()
#square([x1, x2,z1, z2])
plt.gca().invert_yaxis()

plt.subplot(3, 3, 8)
plt.contourf(X, Y, up_f_d1_w_sec)
plt.colorbar()
plt.scatter(X[list_indmax[4]],Y[list_indmax[4]], s=70, c='w', marker='v')
plt.xlabel('position (m)',size=20)
plt.ylabel(r'$\Omega^{(q=1)}$',size=20)
plt.grid()
#square([x1, x2,z1, z2])
plt.gca().invert_yaxis()

plt.subplot(3, 3, 9)
plt.contourf(X, Y, up_f_d2_w_sec)
plt.colorbar()
plt.scatter(X[list_indmax[5]],Y[list_indmax[5]], s=70, c='w', marker='v')
plt.xlabel('position (m)',size=20)
plt.ylabel(r'$\Omega^{(q=2)}$',size=20)
plt.grid()
#square([x1, x2,z1, z2])
#if len(model)>1:
#    square([x1, x2,z1, z2])
plt.gca().invert_yaxis()

plt.close

# %% ridges detection
# ridges — RI, RII, and RIII — are computed by finding the zeros, 
# at different altitudes, of the horizontal and vertical derivatives of the field 
# and of the field itself, respectively.

from scipy.interpolate import UnivariateSpline

nlayers=25

xx, yy, distance, p_d1 = gridder.profile(xp, yp, d1, p1, p2, 1000)
len(up_f)

mesh = _makemesh(xp, yp, shape, zmin=0, zmax=maxdepth, nlayers=nlayers)
# This way, if z is not an array, it is now
zp = zp * np.ones_like(xp)
depths = mesh.get_zs()[:-1]

RI_0 = [] # zeros of the first horizontal derivative of the potential field
RII_0 = [] # zeros of the first vertical derivative of the potential field
RIII_0 = [] # zeros of the potential field

RI_minmax = [] # zeros of the first horizontal derivative of the potential field
RII_minmax = [] # zeros of the first vertical derivative of the potential field
RIII_minmax = [] # zeros of the potential field


data = gz # or d1, gz
sec = up_f_sec # or up_f_sec, up_f_d1_sec

import matplotlib.pyplot as plt
from scipy.signal import find_peaks
# from scipy.optimize import fmin
# from scipy import optimize


    
    
for i, depth in enumerate(depths - zp[0]):
    fd1z = []
    fd1x = []
    fd = []
    # continued field calculation
    up_f= transform.upcontinue(xp, yp, data, shape, depth)
    xx, yy, distance, p_up_f = gridder.profile(xp, yp, up_f, p1, p2, 1000)
    # fd= UnivariateSpline(yy,p_up_f, s=0)   
    fd= UnivariateSpline(yy,p_up_f)   

    Max_peaks, _ = find_peaks(p_up_f)
    Min_peaks, _ = find_peaks(-p_up_f)
    plt.plot(p_up_f)
    plt.plot(yy[Max_peaks], p_up_f[Max_peaks], "x")
    plt.plot(np.zeros_like(x), "--", color="gray")
    plt.plot(Min_peaks, p_up_f[Min_peaks], "v")
    plt.show()
 
    if fd.roots().any():
        RIII_0.append([depth, np.array(fd.roots())])
    else:
        RIII_0.append([depth, []])

    if np.array(Max_peaks).any():
        RIII_minmax.append([depth, yy[Max_peaks]])
    else:
        RIII_minmax.append([depth, []])
    
    # 1st vertical derivate of the continued field
    up_f_d1z = transform.derivz(xp, yp, up_f, shape,order=1)
    xx, yy, distance, p_up_f_d1z = gridder.profile(xp, yp, up_f_d1z, p1, p2, 1000)
    fd1z= UnivariateSpline(yy,p_up_f_d1z, s=0)
    if fd1z.roots().any():
        RII_0.append([depth, np.array(fd1z.roots())])
        RII_minmax.append([depth, yy[Max_peaks]])
#    RII_max.append(list(fd1x.roots()))

    # 1st horizontal derivate of the continued field
    up_f_d1x = transform.derivx(xp, yp, up_f, shape,order=1)
    xx, yy, distance, p_up_f_d1x = gridder.profile(xp, yp, up_f_d1x, p1, p2, 1000)
    fd1x= UnivariateSpline(yy,p_up_f_d1x, s=0)
    if fd1x.roots().any():
        RI_0.append([depth, np.array(fd1x.roots())])
        RI_minmax.append([depth, yy[Max_peaks]])
    else:
        RI_0.append([depth, []])


# import pandas as pd
# ## Create panda dataframe to merge all the ridges
# RIpd= pd.DataFrame(data=RI_0[0][1],    # values
#             index=RI_0[0][0] - zp[0],    # depths column as index
#               columns='RI')  # 1st row as the column names

RI_0= np.array(RI_0)
RII_0= np.array(RII_0)
RIII_0= np.array(RIII_0)

RI_minmax= np.array(RI_0)
RII_minmax= np.array(RII_0)
RIII_minmax= np.array(RIII_minmax)

plt.scatter(RIII_minmax[:,1], RIII_minmax[:,0])


ax = plt.figure()
plt.contourf(X, Y, sec)
plt.colorbar()
plt.xlabel('position (m)',size=20)
for i, nl in enumerate(depths):
    plt.scatter(RII_0[:,:][i][1],nl*np.ones(len(RII_0[:,:][i][1])), color='red', label='RII 0')
    plt.scatter(RI_0[:,:][i][1],nl*np.ones(len(RI_0[:,:][i][1])), color='green', label='RI 0')
    plt.scatter(RIII_0[:,:][i][1],nl*np.ones(len(RIII_0[:,:][i][1])), color='blue', label='RI 0')
    plt.ylim(min(depths),max(depths))
    # plt.plot(RIII_0[:,:],depths - zp[0], color='red', label='RIII 0')
#plt.plot(RII_max[:,:],depths - zp[0], color='red', label='RII_max')
# plt.grid()
plt.axis('equal')
plt.ylim(min(depths),max(depths))

ax = plt.figure()
plt.contourf(X, Y, sec)
plt.colorbar()
plt.xlabel('position (m)',size=20)
for i, nl in enumerate(depths):
    plt.scatter(RI_minmax[:,:][i][1],nl*np.ones(len(RI_minmax[:,:][i][1])), color='red', label='RII 0')
    plt.scatter(RII_minmax[:,:][i][1],nl*np.ones(len(RII_minmax[:,:][i][1])), color='green', label='RI 0')
    plt.scatter(RIII_minmax[i][1],nl*np.ones(len(RIII_minmax[:,:][i][1])), color='blue', label='RI 0')
    plt.ylim(min(depths),max(depths))
    # plt.plot(RIII_0[:,:],depths - zp[0], color='red', label='RIII 0')
#plt.plot(RII_max[:,:],depths - zp[0], color='red', label='RII_max')
# plt.grid()
plt.axis('equal')
plt.ylim(min(depths),max(depths))


# plt.tight_layout()
 
# 
#    # 2nd vertical derivate of the continued field
#    up_f_d2 = transform.derivz(xp, yp, up_f, shape,order=2)
#    xx, yy, distance, p_up_f_d2 = gridder.profile(xp, yp, up_f_d2, p1, p2, 1000)
#    up_f_d2_sec = np.concatenate([up_f_d2_sec, p_up_f_d2])



# %% compute the scaling function

# Extract ridge from up_f_d1_sec (2nd vertical derivate of the continued field)

Xlin=np.reshape(X,[len(up_f_d2_w_sec_lin)])
Ylin=np.reshape(Y,[len(up_f_d2_w_sec_lin)])
D = np.reshape(up_f_sec,[len(up_f_d2_w_sec_lin)])
D1 = np.reshape(up_f_d1_sec,[len(up_f_d2_w_sec_lin)])
D2 = np.reshape(up_f_d2_sec,[len(up_f_d2_w_sec_lin)])

x_r, z_r, d_r, up_f_Centralridge = gridder.profile(Xlin, Ylin, D, [0,6], [0,10.5], 1000)
x_r, z_r, d_r, up_f_d1_Centralridge = gridder.profile(Xlin, Ylin, D1, [0,6], [0,10.5], 1000)
x_r, z_r, d_r, up_f_d2_Centralridge = gridder.profile(Xlin, Ylin, D2, [0,6], [0,10.5], 1000)

plt.figure()
plt.subplot(2, 1, 1)
plt.plot(z_r, up_f_d1_Centralridge, '.k')
plt.xlim(distance.min(), distance.max())
plt.grid()
plt.subplot(2, 1, 2)
plt.title( r'$\frac{\partial u}{\partial z}$',size=20)
plt.plot(x_r, z_r, '-k', label='Profile', linewidth=2)
#scale = np.abs([gz.min(), gz.max()]).max()
plt.tricontourf(Xlin, Ylin, D, 15)
plt.colorbar(orientation='horizontal', aspect=50)
plt.legend(loc='lower right')
plt.tight_layout()
plt.show()

#  Scaling function for the central ridge


#Tau_d1 = np.gradient(np.log(up_f_d1_Centralridge), np.gradient(np.log(z_r)))
#Tau_d1 = np.gradient(np.log(up_f_d1_Centralridge), np.diff(np.log(z_r)))

# Initial z0 guess
z0 = -13
#z_r = z_r - z0
#max(z_r)

plt.figure()
plt.plot(up_f_Centralridge)
plt.plot(up_f_d1_Centralridge)
plt.plot(up_f_d2_Centralridge)
plt.legend(['f',r'$f_{1}=\frac{\partial u}{\partial z}$', r'$f_{2}=\frac{\partial^{2} u}{\partial^{2} z}$'])
plt.title('Extracted ridges')

#Tau = np.log(up_f_Centralridge)/ np.log(z_r)
Tau = np.gradient(np.log(up_f_Centralridge)) / np.gradient(np.log(z_r))
Tau_d1 = np.gradient(np.log(up_f_d1_Centralridge)) / np.gradient(np.log(z_r))
Tau_d2 = np.gradient(np.log(up_f_d2_Centralridge)) / np.gradient(np.log(z_r))

#plt.figure()
#plt.plot(z_r)
q = 1./z_r

factor = (z_r - z0)/z_r
Tau = Tau*factor
Tau_d1 = Tau_d1*factor
Tau_d2 = Tau_d2*factor

#q = q - 1/z0


np.isnan(np.log(up_f_d1_Centralridge)).any()
np.isnan(np.log(z_r)).any()


#Tau_d1 = np.diff(np.log(up_f_d1_Centralridge))/ np.log(z_r)

#Tau_d1 = np.diff(up_f_d1_Centralridge*z0)/ dz
#Tau_d1 = np.gradient(up_f_d1_Centralridge*z0 ,dz)

from scipy.optimize import curve_fit
def f(x, A, B): # this is your 'straight line' y=f(x)
    return A*x + B

popt, pcov = curve_fit(f, q,Tau) # your data x, y to fit
poptd1, pcov = curve_fit(f, q,Tau_d1) # your data x, y to fit
poptd2, pcov = curve_fit(f, q,Tau_d2) # your data x, y to fit

x_min = 0  
x_max = max(q)                                #min/max values for x axis
x_fit = np.linspace(x_min, x_max, 100)   #range of x values used for the fit function

#plt.plot(z_r)

plt.figure(figsize=[15,5])
plt.subplot(1,3,1)
plt.plot(x_fit, f(x_fit, *popt), 'g--',
         label='fit')
plt.scatter(q,Tau,marker='*')
plt.xlim([0,max(q)])
plt.ylim([-5,5])
plt.xlabel('q (m)', size=20)
plt.ylabel('$\\tau_{f}$', size=20)
plt.title(r'$\frac{\partial log(f)}{\partial log(z)}$', size=20)
plt.grid()

plt.subplot(1,3,2)
plt.plot(x_fit, f(x_fit, *poptd1), 'g--',
         label='fit')
plt.scatter(q,Tau_d1,marker='*')
plt.xlim([0,max(q)])
plt.xlabel('q (m)', size=20)
plt.ylim([-5,5])
plt.ylabel('$\\tau_{f_{1}}$', size=20)
plt.title(r'$\frac{\partial log(f_{1})}{\partial log(z)}$', size=20)
plt.grid()

plt.subplot(1,3,3)
plt.plot(x_fit, f(x_fit, *poptd2), 'g--',
         label='fit')
plt.scatter(q,Tau_d2,marker='*')
plt.xlim([0,max(q)])
plt.ylim([-5,5])
plt.xlabel('q (=1/z) (m)', size=20)
plt.ylabel('$\\tau_{f_{2}}$', size=20)
plt.title(r'$\frac{\partial log(f_{2})}{\partial log(z)}$', size=20)
plt.grid()

#font = {'family' : 'normal',
#        'weight' : 'bold',
#        'size'   : 22}
#
#mpl.rc('font', **font)

# %% ------------------------------- Compute ratio of derivative

#order=1
#xl = np.log(xp)
#yl = np.log(yp)
#gzl = np.log(np.abs(gz))
#
##xl = xp
##yl = yp
##gzl = gz
#
#Tn = transform.derivz(xl, yl, gz, shape, order=1, method='fft')
#Tm = transform.derivz(xl, yl, gz, shape, order=0, method='fft')
#ratio = Tm/Tn
##
##tga = transform.tga(x, y, gz, shape)
##tilt = transform.tilt(x, y, gz, shape)
#
#data2plot = ratio
## Plot the profile and the original map data
#plt.figure()
#scale = np.abs([data2plot.min(), data2plot.max()]).max()
#plt.tricontourf(xp, yp, data2plot, 50, cmap='RdBu_r')
#plt.colorbar(orientation='horizontal', aspect=50)
#plt.legend(loc='lower right')
#plt.tight_layout()
#plt.show()


# %%

# plt.plot(x_fit, f(x_fit, *poptd2), 'g--',
#          label='fit')
# plt.scatter(q,Tau_d2,marker='*')
# plt.xlim([0,max(q)])
# plt.ylim([-5,5])
# plt.xlabel('q (=1/z) (m)', size=20)
# plt.ylabel('$\\tau_{f_{2}}$', size=20)
# plt.title(r'$\frac{\partial log(f_{2})}{\partial log(z)}$', size=20)
# plt.grid()

# plt.plot(up_f_d1_sec[0,:],'*')
# https://stackoverflow.com/questions/9148927/matplotlib-extended-line-over-2-control-points

a= 0
b = up_f_d1_sec[0,:]

d=np.sort(abs(b-a));
if d[1] == d[2]:
      vals = np.where(abs(b-a)==d[1])
      lowest = vals(1)
      second_lowest = vals(2)
else:
      lowest=np.where(abs(b-a)==d[1])
      sec_lowest=np.where(abs(b-a)==d[2])

point1= [-20,12]
point2= [X[ind_0_L1[0]],Y[ind_0_L1[0]]]


def slope_from_points(point1, point2):
    return (point2[1] - point1[1])/(point2[0] - point1[0])

def plot_secant(point1, point2, ax):
    # plot the secant
    slope = slope_from_points(point1, point2)
    intercept = point1[1] - slope*point1[0] 
    # update the points to be on the axes limits
    x = ax.get_xlim()
    y = ax.get_ylim()
    data_y = [x[0]*slope+intercept, x[1]*slope+intercept]
    line = plt.plot(x, data_y, color='red')
    # ax.add_line(line)
    return x, data_y

fig = plt.subplots()
ax = plt.gca()
# plt.scatter(X[ind_0_L1[0]],Y[ind_0_L1[0]], s=70, c='w', marker='v')
x, data_y = plot_secant(point1,point2,ax)
# plt.plot(x, data_y, color='red')
x, data_y = plot_secant([20,12],[15,0],ax)
# plt.plot(x, data_y, color='red')
# plt.scatter(X[ind_0_Lend],Y[ind_0_Lend], s=70, c='w', marker='v')
