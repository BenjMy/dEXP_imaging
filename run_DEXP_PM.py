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
import examples_in_prep.load_MALM_PM as MALM

import set_parameters as para

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 15


#%%
MainPath= r'E:\Padova\Software\SourceInversion\Potential_field_imaging\dEXP_imaging\examples_in_prep\\'
os.chdir(MainPath)

## --------- read MALM measured data file (with electrode positions) --------- ##

# RealData = np.loadtxt("./1_Data_2_plot/to_plot.dat",skiprows=0,delimiter='\t') 
out = MALM.load_MALM_Porto_real(MainPath + '/malm_models/',
                          './malm_models/XYObs_real_f_m3.txt',
                          shape=(100,100),
                          radius=100,
                          rcor=100)

coord_xyz, coord_xyz_int = out[0:2]
Uload = out[2]
coords_liner = out[3]
shape, max_elevation = out[4:6]

p = out[6]         # line points                                       
# set imaging pseudo-inversion parameters                                                                        
parameters = para.set_par(shape=shape,max_elevation=max_elevation)

scaled = parameters[0]
SI = parameters[1]
zp, qorder, nlay = parameters[2:5]
minAlt_ridge, maxAlt_ridge = parameters[5:7]

max_elevation = 50
# nlay = 50

xp, yp, zp = coord_xyz_int
# xp, yp, zp = coord_xyz
# len(xp)
U = Uload[3] # U_raw, Ucor, U_int, Ucor_int
p1 , p2 = p

interp=True
smooth = False
# len(U)
# len(xp)

# %% change p1p2 axis 

# p1,p2 = uEXP.perp_p1p2(p1,p2, offset=0)

# %% zeros sym

#%% select part of the grid
# uEXP.mirror(p[:,0],p[:,1],data,a,b)
# uEXP.mirror(xp,yp,data,p1,p2)

#%% ------------------------------- smooth the data 

# U = dEXP.smooth2d(xp, yp, U, sigma=10)

#%% ------------------------------- Plot the data 

# _, p1, p2, _ = MALM.definep1p2(path=Main,radius=300)
xx, yy, distance, profile = pEXP.plot_line(xp, yp, U ,p1,p2, interp=interp)
# xx, yy, distance, profile = pEXP.plot_line(xp, yp, U_f ,p1,p2, interp=interp)

#%% ------------------------------- Pad the edges of grids

# xp,yp,U, shape = dEXP.pad_edges(xp,yp,U,shape,pad_type=0) # reflexion=5
# pEXP.plot_line(xp, yp,U,p1,p2, interp=interp)

#%% ------------------------------- Plot the derivatives

# xderiv = transform.derivx(xp, yp, U, shape,order=1)
# yderiv = transform.derivy(xp, yp, U, shape,order=1)
# zderiv = transform.derivz(xp, yp, U, shape,order=1)

# # interp = True
# pEXP.plot_line(xp, yp, xderiv ,p1,p2,title='xderiv',savefig=False, interp=interp)
# pEXP.plot_line(xp, yp, yderiv ,p1,p2,title='yderiv',savefig=False, interp=interp)
# pEXP.plot_line(xp, yp, zderiv ,p1,p2,title='zderiv',savefig=False, interp=interp)

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
                                      method_peak='find_peaks'
                                      show)  

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
                                            minAlt_ridge,maxAlt_ridge,
                                            minlength=8,rmvNaN=True
                                            )

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

