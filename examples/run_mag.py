import os

from fatiando.vis import mpl #, myv
from fatiando import gridder, mesher, utils
from fatiando.gravmag import prism, imaging, transform
from fatiando.vis.mpl import square

# my own functions
import dEXP as dEXP
from dEXP import _fit
import plot_dEXP as pEXP
import set_parameters as para

# exemples
import examples.fwd_mag_sphere as magfwd

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 15


#%% ------------------------------- MAG DATA
# -------------------------------  Model
xp, yp, zp, U, shape, p1, p2, coord= magfwd.load_mag_synthetic()
max_elevation=2*max(coord[:,2])
scaled, SI, zp, qorder, nlay, minAlt_ridge, maxAlt_ridge = para.set_par(shape=shape,max_elevation=max_elevation)
interp = True

#%% ------------------------------- Plot the data 
pEXP.plot_line(xp, yp, U,p1,p2, interp=interp)

#%% ------- upward continuation of the field data

mesh, label_prop = dEXP.upwc(xp, yp, zp, U, shape, 
                 zmin=0, zmax=max_elevation, nlayers=nlay, 
                 qorder=qorder)

plt, cmap = pEXP.plot_xy(mesh, label=label_prop)
plt.colorbar(cmap)
        

# %% ridges identification

dEXP.ridges_minmax_plot(xp, yp, mesh, p1, p2,
                                      label=label_prop,
                                      method_peak='find_peaks')  

# or  find_peaks or peakdet or spline_roots
dfI,dfII, dfIII = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,
                                      label=label_prop,
                                      method_peak='find_peaks')  

 
#%% ------------------------------- plot ridges over continuated section
    
fig = plt.figure()
ax = plt.gca()
pEXP.plot_xy(mesh, label=label_prop, ax=ax) #, ldg=)
pEXP.plot_ridges_harmonic(dfI,dfII,dfIII,ax=ax)

#%% ------------------------------- filter ridges regionally constrainsted)
   

dfI_f,dfII_f, dfIII_f = dEXP.filter_ridges(dfI,dfII,dfIII,
                                            1,maxAlt_ridge,
                                            minlength=8,rmvNaN=True)
df_f = dfI_f, dfII_f, dfIII_f

#%% ------------------------------- plot ridges fitted over continuated section

fig = plt.figure()
ax = plt.gca()

pEXP.plot_xy(mesh, label=label_prop, ax=ax) #, ldg=)
pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax,label=True)

df_fit = dEXP.fit_ridges(df_f) # fit ridges on filtered data

pEXP.plot_ridges_sources(df_fit, ax=ax, z_max_source=-max_elevation*2,
                          ridge_type=[0,1,2],ridge_nb=None)

