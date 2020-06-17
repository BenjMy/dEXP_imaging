# -*- coding: utf-8 -*-
"""
Imaging methods for potential fields.

Implements the DEXP method described in Fedi and Pilkington (2012). Application on a anomaly of density  (gravimetry).

.. note::

    This is part of a larger project aiming at inverting current sources density (see more at: https://icsd-dev.readthedocs.io/en/latest/)


**References**

Uieda, L., V. C. Oliveira Jr, and V. C. F. Barbosa (2013), Modeling the Earth with Fatiando a Terra, Proceedings of the 12th Python in Science Conference, pp. 91 - 98.

Uieda, L. (2018). Verde: Processing and gridding spatial data using Green's functions. Journal of Open Source Software, 3(29), 957. doi:10.21105/joss.00957

Fedi, M., and M. Pilkington (2012), Understanding imaging methods for potential
field data, Geophysics, 77(1), G13, doi:10.1190/geo2011-0078.1

----
"""

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
import examples_in_prep.load_grav_model as grav

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 15


#%% ------------------------------- GRAV DATA
# -------------------------------  Model
# xp, yp, zp, U, shape,model = gravfwd.fwd_grav_fatiando()

data_struct = grav.load_grav_fatiando(name='grav_models/za_1000zb_1500dens_1200')
# ga, gza = grav.load_grav_pygimli_cylinder()

xp,yp,zp,U = data_struct['xyzg']
shape = data_struct['shape']
model = data_struct['model']
# scaled, SI, zp, qorder, nlay, minAlt_ridge, maxAlt_ridge = para.set_par(shape=shape,max_elevation=max_elevation)


x1, x2, y1, y2, z1, z2 = np.array(model[0].get_bounds())

p1 =[min(yp),0]
p2 =[max(yp),0]

max_elevation=8000
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
                                      fix_peak_nb=2,
                                      method_peak='find_peaks')  

# or  find_peaks or peakdet or spline_roots
dfI,dfII, dfIII = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,
                                      label=label_prop,
                                      fix_peak_nb=2,
                                      method_peak='find_peaks')  

 
#%% ------------------------------- plot ridges over continuated section
    
fig = plt.figure()
ax = plt.gca()
pEXP.plot_xy(mesh, label=label_prop, ax=ax) #, ldg=)
pEXP.plot_ridges_harmonic(dfI,dfII,dfIII,ax=ax)

#%% ------------------------------- filter ridges regionally constrainsted)
   

dfI_f,dfII_f, dfIII_f = dEXP.filter_ridges(dfI,dfII,dfIII,
                                            minAlt_ridge,maxAlt_ridge,
                                            minlength=8,rmvNaN=True)
df_f = dfI_f, dfII_f, dfIII_f

#%% ------------------------------- plot ridges fitted over continuated section

fig = plt.figure()
ax = plt.gca()

pEXP.plot_xy(mesh, label=label_prop, ax=ax) #, ldg=)
pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax,label=True)

df_fit = dEXP.fit_ridges(df_f) # fit ridges on filtered data

pEXP.plot_ridges_sources(df_fit, ax=ax, z_max_source=-max_elevation*1.2,
                          ridge_type=[0,1,2],ridge_nb=None)


# #reload object from file
# file2 = open(r'test.pkl', 'rb')
# new_d = pickle.load(file2)
# file2.close()
