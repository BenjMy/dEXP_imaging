"""
Example of gravimetric potential field data analysis using pyDEXP
-----------------------------------------------------------------

This code shows a step-by-step processing of potential field imaging aiming at giving an estimate of the depth of the anomaly depth using the dEXP tranformation method.
dEXP method implementation from Fedi et al. 2012. 
Calculations used :mod:`dEXP`, while plotting use the :mod:`plot_dEXP` module.

The gravimetric model data was created using geometric objects from :mod:`fatiando.mesher`. The forward simulation of the data was done using :mod:`fatiando.gravmag` module.

Sources locations:
    Center of mass = [,,] # xyz coordinates
    l.w.h = ,,, # length, width, height (in m)
Sources properties: 
    density contrast= ?

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
import examples.gravimetry.loadgrav.grav_models as grav

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 15


#%% 
# load model previously generated using Fatiando a terra package

data_struct = grav.load_grav_fatiando(name='loadgrav/za3000_zb3500_l500_ofs0_dens1200')
xp,yp,zp,U = data_struct['xyzg']
shape = data_struct['shape']
model = data_struct['model']
dens  = data_struct['density']
# scaled, SI, zp, qorder, nlay, minAlt_ridge, maxAlt_ridge = para.set_par(shape=shape,max_elevation=max_elevation)

x1, x2, y1, y2, z1, z2 = np.array(model[0].get_bounds())

p1 =[min(yp),0]
p2 =[max(yp),0]

max_elevation=z2*1.2
scaled, SI, zp, qorder, nlay, minAlt_ridge, maxAlt_ridge = para.set_par(shape=shape,max_elevation=max_elevation)
interp = True
qorder = 0

#%% 
# Plot the data 
pEXP.plot_line(xp, yp, U,p1,p2, interp=interp)

#%% 
# Pad the edges of grids (if necessary)

xp,yp,U, shape = dEXP.pad_edges(xp,yp,U,shape,pad_type=0) # reflexion=5
p1 =[min(yp),0]
p2 =[max(yp),0]

pEXP.plot_line(xp, yp,U,p1,p2, interp=interp)

#%% 
# Upward continuation of the field data

mesh, label_prop = dEXP.upwc(xp, yp, zp, U, shape, 
                 zmin=0, zmax=max_elevation, nlayers=nlay, 
                 qorder=qorder)

# plt, cmap = pEXP.plot_xy(mesh, label=label_prop)
# plt.colorbar(cmap)


#%%
# ridges identification
# dEXP.ridges_minmax_plot(xp, yp, mesh, p1, p2,
#                                       label=label_prop,
#                                       fix_peak_nb=2,
#                                       method_peak='find_peaks')  

# or  find_peaks or peakdet or spline_roots
dfI,dfII, dfIII = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,
                                      label=label_prop,
                                      fix_peak_nb=2,
                                      method_peak='find_peaks')  

#%% 
# Plot ridges over continuated section

fig = plt.figure()
ax = plt.gca()
pEXP.plot_xy(mesh, label=label_prop, ax=ax) #, ldg=)
pEXP.plot_ridges_harmonic(dfI,dfII,dfIII,ax=ax)

#%% 
# Filter ridges regionally constrainsted)
   

dfI_f,dfII_f, dfIII_f = dEXP.filter_ridges(dfI,dfII,dfIII,
                                            minDepth=1000,maxDepth=3000,
                                            minlength=3,rmvNaN=True)
df_f = dfI_f, dfII_f, dfIII_f

#%% 
# Plot ridges fitted over continuated section

fig = plt.figure()
ax = plt.gca()

pEXP.plot_xy(mesh, label=label_prop, ax=ax) #, ldg=)
pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax,label=True)

df_fit = dEXP.fit_ridges(df_f, rmvOutliers=True) # fit ridges on filtered data

pEXP.plot_ridges_sources(df_fit, ax=ax, z_max_source=-max_elevation*1.2,
                          ridge_type=[0,1,2],ridge_nb=None)
square([x1, x2, -z1, -z2])
plt.annotate(dens,[(x1 + x2)/2, -(z1+z2)/2])

