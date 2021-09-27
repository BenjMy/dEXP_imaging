"""
Magnetic field data analysis using pyDEXP: a 2-sources case
-----------------------------------------------------------

This code shows a step-by-step processing of potential field imaging aiming at giving an estimate of magnetic sources positions and depth using the dEXP tranformation method.
dEXP method implementation from Fedi et al. 2012. 
Calculations used :mod:`dEXP`, while plotting use the :mod:`plot_dEXP` module.

The model data was created using geometric objects from :mod:`fatiando.mesher`. The forward simulation of the data was done using :mod:`fatiando.gravmag` module.

Sources locations:
    - S_{A} = [10e3,10e3,2e3] # xyz coordinates
    - S_{B} = [25e3,10e3,1e3]

Sources properties: 
    - radius = 1.5e3
    - inc = 50
    - dec = -30

.. note::

    This is part of a larger project aiming at inverting current sources density (see more at: https://icsd-dev.readthedocs.io/en/latest/)

**References**

Uieda, L., V. C. Oliveira Jr, and V. C. F. Barbosa (2013), Modeling the Earth with Fatiando a Terra, Proceedings of the 12th Python in Science Conference, pp. 91 - 98.

Uieda, L. (2018). Verde: Processing and gridding spatial data using Green's functions. Journal of Open Source Software, 3(29), 957. doi:10.21105/joss.00957

Fedi, M., and M. Pilkington (2012), Understanding imaging methods for potential
field data, Geophysics, 77(1), G13, doi:10.1190/geo2011-0078.1

"""
import matplotlib.pyplot as plt
import numpy as np
import lib.dEXP as dEXP
from lib.dEXP import _fit
import lib.plot_dEXP as pEXP
import lib.set_parameters as para
import examples.magnetic.fwdmag.fwd_mag_sphere as magfwd


#%%
# Create a model using geometric objects from fatiando.mesher
xp, yp, zp, U, shape, p1, p2, coord= magfwd.load_mag_synthetic()
max_elevation=2*max(coord[:,2])
scaled, SI, zp, qorder, nlay, minAlt_ridge, maxAlt_ridge = para.set_par(shape=shape,max_elevation=max_elevation)
interp = True
x_axis='y'

#%%
# Plot field data over a 2d line crossing the anomalies
pEXP.plot_line(xp, yp, U,p1,p2, interp=interp,Xaxis=x_axis)

#%%
# Upward continuation of the field data with discretisation in altitude controlled by the number of layers (nlay) and the maximum elevation desired (max_elevation)
mesh, label_prop = dEXP.upwc(xp, yp, zp, U, shape, 
                 zmin=0, zmax=max_elevation, nlayers=nlay, 
                 qorder=qorder)

plt, cmap = pEXP.plot_xy(mesh, label=label_prop,Xaxis=x_axis)
plt.colorbar(cmap)

#%%
# Ridges identification: plot all extremas obtained via find_peaks function (numpy) for a given altitude
dEXP.ridges_minmax_plot(xp, yp, mesh, p1, p2,
                                      label=label_prop,
                                      method_peak='find_peaks',
                                      fix_peak_nb=5,
                                      Xaxis=x_axis,
                                      showfig=True)  

#%%
# Ridges identification at all levels: plot extremas obtained via find_peaks function (numpy) for all 3 types of extremas familly RI, RII and RIII
dfI,dfII, dfIII = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,
                                      label=label_prop,
                                      method_peak='find_peaks',
                                      fix_peak_nb=5,
                                      Xaxis=x_axis,
                                      showfig=True)  

#%%
# filter ridges using a minimum length criterium and and filter for a specific range of altitude
dfI_f,dfII_f, dfIII_f = dEXP.filter_ridges(dfI,dfII,dfIII,
                                            1,maxAlt_ridge,
                                            minlength=8,rmvNaN=True)
df_f = dfI_f, dfII_f, dfIII_f

#%%
# plot filtered ridges fitted over continuated section
fig, ax = plt.subplots(figsize=(15,3))
pEXP.plot_xy(mesh, label=label_prop, ax=ax) #, ldg=)
pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax,label=False)   
df_fit = dEXP.fit_ridges(df_f) # fit ridges on filtered data

pEXP.plot_ridges_sources(df_fit, ax=ax, z_max_source=-max_elevation*1.2,
                          ridge_type=[0,1,2],ridge_nb=None)

