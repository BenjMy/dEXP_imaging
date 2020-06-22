"""
Test example for mag data analysis using pyDEXP
--------------------------------------------------------------------

Plot 2d randomn data 
Plot field data over a 2d line

"""
import matplotlib.pyplot as plt
import numpy as np
import dEXP as dEXP
from dEXP import _fit
import plot_dEXP as pEXP
import set_parameters as para
import examples.sources_mag.fwd_mag_sphere as magfwd


# Create a model using geometric objects from fatiando.mesher
xp, yp, zp, U, shape, p1, p2, coord= magfwd.load_mag_synthetic()
max_elevation=2*max(coord[:,2])
scaled, SI, zp, qorder, nlay, minAlt_ridge, maxAlt_ridge = para.set_par(shape=shape,max_elevation=max_elevation)
interp = True

# Plot field data over a 2d line
pEXP.plot_line(xp, yp, U,p1,p2, interp=interp)

# Upward continuation of the field data
mesh, label_prop = dEXP.upwc(xp, yp, zp, U, shape, 
                 zmin=0, zmax=max_elevation, nlayers=nlay, 
                 qorder=qorder)

plt, cmap = pEXP.plot_xy(mesh, label=label_prop)
plt.colorbar(cmap)

# Ridges identification
dEXP.ridges_minmax_plot(xp, yp, mesh, p1, p2,
                                      label=label_prop,
                                      method_peak='find_peaks')  

# or find_peaks or peakdet or spline_roots
dfI,dfII, dfIII = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,
                                      label=label_prop,
                                      method_peak='find_peaks')  