import os

from fatiando.vis import mpl #, myv
from fatiando import gridder, mesher, utils
from fatiando.gravmag import prism, imaging, transform
from fatiando.vis.mpl import square

# my own functions
import dEXP as dEXP
from dEXP import _fit
import plot_dEXP as pEXP

# exemples
import exemples.fwd_mag_sphere as magfwd
import exemples.load_grav_model as grav
import exemples.load_MALM_model as MALM

import set_parameters as para

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 15

# icsd functions
from icsd3d.importers.read import load_obs, load_geom

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
        

#%% ------------------------------- ridges identification
import Find_peaks_Du_et_al_2006 as fpeak
import Find_peaks as fp
from Find_peaks_Du_et_al_2006 import _boolrelextrema, _identify_ridge_lines, _filter_ridge_lines

upw_u_test = mesh.props[label_prop]
upw_u_test = np.reshape(upw_u_test, [mesh.shape[0], mesh.shape[1]*mesh.shape[2]])   

peak_width = np.linspace(400,1000)

fp.find_peaks_U(upw_u_test,widths=peak_width, max_distances=None,gap_thresh=None)


