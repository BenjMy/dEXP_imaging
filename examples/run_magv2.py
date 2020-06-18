"""
Imaging methods for potential fields.

Implements the DEXP method described in Fedi and Pilkington (2012). Application on a 2-sources mag.

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
import plot_dEXP as pEXP
import set_parameters as para

# exemples
import examples.sources_mag.fwd_mag.fwd_mag_sphere as magfwd

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


# Create a model using geometric objects from fatiando.mesher
plt.figure()
plt.imshow(np.random.random((50,50)))
plt.colorbar()
plt.show()

