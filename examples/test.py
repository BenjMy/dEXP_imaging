"""
Forward modeling gravity data using spheres in Cartesian coordinates
--------------------------------------------------------------------

The :mod:`fatiando.gravmag` has many functions for forward modeling gravity and
magnetic data. Here we'll show how to build a model out of spheres and
calculate the gravitational attraction and it's gradients in Cartesian
coordinates.

"""
import matplotlib.pyplot as plt
import numpy as np
import dEXP as dEXP
from dEXP import _fit
import plot_dEXP as pEXP
import set_parameters as para


# Create a model using geometric objects from fatiando.mesher
plt.figure()
plt.imshow(np.random.random((50,50)))
plt.colorbar()
plt.show()

xp, yp, zp, U, shape, p1, p2, coord= magfwd.load_mag_synthetic()
max_elevation=2*max(coord[:,2])
scaled, SI, zp, qorder, nlay, minAlt_ridge, maxAlt_ridge = para.set_par(shape=shape,max_elevation=max_elevation)
interp = True

pEXP.plot_line(xp, yp, U,p1,p2, interp=interp)

