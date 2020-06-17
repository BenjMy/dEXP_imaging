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

# Create a model using geometric objects from fatiando.mesher
plt.figure()
plt.imshow(np.random.random((50,50)))
plt.colorbar()
plt.show()
