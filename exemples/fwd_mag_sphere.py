"""
Forward modeling magnetic data using spheres in Cartesian coordinates
-----------------------------------------------------------------------

The :mod:`fatiando.gravmag` has many functions for forward modeling gravity and
magnetic data. Here we'll show how to build a model out of spheres and
calculate the total field magnetic anomaly and the 3 components of the magnetic
induction.
"""
from __future__ import division, print_function
from fatiando import mesher, gridder, utils
from fatiando.gravmag import sphere
import matplotlib.pyplot as plt
import numpy as np

# %%
# Create a model using geometric objects from fatiando.mesher
# Each model element has a dictionary with its physical properties.
# The spheres have different total magnetization vectors (total = induced +
# remanent + any other effects). Notice that the magnetization has to be a
# vector. Function utils.ang2vec converts intensity, inclination, and
# declination into a 3 component vector for easier handling.

def anomag_model(coords,radii=1.5e3,inc=50, dec=-30):
    
    model =[]
    for coord, radius in zip(coords, radii):
  
        model_tmp = mesher.Sphere(x=coord[0], y=coord[1], z=coord[2], radius=radius,
                      props={'magnetization': utils.ang2vec(1, inc=50, dec=-30)})
            
        model.append(model_tmp)
    
    return model

def fwd_model(model, shape = (300, 300), area = [0, 30e3, 0, 30e3]):

    # Set the inclination and declination of the geomagnetic field.
    inc, dec = -10, 13
    
    # Create a regular grid at a constant height
    shape = shape
    area = area
    x, y, z = gridder.regular(area, shape, z=-10)
    
    field = ['Total field Anomaly (nt)', sphere.tf(x, y, z, model, inc, dec)]
    
    return x, y, z, field, shape

def plot_model(x, y, field, shape):

    # Make maps of all fields calculated
    fig = plt.figure()
    ax = plt.gca()
    plt.rcParams['font.size'] = 10
    X, Y = x.reshape(shape)/1000, y.reshape(shape)/1000     
    field_name, data = field
    scale = np.abs([data.min(), data.max()]).max()
    ax.set_title(field_name)
    plot = ax.pcolormesh(Y, X, data.reshape(shape), cmap='RdBu_r',
                         vmin=-scale, vmax=scale)
    plt.colorbar(plot, ax=ax, aspect=30, pad=0)
    ax.set_xlabel('y (km)')
    ax.set_ylabel('x (km)')
    plt.tight_layout(pad=0.5)
    plt.show()

def load_mag_synthetic():
    
    A= [10e3,10e3,2e3]
    B= [25e3,10e3,1e3]
    coord= np.array([A,B])
    radius = 1.5e3*np.ones(len(coord))
    modelmag = anomag_model(coord,radii=radius,inc=50, dec=-30)
    xp, yp, zp, field2d, shape = fwd_model(modelmag, shape = (300, 300),area = [0, 30e3, 0, 30e3])
    plot_model(xp, yp, field2d, shape)
    
    
    U = field2d[1]
    p1, p2 = [min(xp), 10e3], [max(xp), 10e3]
        
    return xp, yp, zp, U, shape, p1, p2, coord
# %% pygimli exemple
