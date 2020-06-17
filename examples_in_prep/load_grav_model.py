# -*- coding: utf-8 -*-
"""
Created on Thu May 28 08:01:49 2020

@author: Benjamin
"""
import os
import numpy as np 
from fatiando.vis import mpl #, myv
from fatiando import gridder, mesher, utils
from fatiando.gravmag import prism, imaging, transform
from fatiando.vis.mpl import square
import dEXP as dEXP


# %%

def load_grav_fatiando(ZZ=-3.75,shape=(30,30)):

    name = '1000_zbot3000_data'
    xp, yp, zp, U= np.loadtxt('./grav_models/' + name + '.txt', unpack=True)
    shape = (25,25) # data interpolation
    maxdepth=10000 # max depth for upward continuation
    minAlt_ridge = 1000
    maxAlt_ridge = 5000
    SI = 2 # structural index
    zp=0  # initial depth (conditionned upward or downward)
    
    return


# %%
# https://www.pygimli.org/_examples_auto/4_gravimetry_magnetics/plot_mod-gravimetry-integration-2d.html#sphx-glr-examples-auto-4-gravimetry-magnetics-plot-mod-gravimetry-integration-2d-py

import numpy as np
import pygimli as pg
from pygimli.meshtools import createCircle
from pygimli.physics.gravimetry import (gradGZCylinderHoriz,
                                        gradGZHalfPlateHoriz,
                                        gradUCylinderHoriz,
                                        gradUHalfPlateHoriz, solveGravimetry)
    


def load_grav_pygimli_cylinder():

    #!/usr/bin/env python
    # -*- coding: utf-8 -*-
    r"""
    Semianalytical Gravimetry and Geomagnetics in 2D
    ------------------------------------------------
    
    Simple gravimetric and magnetostatic field caluculation using integration approach after :cite:`WonBev1987`.
    
    """

    
    radius = 2.
    depth = 5.
    rho = 1000.0
    
    x = np.arange(-20, 20, 1)
    pnts = np.zeros((len(x), 2))
    pnts[:, 0] = x
    pos = [0, -depth]
    
    fig, ax = pg.plt.subplots(nrows=3, ncols=1, figsize=(12,8), sharex=True)

    # Horizontal cylinder
    circ = createCircle([0, -depth], radius=radius, marker=2, area=0.1,
                        segments=32)
    
    pg.show(circ, ax=ax[0], fillRegion=False)
    
    ga = gradUCylinderHoriz(pnts, radius, rho, pos=pos)
    gza = gradGZCylinderHoriz(pnts, radius, rho, pos=pos)
    g, gz = solveGravimetry(circ, rho, pnts, complete=True)
    
    plot(x, ax[1], ga, gza, ax[2], g, gz, legend=False)
    # ax[0].set_ylim(bottom=-depth*2, top=1)

    labels = ["Horizontal cylinder", "Half plate"]
    for ax, label in zip(ax[:], labels):
        ax.set_title(label)
        ax.set_aspect("equal")
        ax.set_xlim(left=x[0], right=x[-1])
        ax.set_ylim(bottom=-depth*2, top=1)

    # pg.wait()
    
    return ga, gza

def load_grav_pygimli_plate():

    # Half plate
    thickness = 1
    mesh = pg.createGrid(x=np.linspace(0,5000),
                          y=[-depth-thickness/2., -depth+thickness/2.0])
    pg.show(mesh, ax=ax[0,1])
    
    ga = gradUHalfPlateHoriz(pnts, thickness, rho, pos=[0, -depth])
    gza = gradGZHalfPlateHoriz(pnts, thickness, rho, pos=[0, -depth])
    g, gz = solveGravimetry(mesh, rho, pnts, complete=True)
    
    plot(x, ax[1,1], ga, gza, ax[2,1], g, gz)
    
    labels = ["Horizontal cylinder", "Half plate"]
    for ax, label in zip(ax[0], labels):
        ax.set_title(label)
        ax.set_aspect("equal")
        ax.set_xlim(left=x[0], right=x[-1])
        ax.set_ylim(bottom=-depth*2, top=1)
    
    fig.tight_layout()
    
    return g

def plot(x, a1, ga, gza, a2, g, gz, legend=True):
    a1.plot(x, ga[:, 0],  label=r'Analytical $\frac{\partial u}{\partial x}$', c="red")
    a1.plot(x, ga[:, 1],  label=r'Analytical $\frac{\partial u}{\partial z}$', c="blue")

    a1.plot(x, g[:, 0], label=r'Won & Bevis: $\frac{\partial u}{\partial x}$',
            marker='o', linewidth=0, c="red")
    a1.plot(x, g[:, 2], label=r'Won & Bevis: $\frac{\partial u}{\partial z}$',
            marker='o', linewidth=0, c="blue")

    a2.plot(x, gza[:, 0],
            label=r'Analytical $\frac{\partial^2 u}{\partial z,x}$', c="red")
    a2.plot(x, gza[:, 1],
            label=r'Analytical $\frac{\partial^2 u}{\partial z,z}$', c="blue")

    a2.plot(x, gz[:, 0], marker='o', linestyle='',
            label=r'Won & Bevis: $\frac{\partial^2 u}{\partial z,x}$', c="red")
    a2.plot(x, gz[:, 2], marker='o', linestyle='',
            label=r'Won & Bevis: $\frac{\partial^2 u}{\partial z,z}$', c="blue")
    a1.set_xlabel('$x$-coordinate [m]')
    a1.set_ylabel(r'$\frac{\partial u}{\partial (x,z)}$ [mGal]')

    a2.set_xlabel('$x$-coordinate [m]')

    if legend:
        a1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        a2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    

        
    