# -*- coding: utf-8 -*-
"""
Miscellaneous utility functions.

@author: Benjamin
"""
import numpy as np 
import matplotlib.pyplot as plt
import math 
from numpy import inf
from fatiando import gridder, mesher, utils
from matplotlib.backends.backend_pdf import PdfPages

#%% Functions SPECIFIC FOR MISE-A-LA-MASSE data preparation to dEXP transformation 

def cor_field_B(x,y,z,u,B,rho=100,**kwargs):
    """
    Calculates the potential field (electric) produced by a current injection in B (return electrode) for a
    given homogeneous electrical resistivity rho.

    Parameters
    ----------
    * x, y : 1D-arrays
        The x and y coordinates of the grid points
    * z : float or 1D-array
        The z coordinate of the grid points
    * u : 1D-array
        The potential field at the grid points
    * Bu : 1D-array
        The position of the return electrode B
    * rho : int
        The estimated electrical resistivity of the medium
    Returns
    -------
    u_cor : 1D-arrays
        The corrected from the B return electrode influence potential field at the grid points
    """
    I = 1 # injected current (A)
    num = rho*I
    dist = np.sqrt((x-B[0])**2 + (y-B[1])**2 + (z-B[2])**2) # distance between B and potential electrodes
    den = 2*math.pi*dist
    u_B = num/den
    
    print(u_B)
    u_B[u_B == inf] = 0
    
    u_cor = u + u_B # correct total measured potential from influence of B


    for key, value in kwargs.items():
        if key == 'vmin':
            vmin = value
    for key, value in kwargs.items():
        if key == 'vmax':
            vmax = value  

    plt_2 = None
    for key, value in kwargs.items():
        if key == 'plt_2':
            plt_2 = value
    print(plt_2)
            
    # norm = plt.Normalize(vmax=vmax, vmin=vmin)
    # v = np.linspace(vmin, vmax, 15, endpoint=True)

    vmin = min(dist)
    vmax = max(dist)
    plt.figure() # figsize=(20,10)
    plt.subplot(2,2,1)
    plt.scatter(x, y, c=dist, cmap='RdBu_r')
    cbar = plt.colorbar(orientation='vertical')
    cbar.set_label('Distance from B (m)')
    # plt.tight_layout()
    plt.axis('square')
    
    vmin = min(u_B)
    vmax = max(u_B)
    plt.subplot(2,2,2)
    plt.scatter(x, y, c=u_B, cmap='RdBu_r')
    cbar = plt.colorbar(orientation='vertical')
    cbar.set_label('$u_{B}$ (V)')
    # plt.scatter(B[0],B[1], marker='v', color='red')
    # plt.tight_layout()
    if plt_2 is not None:
        plt.plot(plt_2[0,:],plt_2[1,:],'*-')
    plt.axis('square')

    # vmin = min(u)
    # vmax = max(u)
    plt.subplot(2,2,3)
    plt.scatter(x, y, c=u, cmap='RdBu_r')
    cbar = plt.colorbar(orientation='vertical')
    cbar.set_label('$u = U_{T} $(V)')
    # plt.scatter(B[0],B[1], marker='v', color='red')
    # plt.tight_layout()
    plt.axis('square')
    
    plt.subplot(2,2,4)
    plt.scatter(x, y, c=u_cor, cmap='RdBu_r')
    cbar = plt.colorbar(orientation='vertical')
    cbar.set_label('$u = U_{T} - u_{B}$ (V)')
    # plt.scatter(B[0],B[1], marker='v', color='red')
    # plt.tight_layout()
    plt.axis('square')

    plt.figure()
    plt.scatter(x, y, c=u-u_cor, cmap='RdBu_r')
    cbar = plt.colorbar(orientation='vertical')
    cbar.set_label('$u = U_{T} - u_{u_cor}$ (V)')
    plt.axis('square')
    
    return u_cor   


#%% Geometrical correction/filtering

def _point_against_line(xp,yp,p1,p2):
# https://stackoverflow.com/questions/45766534/finding-cross-product-to-find-points-above-below-a-line-in-matplotlib

    isabove = lambda p, p1,p2: np.cross(p-p1, p2-p1) < 0
    
    # p1 = np.array([1,1])
    # p2 = np.array([4,3])
    
    fig, (ax,ax2) = plt.subplots(ncols=2, sharex=True, sharey=True)

    # p = np.random.rand(10,2)*5
    p = np.vstack([xp,yp]).T
    print(p)
    
    ax2.plot([p1[0],p2[0]],[p1[1],p2[1]], marker="o", color="k")
    ax2.scatter(p[:,0],p[:,1], c=isabove(p,p1,p2), cmap="bwr", vmin=0, vmax=1)
    
    
    # ax.set_xlim(0,6)
    # ax.set_ylim(0,6)
    
    plt.show()

    return

def mirror(xp,yp, data, p1, p2, mirror='above'):
    
    # put zero over/under a certain line
    _point_against_line(xp,yp, p1, p2)
    
    # replicate values to mirror
    
    # put zero to 
    # data_sym = np.zeros()
    
    # for row in range(data.shape[0]):
    #     for col in range(data.shape[1]):
    #         data_sym = data()
            
    return


def perp_p1p2(p1,p2, offset=0):
    
    midX=(p1[0]+p2[0])/2
    midY=(p1[1]+p2[1])/2 
    
    plt.scatter(midX,midY,c='red')

    # midX=(p1[0]+offset+p2[0]+offset)/2
    # midY=(p1[1]+offset+p2[1]+offset)/2 
    
    # plt.scatter(midX,midY,c='green')


    new_p2 = [midX-p2[1]+p1[1], midY + p2[0]-p1[0]]
    new_p1 = [midX+p2[1]-p1[1], midY - p2[0]+p1[0]]

    # new_p2 = [midX-p2[1]+p1[1]/2, midY + p2[0]-p1[0]/2]
    # new_p1 = [midX+p2[1]-p1[1]/2, midY - p2[0]+p1[0]/2]
    
    p12x=[p1[0],p2[0]]
    p12y=[p1[1],p2[1]]
    
    p12x_new=[new_p1[0],new_p2[0]]
    p12y_new=[new_p1[1],new_p2[1]]
    
    plt.plot(p12x,p12y,c='red')
    plt.plot(p12x_new,p12y_new,c='green')
    plt.axis('square')
    
    
    return new_p1, new_p2


def load_surfer(fname, dtype='float64'):
    """
    Read data from a Surfer ASCII grid file.

    Surfer is a contouring, griding and surface mapping software
    from GoldenSoftware. The names and logos for Surfer and Golden
    Software are registered trademarks of Golden Software, Inc.

    http://www.goldensoftware.com/products/surfer

    Parameters:

    * fname : str
        Name of the Surfer grid file
    * dtype : numpy dtype object or string
        The type of variable used for the data. Default is numpy.float64. Use
        numpy.float32 if the data are large and precision is not an issue.

    Returns:

    * data : dict
        The data in a dictionary with some metadata:

        * ``'file'`` : string
            The name of the original data file
        * ``'shape'`` : tuple
            (nx, ny), the number of grid points in x (North) and y (East)
        * ``'area'`` : tuple
            (x1, x2, y1, y2), the spacial limits of the grid.
        * ``'x'`` : 1d-array
            Value of the North-South coordinate of each grid point.
        * ``'y'`` : 1d-array
            Value of the East-West coordinate of each grid point.
        * ``'data'`` : 1d-array
            Values of the field in each grid point. Field can be for example
            topography, gravity anomaly, etc. If any field values are >=
            1.70141e+38 (Surfers way of marking NaN values), the array will be
            masked at those values (i.e., ``'data'`` will be a numpy masked
            array).

    """
    # Surfer ASCII grid structure
    # DSAA            Surfer ASCII GRD ID
    # nCols nRows     number of columns and rows
    # xMin xMax       X min max
    # yMin yMax       Y min max
    # zMin zMax       Z min max
    # z11 z21 z31 ... List of Z values
    with open(fname) as input_file:
        # DSAA is a Surfer ASCII GRD ID (discard it for now)
        input_file.readline()
        # Read the number of columns (ny) and rows (nx)
        ny, nx = [int(s) for s in input_file.readline().split()]
        shape = (nx, ny)
        # Our x points North, so the first thing we read is y, not x.
        ymin, ymax = [float(s) for s in input_file.readline().split()]
        xmin, xmax = [float(s) for s in input_file.readline().split()]
        area = (xmin, xmax, ymin, ymax)
        dmin, dmax = [float(s) for s in input_file.readline().split()]
        field = np.fromiter((float(s)
                             for line in input_file
                             for s in line.split()),
                            dtype=dtype)
        nans = field >= 1.70141e+38
        if np.any(nans):
            field = np.ma.masked_where(nans, field)
        err_msg = "{} of data ({}) doesn't match one from file ({})."
        assert np.allclose(dmin, field.min()), err_msg.format('Min', dmin,
                                                              field.min())
        assert np.allclose(dmax, field.max()), err_msg.format('Max', dmax,
                                                              field.max())
    x, y = gridder.regular(area, shape)
    data = dict(file=fname, shape=shape, area=area, data=field, x=x, y=y)
    return data


def multipage(filename, figs=None, dpi=200):
    pp = PdfPages(filename)
    if figs is None:
        figs = [plt.figure(n) for n in plt.get_fignums()]
    for fig in figs:
        fig.savefig(pp, format='pdf')
    pp.close()