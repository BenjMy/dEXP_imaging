"""
Imaging methods for potential fields.

Implements the DEXP method described in Fedi and Pilkington (2012).

.. note::

    This is part of a larger project aiming at inverting current sources density (see more at: https://icsd-dev.readthedocs.io/en/latest/)


**References**

Fedi, M., and M. Pilkington (2012), Understanding imaging methods for potential
field data, Geophysics, 77(1), G13, doi:10.1190/geo2011-0078.1

----
"""

import math
import numpy as np
import matplotlib.pyplot as plt

from fatiando.gravmag import imaging, transform
from fatiando.gravmag.imaging import _makemesh
from fatiando import gridder, mesher, utils

from scipy.interpolate import UnivariateSpline
from scipy.signal import find_peaks # useful for ridges detection


def cor_field_B(x,y,z,u,B,rho=100):
    """
    Calculates the potential field (electric) produced by a current injection in B (return electrode) for a
    given homogeneous electrical resistivity rho
    """
    I = 1 # injected current (A)
    num = rho*I
    dist = np.sqrt((x-B[0])**2 + (y-B[1])**2 + (z-B[2])**2) # distance between B and potential electrodes
    den = 2*math.pi*dist
    u_B = num/den

    u_cor = u - u_B # correct total measured potential from influence of B

    plt.figure()
    plt.subplot(1,3,1)
    plt.tricontourf(x, y, dist, 50, cmap='RdBu_r')
    cbar = plt.colorbar(orientation='vertical', aspect=50)
    cbar.set_label('Distance from B (m)')
    plt.tight_layout()
    plt.axis('square')

    plt.subplot(1,3,2)
    plt.tricontourf(x, y, u_B, 50, cmap='RdBu_r')
    cbar = plt.colorbar(orientation='vertical', aspect=50)
    cbar.set_label('$u_{B}$')
    plt.tight_layout()
    plt.axis('square')

    plt.subplot(1,3,3)
    plt.tricontourf(x, y, u_cor, 50, cmap='RdBu_r')
    cbar = plt.colorbar(orientation='vertical', aspect=50)
    cbar.set_label('$u = U_{T} - u_{B}$ (V)')
    plt.tight_layout()
    plt.axis('square')

    return u_cor   



def ridges_0(x, y, mesh, p1, p2, qorder=0, z=0):
    """
    Text here

    Parameters:

    * mesh
        The upward continuated field mesh (of order-q derivative)

    Returns:

    * BB : 
        Text here - Panda dataframe containing RI, RII and RII

    """
    
    # This way, if z is not an array, it is now
    z = z * np.ones_like(x)
    
    RI_0 = [] # zeros of the first horizontal derivative of the potential field
    RII_0 = [] # zeros of the first vertical derivative of the potential field
    RIII_0 = [] # zeros of the potential field
    
    depths = mesh.get_zs()[:-2]
    upw_u = mesh.props['csd']
    upw_u = np.reshape(upw_u, [mesh.shape[0], mesh.shape[1]*mesh.shape[2]])      

    for i, depth in enumerate(depths - z[0]):
        
        
        fd1z = []
        fd1x = []
        fd = []
        
        upw_u_l = upw_u[i+1,:]    
        xx, yy, distance, p_up_f = gridder.profile(x, y, upw_u_l, p1, p2, 1000)
        
        
        fd= UnivariateSpline(xx,p_up_f, s=0)
        if fd.roots().any():
            RIII_0.append([depth, np.array(fd.roots())])
        else:
            RIII_0.append([depth, []])
        
        # 1st vertical derivate of the continued field
        up_f_d1z = transform.derivz(x, y, upw_u_l,(mesh.shape[1],mesh.shape[1]),order=1)
        xx, yy, distance, p_up_f_d1z = gridder.profile(x, y, up_f_d1z, p1, p2, 1000)
        fd1z= UnivariateSpline(xx,p_up_f_d1z, s=0)
        if fd1z.roots().any():
            RII_0.append([depth, np.array(fd1z.roots())])
        else:
            RII_0.append([depth, []])
            
        # 1st horizontal derivate of the continued field
        up_f_d1x = transform.derivx(x, y, upw_u_l,(mesh.shape[1],mesh.shape[1]),order=1)
        xx, yy, distance, p_up_f_d1x = gridder.profile(x, y, up_f_d1x, p1, p2, 1000)
        fd1x= UnivariateSpline(xx,p_up_f_d1x, s=0)
        if fd1x.roots().any():
            RI_0.append([depth, np.array(fd1x.roots())])
        else:
            RI_0.append([depth, []])
    

    # Once extreme points are determined at different altitudes, 
    # ridges can be obtained by linking each of them, at a given altitude, 
    # to the nearest one computed at the altitude just above.


    # import pandas as pd
    # ## Create panda dataframe to merge all the ridges
    # RIpd= pd.DataFrame(data=RI_0[0][1],    # values
    #             index=RI_0[0][0] - zp[0],    # depths column as index
    #               columns='RI')  # 1st row as the column names
    
    RI_0= np.array(RI_0)
    RII_0= np.array(RII_0)
    RIII_0= np.array(RIII_0)

    return RIII_0


def ridges_minmax(x, y, mesh, p1, p2, qorder=0, z=0):
    """
    Text here

    Parameters:

    * mesh
        The upward continuated field mesh (of order-q derivative)

    Returns:

    * BB : 
        Text here - Panda dataframe containing RI, RII and RII

    """
    
    # This way, if z is not an array, it is now
    z = z * np.ones_like(x)
       
    RI_minmax = [] # minmax of the first horizontal derivative of the potential field
    RII_minmax = [] # minmax of the first vertical derivative of the potential field
    RIII_minmax = [] # minmax of the potential field
   
    depths = mesh.get_zs()[:-1]
    
    if mesh.props['upwc'].any():
        upw_u = mesh.props['upwc']
    else:
        print('Missing upward continuation step')
    upw_u = np.reshape(upw_u, [mesh.shape[0], mesh.shape[1]*mesh.shape[2]])      

    for i, depth in enumerate(depths - z[0]):
        
        
        fd1z = []
        fd1x = []
        fd = []
        
        upw_u_l = upw_u[i,:]    
        xx, yy, distance, p_up_f = gridder.profile(x, y, upw_u_l, p1, p2, 1000)
        
        Max_peaks, _ = find_peaks(p_up_f)
        Min_peaks, _ = find_peaks(-p_up_f)
    
        if np.array(Max_peaks).any():
            RIII_minmax.append([depth, xx[Max_peaks]])
        else:
            RIII_minmax.append([depth, []])
        
        # 1st vertical derivate of the continued field
        up_f_d1z = transform.derivz(x, y, upw_u_l,(mesh.shape[1],mesh.shape[1]),order=1)
        xx, yy, distance, p_up_f_d1z = gridder.profile(x, y, up_f_d1z, p1, p2, 1000)
        Max_peaks, _ = find_peaks(p_up_f_d1z)
        Min_peaks, _ = find_peaks(-p_up_f_d1z)

        if np.array(Max_peaks).any():
            RIII_minmax.append([depth, xx[Max_peaks]])
        else:
            RII_minmax.append([depth, []])
            
        # 1st horizontal derivate of the continued field
        up_f_d1x = transform.derivx(x, y, upw_u_l,(mesh.shape[1],mesh.shape[1]),order=1)
        xx, yy, distance, p_up_f_d1x = gridder.profile(x, y, up_f_d1x, p1, p2, 1000)
        Max_peaks, _ = find_peaks(p_up_f_d1x)
        Min_peaks, _ = find_peaks(-p_up_f_d1x)

        if np.array(Max_peaks).any():
            RI_minmax.append([depth, xx[Max_peaks]])
        else:
            RI_minmax.append([depth, []])
            
    # import pandas as pd
    # ## Create panda dataframe to merge all the ridges
    # RIpd= pd.DataFrame(data=RI_0[0][1],    # values
    #             index=RI_0[0][0] - zp[0],    # depths column as index
    #               columns='RI')  # 1st row as the column names
    R = [np.array(RI_minmax), np.array(RII_minmax), np.array(RIII_minmax)]

    R_fit = _build_ridge(RI_minmax,RII_minmax,RIII_minmax)
    
    return R, R_fit


def geom_Z0():
    """
    Find intersection of ridges

    Parameters:

    * a
        Text here

    Returns:

    * BB : 
        Text here

    """
    return

def scalFUN():
    """
    Analysis (Euler deconvolution of ridges)

    Parameters:

    * a
        Text here

    Returns:

    * BB : 
        Text here

    """
    return Tau, q, SI

def upwc(x, y, z, data, shape, zmin, zmax, nlayers, qorder=0):
    """
    Upward continuation model (Fedi, 2012).

    Calculates the upward continuation for given potential field data on a
    **regular grid**. 
    Parameters:

    * x, y : 1D-arrays
        The x and y coordinates of the grid points
    * z : float or 1D-array
        The z coordinate of the grid points
    * data : 1D-array
        The potential field at the grid points
    * shape : tuple = (ny, nx)
        The shape of the grid
    * zmin, zmax : float
        The top and bottom, respectively, of the region where the physical
        property distribution is calculated
    * nlayers : int
        The number of layers used to divide the region where the physical
        property distribution is calculated
    * qorder : float
        The order of the vertical derivative

    Returns:

    * mesh : :class:`fatiando.mesher.PrismMesh`
        The estimated physical property distribution set in a prism mesh (for
        easy 3D plotting)

    """
    mesh = _makemesh(x, y, shape, zmin, zmax, nlayers)
    # This way, if z is not an array, it is now
    z = z * np.ones_like(x)
    # Remove the last z because I only want depths to the top of the layers
    depths = mesh.get_zs()[:-1]

    upw_f = []
    upw_f_dq = []
    # Offset by the data z because in the paper the data is at z=0
    for depth in (depths - z[0]):

        # continued field calculation
        upw_fhi= transform.upcontinue(x, y, data, shape, depth)

        # qorder vertical derivate of the continued field
        upw_f_dqhi = transform.derivz(x, y, upw_fhi, shape,order=qorder)
        
        upw_f.extend(upw_fhi)
        upw_f_dq.extend(upw_f_dqhi)

    mesh.addprop('upwc', np.array(upw_f))
    mesh.addprop('upw_f_d'+str(qorder), np.array(upw_f_dq))
    return mesh

def dEXP(x, y, z, data, shape, zmin, zmax, nlayers, qorder=0, SI=1):
    """
    DEXP model (Fedi, 2012).

    Calculates a physical property distribution given potential field data on a
    **regular grid**. Uses depth weights.
    Parameters:

    * x, y : 1D-arrays
        The x and y coordinates of the grid points
    * z : float or 1D-array
        The z coordinate of the grid points
    * data : 1D-array
        The potential field at the grid points
    * shape : tuple = (ny, nx)
        The shape of the grid
    * zmin, zmax : float
        The top and bottom, respectively, of the region where the physical
        property distribution is calculated
    * nlayers : int
        The number of layers used to divide the region where the physical
        property distribution is calculated
    * qorder : float
        The order of the vertical derivative
    * SI : float
        The structural index

    Returns:

    * mesh : :class:`fatiando.mesher.PrismMesh`
        The estimated physical property distribution set in a prism mesh (for
        easy 3D plotting)

    """
    mesh = _makemesh(x, y, shape, zmin, zmax, nlayers)
    # This way, if z is not an array, it is now
    z = z * np.ones_like(x)
    # Remove the last z because I only want depths to the top of the layers
    depths = mesh.get_zs()[:-1]
    weights = (np.abs(depths)) ** ((SI+qorder)/2)
    csd = []
    # Offset by the data z because in the paper the data is at z=0
    for depth, weight in zip(depths - z[0], weights):

        # continued field calculation
        upw_f= transform.upcontinue(x, y, data, shape, depth)

        # qorder vertical derivate of the continued field
        upw_f_dq = transform.derivz(x, y, upw_f, shape,order=qorder)
        
        # the continued field weigted (=DEXP)
        upw_f_dq_w= upw_f_dq*weight
        csd.extend(upw_f_dq_w)

    mesh.addprop('dexp', np.array(csd))
    return mesh, np.array(csd)

# def auto_dEXP():


def _build_ridge(RI_minmax,RII_minmax,RIII_minmax):
    # Once extreme points are determined at different altitudes, 
    # ridges can be obtained by linking each of them, at a given altitude, 
    # to the nearest one computed at the altitude just above.
    
    return