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

import pandas as pd
from scipy.optimize import curve_fit


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

    plt.figure(figsize=(20,10))
    plt.subplot(2,2,1)
    plt.tricontourf(x, y, dist, 50, cmap='RdBu_r')
    cbar = plt.colorbar(orientation='vertical', aspect=50)
    cbar.set_label('Distance from B (m)')
    plt.tight_layout()
    plt.axis('square')

    plt.subplot(2,2,2)
    plt.tricontourf(x, y, u_B, 50, cmap='RdBu_r')
    cbar = plt.colorbar(orientation='vertical', aspect=50)
    cbar.set_label('$u_{B}$ (V)')
    plt.tight_layout()
    plt.axis('square')
    
    plt.subplot(2,2,3)
    plt.tricontourf(x, y, u, 50, cmap='RdBu_r')
    cbar = plt.colorbar(orientation='vertical', aspect=50)
    cbar.set_label('$u = U_{T} $(V)')
    plt.tight_layout()
    plt.axis('square')
    
    plt.subplot(2,2,4)
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


def ridges_minmax(x, y, mesh, p1, p2, qorder=0, z=0, label='upwc', **kwargs):
    """
    Form a multiridge set
    RI and RII : zeros of the first horizontal and first vertical derivatives of the potential field
    RIII :zeros of the potential field itself

    .. note:: ridges generated by isolated simple sources point, line, sheet, and contact are straight lines defined by the zeros of a potential 
    field and its horizontal and vertical derivatives at all measured or computed levels

    Parameters:

    * mesh
        The upward continuated field mesh (of order-q derivative). 
        Read the property label 'upwc' of the mesh by default
    
    **kwargs
        prominence for peak detection

    Returns:

    * BB : 
        Panda dataframe containing RI, RII and RII

    """
    prom = 0.1
    # This way, if z is not an array, it is now
    z = z * np.ones_like(x)
       
    RI_minmax = [] # minmax of the first horizontal derivative of the potential field
    RII_minmax = [] # minmax of the first vertical derivative of the potential field
    RIII_minmax = [] # minmax of the potential field
   
    depths = mesh.get_zs()[:-1]
    
    if label not in mesh.props:
        raise ValueError("mesh doesn't have a '%s' property." % (label))
    else:
        upw_u = mesh.props[label]
    upw_u = np.reshape(upw_u, [mesh.shape[0], mesh.shape[1]*mesh.shape[2]])      

    # for i, depth in enumerate(depths - z[0]):
    #     upw_u_l = upw_u[i,:]    # analysing extrema layers by layers
    #     xx, yy, distance, p_up_f = gridder.profile(x, y, upw_u_l, p1, p2, 1000)
   
        
    for i, depth in enumerate(depths - z[0]): # Loop over RIII extremas
        upw_u_l = upw_u[i,:]    # analysing extrema layers by layers      
        xx, yy, distance, p_up_f = gridder.profile(x, y, upw_u_l, p1, p2, 1000)
        Max_peaks, _ = find_peaks(p_up_f)
        Min_peaks, _ = find_peaks(-p_up_f)
        MinMax_peaks= np.hstack([Min_peaks,Max_peaks])


        if np.array(MinMax_peaks).any():
            RIII_minmax.append(np.hstack([[depth], xx[MinMax_peaks]]))
        else:
            RIII_minmax.append(np.hstack([[depth],[]]))
        
        if i == 3:
            plt.figure()
            plt.subplot(3,1,1)
            plt.plot(xx,p_up_f,label='u')
            for ind in range(len(Max_peaks)):
                plt.scatter(xx[Max_peaks[ind]],p_up_f[Max_peaks[ind]],color='g')
            plt.legend()
 
    for i, depth in enumerate(depths - z[0]): # Loop over RII extremas
        upw_u_l = upw_u[i,:]    # analysing extrema layers by layers        

        # 1st vertical derivate of the continued field
        up_f_d1z = transform.derivz(x, y, upw_u_l,(mesh.shape[1],mesh.shape[1]),order=1)
        xx, yy, distance, p_up_f_d1z = gridder.profile(x, y, up_f_d1z, p1, p2, 1000)
        Max_peaks, _ = find_peaks(p_up_f_d1z)
        Min_peaks, _ = find_peaks(-p_up_f_d1z)

        MinMax_peaks= np.hstack([Min_peaks,Max_peaks])

        if np.array(MinMax_peaks).any():
            RII_minmax.append(np.hstack([[depth], xx[MinMax_peaks]]))
        else:
            RII_minmax.append(np.hstack([[depth],[]]))

        if i == 3:
            plt.subplot(3,1,2)
            plt.plot(xx,p_up_f_d1z,label='dz')
            for ind in range(len(Max_peaks)):
                plt.scatter(xx[Max_peaks[ind]],p_up_f_d1z[Max_peaks[ind]],color='b')
            plt.legend()

    for i, depth in enumerate(depths - z[0]): # Loop for RII extremas
        upw_u_l = upw_u[i,:]    # analysing extrema layers by layers   
        # 1st horizontal derivate of the continued field
        up_f_d1x = transform.derivx(x, y, upw_u_l,(mesh.shape[1],mesh.shape[1]),order=1)
        xx, yy, distance, p_up_f_d1x = gridder.profile(x, y, up_f_d1x, p1, p2, 1000)
        Max_peaks, _ = find_peaks(p_up_f_d1x)
        Min_peaks, _ = find_peaks(-p_up_f_d1x)
        # Max_peaks, _ = find_peaks(p_up_f_d1x)
        # Min_peaks, _ = find_peaks(-p_up_f_d1x)       
        MinMax_peaks= np.hstack([Min_peaks,Max_peaks])
        print(depth,len(Min_peaks),len(Max_peaks))
        
        if np.array(MinMax_peaks).any():
            RI_minmax.append(np.hstack([[depth], xx[MinMax_peaks]]))
        else:
            RI_minmax.append(np.hstack([[depth],[]]))

        if i == 3:
            plt.subplot(3,1,3)
            plt.plot(xx,p_up_f_d1x,label='dx')
            for ind in range(len(Max_peaks)):
                plt.scatter(xx[Max_peaks[ind]],p_up_f_d1x[Max_peaks[ind]],color='r')
                plt.scatter(xx[Min_peaks[ind]],p_up_f_d1x[Min_peaks[ind]],color='r')
            plt.legend()
            
    # R = [np.array(RI_minmax), np.array(RII_minmax), np.array(RIII_minmax)]
    dfI,dfII, dfIII = _ridges_2_df(RI_minmax, RII_minmax, RIII_minmax)
    # R_fit = _build_ridge(RI_minmax,RII_minmax,RIII_minmax)
    
    return dfI,dfII, dfIII #, R, R_fit


def filter_ridges(dfI,dfII,dfIII,minDepth,maxDepth, minlength=3, rmvNaN=False):
    """
    Filter non-rectiligne ridges (denoising)

    Parameters:

    * minDepth
        Text here
    * maxDepth

    Returns:

    * BB : 
        Text here

    """
    
    if rmvNaN == True:
        if dfI.isnull().values.any():
            print('NaN or Inf detected - trying to remove')
            dfI.dropna(axis=1, inplace=True) # remove collumns
            dfI = dfI[~dfI.isin([np.nan, np.inf, -np.inf]).any(1)] #remove lines
        if dfII.isnull().values.any():
            dfII.dropna(axis=1, inplace=True) # remove collumns
            dfII = dfII[~dfII.isin([np.nan, np.inf, -np.inf]).any(1)] #remove lines
        if dfIII.isnull().values.any():
            dfIII.dropna(axis=1, inplace=True) # remove collumns
            dfIII = dfIII[~dfIII.isin([np.nan, np.inf, -np.inf]).any(1)] #remove lines
        
    
    if dfI is not None:
        dfI = dfI.loc[(dfI['elevation'] > minDepth) & (dfI['elevation'] < maxDepth)]
    if dfII is not None:
        dfII = dfII.loc[(dfII['elevation'] > minDepth) & (dfII['elevation'] < maxDepth)]
    if dfI is not None:
        dfIII = dfIII.loc[(dfIII['elevation'] > minDepth) & (dfIII['elevation'] < maxDepth)]
    
    # check length of ridges (remove column if less than 3 points)
    if dfI is not None:
        smallCol = dfI.count()
        idrmv = np.where(smallCol<minlength)[0].tolist()
        dfI = dfI.drop(dfI.columns[idrmv], axis=1)
    if dfII is not None:
        smallCol = dfII.count()
        idrmv = np.where(smallCol<minlength)[0].tolist()
        dfII = dfII.drop(dfII.columns[idrmv], axis=1)
    if dfIII is not None:
        smallCol = dfIII.count()
        idrmv = np.where(smallCol<minlength)[0].tolist()
        dfIII = dfIII.drop(dfIII.columns[idrmv], axis=1)


    return dfI, dfII, dfIII


def fit_ridges(df):
    """
    Fit ridges and return points and fit equations to plot

    Parameters:

    * df
        dataframe including all tree types of ridges

    Returns:

    * BB : 
        points and fit equations to plot

    """
    if len(df)==1:
        df = [df]
    
    if len(df[0])==0:
        raise ValueError("No data to fit")
    
    
    df_Rfit = []
    # plt.figure()
    for r_type in range(len(df)): # loop over ridges type I, II, III
        fit_ridges_all = []
        lable = []
        cols = []
    
        for k in enumerate(df[r_type].columns[1:]): # loop over ridges of the same familly
            if np.count_nonzero(np.diff(df[r_type][k[1]]))<5: # check if ridge is vertical
                # print('vertical ridge type:' + str(r_type) + ' / ridgenb:' + k[1])
                fit_name = 'R'+ str(r_type) + ' Vert.' +  k[1]
                y_fit = np.linspace(-max(df[r_type]['elevation'])*2,max(df[r_type]['elevation']),100)                                       
                x_fit = df[r_type][k[1]].iloc[[0]].to_numpy()*np.ones(len(y_fit))
            
            else:
                # print('oblique ridge type:' + str(r_type) + ' / ridgenb:' + k[1])
                sign = np.mean(np.diff(df[r_type][k[1]])) # slope sign
                fit_name = 'R'+ str(r_type) + ' Obl.' +  k[1]
                if sign < 0:
                    x_min = min(df[r_type][k[1]])  + 2*np.abs(max(df[0][k[1]]))
                    x_max = max(df[r_type][k[1]])
                if sign > 0:
                    x_min = min(df[r_type][k[1]])  - 2*max(df[r_type][k[1]])
                    x_max = max(df[r_type][k[1]])  #+ 2*np.abs(max(df[0][k[1]])  )
                    
                print(r_type)
                x_fit, y_fit, _ = _fit(df[r_type][k[1]],df[r_type]['elevation'],xmin=x_min, xmax=x_max) # fit function
                
               
                # plt.plot(df[r_type][k[1]].to_numpy(), m*xx + c, 'r', label='Fitted line')
            # plt.plot(x_fit, y_fit, '--')
            # plt.plot(df[i][k[1]],df[i]['depth'], 'b*')
            # print(len(x_fit))
    
            fit_xy = np.array([x_fit,y_fit]).T 
            lable_xy = np.array(['x', 'y'])
            lable = np.array([fit_name])
    
            cols = pd.MultiIndex.from_product([lable, lable_xy])   
            fit_tmp = pd.DataFrame(fit_xy, columns=cols)
            
            if k[0] == 0:
                fit_ridges_all = fit_tmp
            else:
                fit_ridges_all = pd.concat([fit_ridges_all,fit_tmp], axis=1)
    
        df_Rfit.append(fit_ridges_all) # merge ridges from different fanilly
    
    return df_Rfit


def ridges_intersection_Z0(fit, ax=None,ridge_nb=None):
    """
    Find intersection of ridges

    Parameters:

    * a
        Text here

    Returns:

    * BB : 
        return intersection by pairs

    """
    # https://stackoverflow.com/questions/28766692/intersection-of-two-graphs-in-python-find-the-x-value
    # if ax == None:
    #     fig = plt.subplots()
    #     ax = plt.gca()
    # plt.rcParams['font.size'] = 15

    # my_list = [1,2,3,4]
    # for pair in itertools.combinations(ridge_nb, r=2):
    #     print(pair)
    
    # if ridge_nb is None:
    #     ridge_nb = np.arange(0,len(fit))
        
    # for i in enumerate(ridge_nb):
        
    #     idx = np.argwhere(np.diff(np.sign(f - g))).flatten()
    #     plt.plot(x[idx], f[idx], 'ro')

    return

def scalFUN(df, EXTnb=[1], z0=0):
    """
    Analysis of ridges

    Parameters:

    * a
        Text here

    Returns:

    * BB : 
        Text here

    """
    
    
    SI = []
    
    if df.isnull().values.any():
        print('NaN or Inf detected - trying to remove')
        df.dropna(axis=1, inplace=True) # remove collumns
        df = df[~df.isin([np.nan, np.inf, -np.inf]).any(1)] #remove lines
        
        
    for i in enumerate(EXTnb):
        print(df)
        print(df['EX_xpos'+str(i[1])])
        num = np.gradient(np.log(np.abs(df['EX_xpos'+str(i[1])])))
        den = np.gradient(np.log(df['elevation']))
        # print(den)
        Tau = num/den
        q = 1./df['elevation']
        
    factor = (df['elevation'] - z0)/df['elevation']
    Tau = Tau*factor
    
    points = np.array([q,Tau]).T
    


   
    # df.isnull().values.any()
    
    x_fit, f, SI  = _fit(q,Tau,xmin=0)
    fit = np.array([x_fit, f]).T
    
    SI.append(SI_tmp)
    
    return  points, fit, SI

def scalEULER(df, EXTnb=[1], z0=0):
    """
    Analysis (Euler deconvolution of ridges)

    Parameters:

    * a
        Text here

    Returns:

    * BB : 
        Text here

    """
    

    return  points, fit #, SI


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

    label_prop = 'upwc_q' + str(qorder)
    # mesh.addprop('upwc', np.array(upw_f))
    mesh.addprop(label_prop, np.array(upw_f_dq))
    return mesh, label_prop

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

    label_prop = 'dexp_q' + str(qorder)
    mesh.addprop(label_prop, np.array(csd))
    return mesh, label_prop

# def auto_dEXP():

def _fit(x,y,**kwargs):
    """
    Curve least square fit.

    Parameters:

    * 

    Returns:

    * Intersect y(0) (with the y-axis, useful for scalFUN function)

    """
    def f(x, A, B): # this is your 'straight line' y=f(x)
        return A*x + B
   
    if np.count_nonzero(~np.isnan(x)) <2:
        raise ValueError("Need at least 3 points to fit the data")

    else:
        try:
            popt, pcov = curve_fit(f,x,y) # your data x, y to fit
            x_min = min(x) 
            x_max = max(x)                                #min/max values for x axis
        
            for key, value in kwargs.items():
                if key == 'xmin':
                   x_min = value 
                if key == 'xmax':
                   x_max = value 
            x_fit = np.linspace(x_min, x_max, 100)   #range of x values used for the fit function
            y_fit = f(x_fit, *popt)
            
            # evaluate function at y(0)
            intersect = y_fit[0]
        
        except RuntimeError:
            print("Can't fit this ridge - go to the next")
            x_fit = []   
            y_fit = []
            
        
    return x_fit, y_fit, intersect

def _build_ridge(RI_minmax,RII_minmax,RIII_minmax):
    # Once extreme points are determined at different altitudes, 
    # ridges can be obtained by linking each of them, at a given altitude, 
    # to the nearest one computed at the altitude just above.
    fit_ridges()
    
    return

def _ridges_2_df(RI_minmax, RII_minmax, RIII_minmax, **kwargs):
    
    # kwargs
    # prefix
    
    dfI = pd.DataFrame(RI_minmax)
    # df[0] = ['layer']
    dfI = dfI.add_prefix('EX_xpos')
    dfI = dfI.rename(columns={"EX_xpos0": "elevation"})
    dfI.head(5)
    
    dfII = pd.DataFrame(RII_minmax)
    # df[0] = ['layer']
    dfII = dfII.add_prefix('EX_xpos')
    dfII = dfII.rename(columns={"EX_xpos0": "elevation"})
    dfII.head(5)
    
    dfIII = pd.DataFrame(RIII_minmax)
    # df[0] = ['layer']
    dfIII = dfIII.add_prefix('EX_xpos')
    dfIII = dfIII.rename(columns={"EX_xpos0": "elevation"})
    dfIII.head(5)
    
    return dfI,dfII, dfIII
