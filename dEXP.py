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
import matplotlib.pylab as pl

from fatiando.gravmag import imaging, transform
from fatiando.gravmag.imaging import _makemesh
from fatiando import gridder, mesher, utils


from scipy.interpolate import sproot, CubicSpline, UnivariateSpline, InterpolatedUnivariateSpline, BSpline, LSQUnivariateSpline
from scipy.signal import find_peaks, peak_prominences, peak_widths # useful for ridges detection
from scipy.signal import savgol_filter
from scipy.signal import butter,filtfilt
from scipy.ndimage import gaussian_filter

import pandas as pd
from scipy.optimize import curve_fit

# class DEXP():
    
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


def ridges_minmax_plot(x, y, mesh, p1, p2, qorder=0, z=0, label='upwc',interp=True,smooth=False, **kwargs):
    plt.figure()
    
    method_peak = 'find_peaks'
    # --------------------------------------------
    # parameters to parse into find_peaks function
    for key, value in kwargs.items():
        if key == 'fix_peak_nb':
           fix_peak_nb = value
        if key == 'method_peak':
             method_peak = value  
    # prom = 0.1 #
    # fix_nb_peaks = 3
    # --------------------------------------------
    
    # This way, if z is not an array, it is now
    z = z * np.ones_like(x)
       
    RI_minmax = [] # minmax of the first horizontal derivative of the potential field   
    depths = mesh.get_zs()[:-1]
    
    if label not in mesh.props:
        raise ValueError("mesh doesn't have a '%s' property." % (label))
    else:
        upw_u = mesh.props[label]
    upw_u = np.reshape(upw_u, [mesh.shape[0], mesh.shape[1]*mesh.shape[2]])      

    for key, value in kwargs.items():
        if key == 'reverse':
           reverse = True
           print('fix_peak_nb =' + str(value))  
           print('to test')
           depths = np.reverse(depths)
           upw_u = np.reverse(upw_u)

    for i, depth in enumerate(depths - z[0]): # Loop for RII extremas
        upw_u_l = upw_u[i,:]     # analysing extrema layers by layers from top to bottom  
        # 1st horizontal derivate of the continued field
        up_f_d1x = transform.derivx(x, y, upw_u_l,(mesh.shape[1],mesh.shape[1]),order=1)
        
        if interp == True:
            xx, yy, distance, p_up_f_d1x = gridder.profile(x, y, up_f_d1x, p1, p2, 1000)
        else:
            xx, yy, distance, p_up_f_d1x, p_up_f_d1x_smooth  = profile_noInter(x, y, up_f_d1x, p1, p2, 1000)

        if smooth == True:
            p_up_f_d1x = _smooth_lowpass(xx,p_up_f_d1x)
            
        # peak analysis
        MinMax_peaks = _peaks_analysis(xx,p_up_f_d1x,fix_peak_nb=fix_peak_nb,
                                       method_peak=method_peak,proxy=None)     
        
        colors = pl.cm.viridis(np.linspace(0,1,len(depths)))
        plt.plot(xx,p_up_f_d1x, color=colors[i], label=str(int(depth)))
        for ind in range(len(MinMax_peaks)):
            plt.scatter(xx[MinMax_peaks[ind]],p_up_f_d1x[MinMax_peaks[ind]],color= 'r')
            # plt.scatter(MinMax_peaks[ind],0,color= 'r')
        plt.legend()
            

    

    
def ridges_minmax(x, y, mesh, p1, p2, qorder=0, z=0, label='upwc',interp=True,smooth=True, **kwargs):
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
        
        reverse: start peak analysis from bottom to top

    Returns:

    * BB : 
        Panda dataframe containing RI, RII and RII

    """
    # --------------------------------------------
    # parameters to parse into find_peaks function
    for key, value in kwargs.items():
        if key == 'fix_peak_nb':
           fix_peak_nb = value
                  
    # prom = 0.1 #
    # fix_nb_peaks = 3
    # --------------------------------------------
    
    # This way, if z is not an array, it is now
    z = z * np.ones_like(x)
       
    RI_minmax = [] # minmax of the first horizontal derivative of the potential field
    RII_minmax = [] # minmax of the first vertical derivative of the potential field
    RIII_minmax = [] # minmax of the potential field
   

    # --------------------------------------------
    # select depths
    depths = mesh.get_zs()[:-1]

    for key, value in kwargs.items():
        if key == 'minAlt_ridge':
            minAlt_ridge = value 
            depths = depths[depths>minAlt_ridge]    
        if key == 'maxAlt_ridge':
             maxAlt_ridge = value    
             depths = depths[depths<maxAlt_ridge]
        if key == 'method_peak':
             method_peak = value                                 

    # if minAlt_ridge is not None:
    #     print('test')
    #     # depths = 
           
    # --------------------------------------------
    
    
    
    
    if label not in mesh.props:
        raise ValueError("mesh doesn't have a '%s' property." % (label))
    else:
        upw_u = mesh.props[label]
    upw_u = np.reshape(upw_u, [mesh.shape[0], mesh.shape[1]*mesh.shape[2]])      

        
    for i, depth in enumerate(depths - z[0]): # Loop over RIII extremas
        upw_u_l = upw_u[i,:]     # analysing extrema layers by layers from top to bottom     
        
        if interp == True:
            xx, yy, distance, p_up_f = gridder.profile(x, y, upw_u_l, p1, p2, 1000)
        else:
            xx, yy, distance, p_up_f, p_up_f_smooth = profile_noInter(x, y, upw_u_l, p1, p2, 1000)

        if smooth == True:
            p_up_f = _smooth_lowpass(xx,p_up_f)

        # peak analysis
        MinMax_peaks = _peaks_analysis(xx, p_up_f,fix_peak_nb=fix_peak_nb,
                                       method_peak=method_peak,proxy=None)
        if np.array(MinMax_peaks).any():
            # RIII_minmax.append(np.hstack([[depth], MinMax_peaks]))
            RIII_minmax.append(np.hstack([[depth], xx[MinMax_peaks]]))
        else:
            RIII_minmax.append(np.hstack([[depth],[]]))
        
        if i == 3:

            plt.figure()
            plt.subplot(3,1,1)
            plt.plot(xx,p_up_f,label='u')
            for ind in range(len(MinMax_peaks)):
                plt.scatter(xx[MinMax_peaks[ind]],p_up_f[MinMax_peaks[ind]],color='g')
                # plt.scatter([MinMax_peaks[ind]],0,color='g')
            plt.legend()
 
    for i, depth in enumerate(depths - z[0]): # Loop over RII extremas
        upw_u_l = upw_u[i,:]    # analysing extrema layers by layers from top to bottom  

        # 1st vertical derivate of the continued field
        up_f_d1z = transform.derivz(x, y, upw_u_l,(mesh.shape[1],mesh.shape[1]),order=1)
        
        if interp == True:
            xx, yy, distance, p_up_f_d1z = gridder.profile(x, y, up_f_d1z, p1, p2, 1000)
        else:
            xx, yy, distance, p_up_f_d1z, p_up_f_d1z_smooth = profile_noInter(x, y, up_f_d1z, p1, p2, 1000)
            p_up_f_d1z = p_up_f_d1z_smooth

        if smooth == True:
            p_up_f_d1z = _smooth_lowpass(xx,p_up_f_d1z)
            
        # peak analysis
        MinMax_peaks = _peaks_analysis(xx,p_up_f_d1z,fix_peak_nb=fix_peak_nb,
                                       method_peak=method_peak,proxy=None)

        if np.array(MinMax_peaks).any():
            # RII_minmax.append(np.hstack([[depth], MinMax_peaks]))
            RII_minmax.append(np.hstack([[depth], xx[MinMax_peaks]]))
        else:
            RII_minmax.append(np.hstack([[depth],[]]))

        if i == 3:
            plt.subplot(3,1,2)
            plt.plot(xx,p_up_f_d1z,label='dz')
            for ind in range(len(MinMax_peaks)):
                plt.scatter(xx[MinMax_peaks[ind]],p_up_f_d1z[MinMax_peaks[ind]],color='b')
                # plt.scatter([MinMax_peaks[ind]],0,color='g')
            plt.legend()

    for i, depth in enumerate(depths - z[0]): # Loop for RII extremas
        upw_u_l = upw_u[i,:]     # analysing extrema layers by layers from top to bottom  
        # 1st horizontal derivate of the continued field
        up_f_d1x = transform.derivx(x, y, upw_u_l,(mesh.shape[1],mesh.shape[1]),order=1)
        
        if interp == True:
            xx, yy, distance, p_up_f_d1x = gridder.profile(x, y, up_f_d1x, p1, p2, 1000)
        else:
            xx, yy, distance, p_up_f_d1x, p_up_f_d1x_smooth = profile_noInter(x, y, up_f_d1x, p1, p2, 1000)
            p_up_f_d1x = p_up_f_d1x_smooth

        if smooth == True:
            p_up_f_d1x = _smooth_lowpass(xx,p_up_f_d1x)
            
        # peak analysis
        MinMax_peaks = _peaks_analysis(xx,p_up_f_d1x,fix_peak_nb=fix_peak_nb,
                                       method_peak=method_peak,proxy=None)
        
        if np.array(MinMax_peaks).any():
            # RI_minmax.append(np.hstack([[depth], MinMax_peaks]))
            RI_minmax.append(np.hstack([[depth], xx[MinMax_peaks]]))
        else:
            RI_minmax.append(np.hstack([[depth],[]]))

        if i == 3:
            plt.subplot(3,1,3)
            plt.plot(xx,p_up_f_d1x,label='dx')
            for ind in range(len(MinMax_peaks)):
                plt.scatter(xx[MinMax_peaks[ind]],p_up_f_d1x[MinMax_peaks[ind]],color='r')
                # plt.scatter([MinMax_peaks[ind]],0,color='g')
            plt.legend()
            
    # R = [np.array(RI_minmax), np.array(RII_minmax), np.array(RIII_minmax)]
    dfI,dfII, dfIII = _ridges_2_df(RI_minmax, RII_minmax, RIII_minmax)
    # R_fit = _build_ridge(RI_minmax,RII_minmax,RIII_minmax)
    
    return dfI,dfII, dfIII #, R, R_fit

def _peaks_analysis(x_axis, p_up_f, fix_peak_nb=None, 
                    method_peak='spline_roots', **kwargs):

    for key, value in kwargs.items():
        if key == 'proxy':
           proxy = value 
           if value == 'width':
               pxy = 2
           else:
               pxy = 1

    if method_peak == 'spline_roots':
        MinMax_peaks  = _spline_roots(x_axis,p_up_f)
    
    else: 
            
        if method_peak == 'find_peaks':
            Max_peaks, _ = find_peaks(p_up_f)
            
            # proxies to evaluate the peak 
            prominences = peak_prominences(p_up_f, Max_peaks)[0]
            results_half = peak_widths(p_up_f, Max_peaks, rel_height=0.5)[0]
            p_max = np.array([Max_peaks, prominences, results_half]).T
            
            if p_max.shape[0]>2:
                p_max = p_max[p_max[:,pxy].argsort()[::-1]]
        
            # --- repeat for the min --------
            Min_peaks, _ = find_peaks(-p_up_f)
            # proxies to evaluate the peak 
            prominences = peak_prominences(-p_up_f, Min_peaks)[0]
            results_half = peak_widths(-p_up_f, Min_peaks, rel_height=0.5)[0]
    
            p_min = np.array([Min_peaks, prominences, results_half]).T 
            
            if p_min.shape[0]>2:
                p_min = p_min[p_min[:,pxy].argsort()[::-1]]
                
            # MinMax_peaks= np.append(x_axis[Max_peaks],x_axis[Min_peaks])
            MinMax_peaks= np.append(Max_peaks,Min_peaks)
        
            if fix_peak_nb is not None:
                MinMax_peaks = _select_ridges_nb(fix_peak_nb,p_max,p_min)
            
        elif method_peak == 'peakdet':
            delta = 0.01
            Max_peaks, Min_peaks  = peakdet(p_up_f, delta)
            
            # MinMax_peaks= np.append(x_axis[Max_peaks],x_axis[Min_peaks])
            MinMax_peaks= np.append(Max_peaks,Min_peaks)
     
    return MinMax_peaks

def _select_ridges_nb(fix_peak_nb,p_max,p_min):

    # --- select a fixed number  --------
    if fix_peak_nb<p_max.shape[0]:
        Max_peaks_select =  p_max[0:fix_peak_nb,0].astype(int)
    else:
        Max_peaks_select = p_max[:,0].astype(int)

    warn_peak = 0
    if fix_peak_nb<p_min.shape[0]:
        Min_peaks_select = p_min[0:fix_peak_nb,0].astype(int)
    else:
        if p_min.shape[1]==0:
            warn_peak = 1
        else:
            Min_peaks_select = p_min[:,0].astype(int)

    if warn_peak == 0:
        MinMax_peaks= np.hstack([Min_peaks_select,Max_peaks_select])
    else:
        MinMax_peaks= np.hstack([Max_peaks_select])
        
    return MinMax_peaks
        

def _spline_roots(x_axis, p_up_f):
    # print('spline root analysis')
    
    # smooth = UnivariateSpline(x_axis, p_up_f, s=1)
    spl = CubicSpline(x_axis, p_up_f)
    # find function roots
    # Max_peaks = spl.roots()
    zeros_der = spl.derivative().roots()

    # index_Max_peaks=[]
    # for mp in Max_peaks:
    #     index_Max_peaks.append(abs(smooth(x_axis)-mp).argmin())
    index_Max_peaks_der=[]
    for mpder in zeros_der:
        index_Max_peaks_der.append(abs(x_axis-mpder).argmin())

    # MinMax_peaks= np.append(index_Max_peaks,index_Max_peaks_der)
    # MinMax_peaks = np.array(MinMax_peaks.astype(int))
    
    MinMax_peaks = np.array(zeros_der)
    
    # plt.scatter(x_axis, p_up_f)
    # plt.figure()
    # plt.plot(x_axis, spl(x_axis), 'g', lw=3, alpha=0.7)
    # plt.scatter(MinMax_peaks, np.zeros(MinMax_peaks.shape), color='red')
    # plt.scatter(x_axis[index_Max_peaks_der], p_up_f[index_Max_peaks_der], color='blue')
    # # plt.plot(Min_peaks, p_up_f[Min_peaks], "v")
    # plt.show()
    
    return MinMax_peaks
        
def _connect_ridges_lines(df,max_distances,gap_thresh):
    """
    Identify ridges in the 2D matrix. Expect that the width of
    the wavelet feature increases with increasing row number.

    Parameters
    ----------
    matr: 2-D ndarray
        Matrix in which to identify ridge lines.
    max_distances: 1-D sequence
        At each row, a ridge line is only connected
        if the relative max at row[n] is within
        `max_distances`[n] from the relative max at row[n+1].
    gap_thresh: int
        If a relative maximum is not found within `max_distances`,
        there will be a gap. A ridge line is discontinued if
        there are more than `gap_thresh` points without connecting
        a new relative maximum.

    Returns
    -------
    ridge_lines: tuple
        tuple of 2 1-D sequences. `ridge_lines`[ii][0] are the rows of the ii-th
        ridge-line, `ridge_lines`[ii][1] are the columns. Empty if none found.
        Each ridge-line will be sorted by row (increasing), but the order
        of the ridge lines is not specified

    References
    ----------
    Bioinformatics (2006) 22 (17): 2059-2065.
    doi: 10.1093/bioinformatics/btl355
    http://bioinformatics.oxfordjournals.org/content/22/17/2059.long

    Examples
    --------
    >>> data = np.random.rand(5,5)
    >>> ridge_lines = identify_ridge_lines(data, 1, 1)

    Notes:
    ------
    This function is intended to be used in conjuction with `cwt`
    as part of find_peaks_cwt.
    """
    # -----------------------------------------------------------------------#
    # check ridge position consistency - create a new collumn if necessary
    # dfI_2add = []
    # for k in enumerate(dfI.columns[1:]): # loop over ridges of the same familly
    #     # check presence of drop
    #     diff_test = np.diff(dfI[k[1]])
    #     print(diff_test)
    #     broken_ridge = np.where(abs(diff_test)> 5*abs(np.mean(np.diff(dfI[k[1]]))))[0]
    #     icol_max = len(dfI[k[1]])
    #     for b_r in enumerate(broken_ridge):
    #         dfI_new = dfI[k[1]][b_r[1]+1:-1]
    #         # dfI_new.add_prefix('EX_xpos')
    #         # dfI = dfI.drop(dfI[k[1]].index[b_r[1]+1:-1])
    #     dfI_2add.append(dfI_new)
    
    #     # build new ridge collumn
    #     dfI_new = 
    #     dfI[]
    #     dfI = pd.DataFrame(RI_minmax)
    #     dfI = dfI.add_prefix('EX_xpos')
    
    
def filter_ridges(dfI,dfII,dfIII,minDepth,maxDepth, minlength=3, rmvNaN=False, **kwargs):
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
    

    
    # -----------------------------------------------------------------------#
    # select a range of ridges within x limits
    # for key, value in kwargs.items():
    #     if key == 'xmin':
    #         minx = value 

    #         idfI_col_2rmv = []
    #         for k in enumerate(dfI.columns[1:]): # loop over ridges of the same familly
    #             if dfI[k[1]].any()< minx:
    #                 idfI_col_2rmv.append(k[0]+1)
                    
    #         dfI = dfI.drop(dfI.columns[idfI_col_2rmv], axis=1)
            
    #         idfII_col_2rmv = []
    #         for k in enumerate(dfII.columns[1:]): # loop over ridges of the same familly
    #             if dfII[k[1]].any()< minx:
    #                 idfII_col_2rmv.append(k[0]+1)
            
    #         dfII = dfII.drop(dfII.columns[idfII_col_2rmv], axis=1)
            
            
    #         idfIII_col_2rmv = []
    #         for k in enumerate(dfIII.columns[1:]): # loop over ridges of the same familly
    #             if dfIII[k[1]].any()< minx:
    #                 idfIII_col_2rmv.append(k[0]+1)
            
    #         dfIII = dfIII.drop(dfIII.columns[idfIII_col_2rmv], axis=1)
       
        # if key == 'xmax':
        #     maxx = value 
             
        #     for k in enumerate(dfI.columns[1:]): # loop over ridges of the same familly
        #         if dfI[k[1]].any()> maxx:
        #             dfI = dfI.drop(dfI.columns[k[0]+1], axis=1)
        #     for k in enumerate(dfII.columns[1:]): # loop over ridges of the same familly
        #         if dfII[k[1]].any()> maxx:
        #             dfII = dfII.drop(dfII.columns[k[0]+1], axis=1)
        #     for k in enumerate(dfIII.columns[1:]): # loop over ridges of the same familly
        #         if dfIII[k[1]].any()> maxx:
        #             dfIII = dfIII.drop(dfIII.columns[k[0]+1], axis=1)

    # -----------------------------------------------------------------------#
    # remove lines NaN (produce when a peak defined only for some elevation levels)
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
     
    # -----------------------------------------------------------------------#
    # regional filtering between two elevations. 
    # Particulary useful to remove noise for data close to the surface
    if dfI is not None:
        dfI = dfI.loc[(dfI['elevation'] > minDepth) & (dfI['elevation'] < maxDepth)]
    if dfII is not None:
        dfII = dfII.loc[(dfII['elevation'] > minDepth) & (dfII['elevation'] < maxDepth)]
    if dfI is not None:
        dfIII = dfIII.loc[(dfIII['elevation'] > minDepth) & (dfIII['elevation'] < maxDepth)]
    
    # -----------------------------------------------------------------------#
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
            # if np.count_nonzero(np.diff(df[r_type][k[1]]))<5: # check if ridge is vertical
            if abs(np.mean(np.diff(df[r_type][k[1]])))>15: # check if ridge is vertical
                print('vertical ridge type:' + str(r_type) + ' / ridgenb:' + k[1])
                fit_name = 'R'+ str(r_type) + ' Vert.' +  k[1]
                y_fit = np.linspace(-max(df[r_type]['elevation'])*2,max(df[r_type]['elevation']),100)                                       
                x_fit = df[r_type][k[1]].iloc[[0]].to_numpy()*np.ones(len(y_fit))
            
            else:
                print('oblique ridge type:' + str(r_type) + ' / ridgenb:' + k[1])
                sign = np.mean(np.diff(df[r_type][k[1]])) # slope sign
                fit_name = 'R'+ str(r_type) + ' Obl.' +  k[1]
                if sign < 0:
                    x_min = min(df[r_type][k[1]])  + 2*np.abs(max(df[r_type][k[1]]))
                    x_max = max(df[r_type][k[1]])
                if sign > 0:
                    x_min = min(df[r_type][k[1]])  - 2*max(df[r_type][k[1]])
                    x_max = max(df[r_type][k[1]])  #+ 2*np.abs(max(df[0][k[1]])  )
                    
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
    FIT = []
    PT = []
    if df.isnull().values.any():
        print('NaN or Inf detected - better to filter data first!')
        df.dropna(axis=1, inplace=True) # remove collumns
        df = df[~df.isin([np.nan, np.inf, -np.inf]).any(1)] #remove lines
        
        Tau = np.gradient(np.log(up_f_Centralridge)) / np.gradient(np.log(z_r))

    for i in enumerate(EXTnb):
        # print(df['EX_xpos'+str(i[1])])
        num = np.gradient(np.log(np.abs(df['EX_xpos'+str(i[1])])))
        den = np.gradient(np.log(df['elevation']))
        # print(den)
        Tau = num/den
        q = 1./df['elevation']
        
        factor = (df['elevation'] - z0)/df['elevation']
        Tau = Tau*factor
        
        points = np.array([q,Tau]).T  
        x_fit, f, si  = _fit(q,Tau,xmin=0)
        fit = np.array([x_fit, f]).T
        
        FIT.append(fit)
        PT.append(points)
        SI.append(si)
    
    return  np.array(PT), np.array(FIT), np.array(SI), EXTnb

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
            intersect = []
            
        
    return x_fit, y_fit, intersect


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

def pad_edges(xp,yp,U,shape,pad_type=2):
       
    padtypes = ['0', 'mean', 'edge', 'lintaper', 'reflection', 'oddreflection',
            'oddreflectiontaper']
    fig = plt.figure()
    ax = plt.gca()
    
    xs = xp.reshape(shape)
    ys = yp.reshape(shape)
    data = U.reshape(shape)
    
    padtype = padtypes[pad_type]
    padded_data, nps = gridder.pad_array(data, padtype=padtype)
    # Get coordinate vectors
    pad_x, pad_y = gridder.pad_coords([xs, ys], shape, nps)
    padshape = padded_data.shape
    ax.set_title(padtype)
    ax.pcolormesh(pad_y.reshape(padshape), pad_x.reshape(padshape),
                  padded_data, cmap='RdBu_r')
    ax.set_xlim(pad_y.min(), pad_y.max())
    ax.set_ylim(pad_x.min(), pad_x.max())
        
    shape = padded_data.shape
    U = padded_data.reshape(shape[0]*shape[1])
    xp = pad_x
    yp = pad_y
        
    return xp, yp, U, shape


def profile_noInter(x, y, v, point1, point2, size):
    # https://stackoverflow.com/questions/7878398/how-to-extract-an-arbitrary-line-of-values-from-a-numpy-array
    """
    Extract a profile between 2 points from spacial data.

    NO interpolation to calculate the data values at the profile points.

    Parameters:

    * x, y : 1D arrays
        Arrays with the x and y coordinates of the data points.
    * v : 1D array
        Array with the scalar value assigned to the data points.
    * point1, point2 : lists = [x, y]
        Lists the x, y coordinates of the 2 points between which the profile
        will be extracted.
    * size : int
        Number of points along the profile.

    Returns:

    * [xp, yp, distances, vp] : 1d arrays
        ``xp`` and ``yp`` are the x, y coordinates of the points along the
        profile. ``distances`` are the distances of the profile points from
        ``point1``. ``vp`` are the data points along the profile.

    """
    x1, y1 = point1
    x2, y2 = point2
    maxdist = np.sqrt((x1 - x2)**2 + (y1 - y2)**2)
    distances = np.linspace(0, maxdist, size)
    angle = np.arctan2(y2 - y1, x2 - x1)
    xp = x1 + distances*np.cos(angle)
    yp = y1 + distances*np.sin(angle)
    
    nodes = np.array([x,y]).T
    points_p = np.array([xp,yp]).T
    # find nearest point
    vp = []
    for p in points_p:
        ind = _closest_node(p, nodes)
        vp.append(v[ind])

    # window_size, poly_order = 101, 3
    # vp_smooth = savgol_filter(vp, window_size, poly_order)
    spl = UnivariateSpline(xp, vp, s=10)
    plt.plot(xp, spl(xp), 'g', lw=3)
    vp_smooth = np.array(spl(xp))
    # vp = interp_at(x, y, v, xp, yp, algorithm=algorithm, extrapolate=True)
    return xp, yp, distances, vp, vp_smooth

def _closest_node(node, nodes):
    nodes = np.asarray(nodes)
    dist_2 = np.sum((nodes - node)**2, axis=1)
    return np.argmin(dist_2)


# ----- peaks detections

from numpy import NaN, Inf, arange, isscalar, asarray, array

def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    
    Returns two arrays
    
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %      
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.
    
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    
    """
    maxtab = []
    mintab = []
       
    if x is None:
        x = arange(len(v))
    
    v = asarray(v)
    
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    
    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    
    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN
    
    lookformax = True
    
    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return np.array(maxtab), np.array(mintab)


def _smooth1d_old(x_axis, p_up):
    import csaps
    p_up_smooth = csaps.CubicSmoothingSpline(x_axis, p_up, smooth=0.25)
    # p_up_smooth = UnivariateSpline(x_axis, p_up, s=3)
    return np.array(p_up_smooth(x_axis))

def _smooth1d(x_axis, p_up ,window_len=11,window='hanning'):
    """smooth the data using a window with requested size."""
    
    from scipy.signal import savgol_filter
    # from scipy.interpolate import interp1d
    
    # itp = interp1d(x,y, kind='linear')
    window_size, poly_order = 201, 2
    yy_sg = savgol_filter(p_up, window_size, poly_order)
    # print(yy_sg)
    # print(yy_sg)

    return yy_sg

def butter_lowpass_filter(data, cutoff, fs, nyq, order):
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y
    
def _smooth_lowpass(x_axis, p_up):
    
    # Filter requirements.
    fs = abs(1/(x_axis[0] - x_axis[1]))      # sample rate, Hz
    cutoff = 0.025      # desired cutoff frequency of the filter, Hz ,      slightly higher than actual 1.2 Hz
    nyq = 0.5 * fs  # Nyquist Frequency
    order = 2      # sin wave can be approx represented as quadratic
        
    filtdata = butter_lowpass_filter(p_up, cutoff, fs, nyq, order)
    return filtdata

def smooth2d(x, y, U, sigma = 10):

    plt.figure()
    plt.subplot(1,2,1)
    plt.scatter(x, y, c=U, cmap='viridis',vmax=0.25)
    plt.colorbar()
    plt.axis('square')
    # plt.show()
    
    U2d = U.reshape(int(np.sqrt(U.shape)),int(np.sqrt(U.shape)))
    U2d_f = gaussian_filter(U2d, sigma=sigma)
    
    U_f = np.copy(U)
    U_f = U2d_f.reshape(U.shape)
    plt.subplot(1,2,2)
    plt.scatter(x, y, c=U_f, cmap='viridis',vmax=0.25)
    plt.colorbar()
    plt.axis('square')
    plt.show()
    
    return U_f