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


from scipy.interpolate import (
    sproot,
    CubicSpline,
    UnivariateSpline,
    InterpolatedUnivariateSpline,
    BSpline,
    LSQUnivariateSpline,
)
from scipy.signal import (
    find_peaks,
    peak_prominences,
    peak_widths,
)  # useful for ridges detection
from scipy.signal import savgol_filter
from scipy.signal import butter, filtfilt
from scipy.ndimage import gaussian_filter
from scipy import interpolate

import pandas as pd
from scipy.optimize import curve_fit

# class DEXP():


def ridges_minmax_plot(
    x,
    y,
    mesh,
    p1,
    p2,
    qorder=0,
    z=0,
    fix_peak_nb=None,
    label="upwc",
    interp=True,
    smooth=False,
    showfig=False,
    **kwargs
):
    """
    * x, y : 1D-arrays
        The x and y coordinates of the grid points
    * mesh : fatiando.mesher.PrismMesh
        The estimated physical property distribution set in a prism mesh (for easy 3D plotting)
    * p1, p2 : 1D-arrays
        The p1 and p2 coordinates of the extracted profile end points
    * qorder : int
        The derivative order
    * z : float or 1D-array
        The z coordinate of the grid points
    * fix_peak_nb : int
        The maximum number of peak to identify
    * label : string
        Label of the estimated physical property distribution
    * interp : True or False
        If True, will interpolate values between the data points.
    * smooth : True or False
        If True, will apply a low-pass filter values.
    * showfig : True or False
        If True, will display the figure.

    Returns:

    """
    if showfig == True:
        ax = plt.figure()


    method_peak = "find_peaks"
    x_resolution = len(x)
    Xaxis = "x"
    # --------------------------------------------
    # parameters to parse into find_peaks function
    for key, value in kwargs.items():
        if key == "fix_peak_nb":
            fix_peak_nb = value
        if key == "method_peak":
            method_peak = value
        if key == "x_resolution":
            x_resolution = value
        if key == "Xaxis":
            Xaxis = value

            # prom = 0.1 #
    # fix_nb_peaks = 3
    # --------------------------------------------

    # This way, if z is not an array, it is now
    z = z * np.ones_like(x)

    RI_minmax = []  # minmax of the first horizontal derivative of the potential field
    depths = mesh.get_zs()[:-1]

    if label not in mesh.props:
        raise ValueError("mesh doesn't have a '%s' property." % (label))
    else:
        upw_u = mesh.props[label]
    upw_u = np.reshape(upw_u, [mesh.shape[0], mesh.shape[1] * mesh.shape[2]])

    for key, value in kwargs.items():
        if key == "reverse":
            reverse = True
            print("fix_peak_nb =" + str(value))
            # print("to test")
            depths = np.reverse(depths)
            upw_u = np.reverse(upw_u)

    for i, depth in enumerate(depths - z[0]):  # Loop for RII extremas
        upw_u_l = upw_u[i, :]  # analysing extrema layers by layers from top to bottom
        # 1st horizontal derivate of the continued field
        # up_f_d1x = transform.derivx(x, y, upw_u_l,(mesh.shape[1],mesh.shape[1]),order=1)
        if Xaxis == "dist":
            raise ValueError(
                "Impossible to calculate the derivative in this profil direction"
            )
        elif Xaxis == "y":
            up_f_d1x = transform.derivx(
                x, y, upw_u_l, (mesh.shape[1], mesh.shape[1]), order=1
            )
        else:
            up_f_d1x = transform.derivy(
                x, y, upw_u_l, (mesh.shape[1], mesh.shape[1]), order=1
            )

        if interp == True:
            xx, yy, distance, p_up_f_d1x = gridder.profile(
                x, y, up_f_d1x, p1, p2, x_resolution
            )
        else:
            xx, yy, distance, p_up_f_d1x, p_up_f_d1x_dict = profile_noInter(
                x, y, up_f_d1x, p1, p2, x_resolution, showfig=False
            )

        if Xaxis == "dist":
            xaxis = distance
        elif Xaxis == "y":
            xaxis = xx
        else:
            xaxis = yy

        if smooth is not False:
            if smooth == True:
                p_up_f_d1x = _smooth_lowpass(xaxis, p_up_f_d1x)
            else:
                p_up_f_d1x = p_up_f_d1x_dict[smooth]

        # if smooth == True:
        #     p_up_f_d1x = _smooth_lowpass(xaxis,p_up_f_d1x)

        # peak analysis
        MinMax_peaks, MinMax_peaks_p = _peaks_analysis(
            xaxis,
            p_up_f_d1x,
            fix_peak_nb=fix_peak_nb,
            method_peak=method_peak,
            proxy=None,
        )

        if showfig == True:
            colors = pl.cm.viridis(np.linspace(0, 1, len(depths)))
            plt.plot(xaxis, p_up_f_d1x, color=colors[i], label=str(int(depth)))
            for ind in range(len(MinMax_peaks)):
                plt.scatter(
                    xaxis[MinMax_peaks[ind]], p_up_f_d1x[MinMax_peaks[ind]], color="r"
                )
                # plt.scatter(MinMax_peaks[ind],0,color= 'r')
            plt.legend()
            
    return ax


def ridges_minmax(
    x,
    y,
    mesh,
    p1,
    p2,
    qorder=0,
    z=0,
    label="upwc",
    fix_peak_nb=None,
    interp=True,
    smooth=False,
    showfig=False,
    **kwargs
):
    """
    Form a multiridge set
    RI and RII : zeros of the first horizontal and first vertical derivatives of the potential field
    RIII :zeros of the potential field itself

    .. note:: ridges generated by isolated simple sources point, line, sheet, and contact are straight lines defined by the zeros of a potential
    field and its horizontal and vertical derivatives at all measured or computed levels

    Parameters:

    * x, y : 1D-arrays
        The x and y coordinates of the grid points
    * mesh : fatiando.mesher.PrismMesh
        The estimated physical property distribution set in a prism mesh (for easy 3D plotting). The upward continuated field mesh (of order-q derivative).
        Read the property label 'upwc' of the mesh by default
    * p1, p2 : 1D-arrays
        The p1 and p2 coordinates of the extracted profile end points
    * qorder : int
        The derivative order
    * z : float or 1D-array
        The z coordinate of the grid points
    * label : string
        Label of the estimated physical property distribution
    * fix_peak_nb : int
        The maximum number of peak to identify
    * interp : True or False
        If True, will interpolate values between the data points.
    * smooth : True or False
        If True, will apply a low-pass filter values.
    * showfig : True or False
        If True, will display the figure.

    **kwargs
        prominence for peak detection
        reverse: start peak analysis from bottom to top

    Returns:

    * MinMax_peaks :
        Panda dataframe containing ridges RI, RII and RII

    """
    # --------------------------------------------
    # parameters to parse into find_peaks function
    for key, value in kwargs.items():
        if key == "fix_peak_nb":
            fix_peak_nb = value
    # --------------------------------------------

    # This way, if z is not an array, it is now
    z = z * np.ones_like(x)

    # indices of the peaks
    RI_minmax = []  # minmax of the first horizontal derivative of the potential field
    RII_minmax = []  # minmax of the first vertical derivative of the potential field
    RIII_minmax = []  # minmax of the potential field

    # dictionnary of the peaks properties
    RI_minmax_p = []
    RII_minmax_p = []
    RIII_minmax_p = []

    # --------------------------------------------
    # select depths
    depths = mesh.get_zs()[:-1]
    x_resolution = len(x)
    Xaxis = "x"
    peakp_out = False
    iplot = 4

    for key, value in kwargs.items():
        if key == "minAlt_ridge":
            minAlt_ridge = value
            depths = depths[depths > minAlt_ridge]
        if key == "maxAlt_ridge":
            maxAlt_ridge = value
            depths = depths[depths < maxAlt_ridge]
        if key == "method_peak":
            method_peak = value
        if key == "x_resolution":
            x_resolution = value
        if key == "Xaxis":
            Xaxis = value
        if key == "returnAmp":
            peakp_out = True
        if key == "iplot":
            iplot = value
    # --------------------------------------------

    if label not in mesh.props:
        raise ValueError("mesh doesn't have a '%s' property." % (label))
    else:
        upw_u = mesh.props[label]
    upw_u = np.reshape(upw_u, [mesh.shape[0], mesh.shape[1] * mesh.shape[2]])

    for i, depth in enumerate(depths - z[0]):  # Loop over RIII extremas
        upw_u_l = upw_u[i, :]  # analysing extrema layers by layers from top to bottom

        if interp == True:
            xx, yy, distance, p_up_f = gridder.profile(
                x, y, upw_u_l, p1, p2, x_resolution
            )
        else:
            xx, yy, distance, p_up_f, p_up_f_dict = profile_noInter(
                x, y, upw_u_l, p1, p2, x_resolution
            )

        if Xaxis == "dist":
            xaxis = distance
        elif Xaxis == "y":
            xaxis = xx
        else:
            xaxis = yy

        # if smooth is not None:
        #     p_up_f = _smooth_lowpass(xaxis,p_up_f)
        #     p_up_f =
        # p_up_f = np.array(p_up_f)
        if smooth is not False:
            if smooth == True:
                p_up_f_d1z = _smooth_lowpass(xaxis, p_up_f)
            else:
                p_up_f = p_up_f_dict[smooth]
        p_up_f = np.array(p_up_f)

        # print(p_up_f)

        # peak analysis
        MinMax_peaks, MinMax_peaks_p = _peaks_analysis(
            xaxis, p_up_f, fix_peak_nb=fix_peak_nb, method_peak=method_peak, proxy=None
        )
        if np.array(MinMax_peaks).any():
            # RIII_minmax.append(np.hstack([[depth], MinMax_peaks]))
            RIII_minmax.append(np.hstack([[depth], xaxis[MinMax_peaks]]))
            RIII_minmax_p.append(np.hstack([[depth], MinMax_peaks_p]))
        else:
            RIII_minmax.append(np.hstack([[depth], []]))
            RIII_minmax_p.append(np.hstack([[depth], MinMax_peaks_p]))

        if showfig == True:
            if i == iplot:
                ax = plt.figure()
                plt.subplot(3, 1, 1)
                plt.plot(xaxis, p_up_f, label="u")
                for ind in range(len(MinMax_peaks)):
                    plt.scatter(
                        xaxis[MinMax_peaks[ind]], p_up_f[MinMax_peaks[ind]], color="g"
                    )
                    # plt.scatter([MinMax_peaks[ind]],0,color='g')
                plt.legend()
                ax_list = plt.gca()
                ax_list.xaxis.set_ticklabels([])
        else:
            ax = None

    for i, depth in enumerate(depths - z[0]):  # Loop over RII extremas
        upw_u_l = upw_u[i, :]  # analysing extrema layers by layers from top to bottom

        # 1st vertical derivate of the continued field
        up_f_d1z = transform.derivz(
            x, y, upw_u_l, (mesh.shape[1], mesh.shape[1]), order=1
        )

        if interp == True:
            xx, yy, distance, p_up_f_d1z = gridder.profile(
                x, y, up_f_d1z, p1, p2, x_resolution
            )
        else:
            xx, yy, distance, p_up_f_d1z, p_up_f_d1z_dict = profile_noInter(
                x, y, up_f_d1z, p1, p2, x_resolution
            )
            # p_up_f_d1z = p_up_f_d1z_smooth

        if Xaxis == "dist":
            xaxis = distance
        elif Xaxis == "y":
            xaxis = xx
        else:
            xaxis = yy

        # if smooth == True:
        #     p_up_f_d1z = _smooth_lowpass(xaxis,p_up_f_d1z)
        # p_up_f = np.array(p_up_f)

        if smooth is not False:
            if smooth == True:
                p_up_f_d1z = _smooth_lowpass(xaxis, p_up_f_d1z)
            else:
                p_up_f_d1z = p_up_f_d1z_dict[smooth]

        # peak analysis
        MinMax_peaks, MinMax_peaks_p = _peaks_analysis(
            xaxis,
            p_up_f_d1z,
            fix_peak_nb=fix_peak_nb,
            method_peak=method_peak,
            proxy=None,
        )

        if np.array(MinMax_peaks).any():
            # RII_minmax.append(np.hstack([[depth], MinMax_peaks]))
            RII_minmax.append(np.hstack([[depth], xaxis[MinMax_peaks]]))
            RII_minmax_p.append(np.hstack([[depth], MinMax_peaks_p]))
        else:
            RII_minmax.append(np.hstack([[depth], []]))
            RII_minmax_p.append(np.hstack([[depth], []]))

        if showfig == True:
            if i == iplot:
                plt.subplot(3, 1, 2)
                plt.plot(xaxis, p_up_f_d1z, label="dz")
                for ind in range(len(MinMax_peaks)):
                    plt.scatter(
                        xaxis[MinMax_peaks[ind]],
                        p_up_f_d1z[MinMax_peaks[ind]],
                        color="b",
                    )
                    # plt.scatter([MinMax_peaks[ind]],0,color='g')
                plt.legend()
                ax_list = plt.gca()
                ax_list.xaxis.set_ticklabels([])


    for i, depth in enumerate(depths - z[0]):  # Loop for RII extremas
        upw_u_l = upw_u[i, :]  # analysing extrema layers by layers from top to bottom
        # 1st horizontal derivate of the continued field

        if Xaxis == "dist":
            raise ValueError(
                "Impossible to calculate the derivative in this profil direction"
            )
        elif Xaxis == "y":
            up_f_d1x = transform.derivx(
                x, y, upw_u_l, (mesh.shape[1], mesh.shape[1]), order=1
            )
        else:
            up_f_d1x = transform.derivy(
                x, y, upw_u_l, (mesh.shape[1], mesh.shape[1]), order=1
            )

        if interp == True:
            xx, yy, distance, p_up_f_d1x = gridder.profile(
                x, y, up_f_d1x, p1, p2, x_resolution
            )
        else:
            xx, yy, distance, p_up_f_d1x, p_up_f_d1x_dict = profile_noInter(
                x, y, up_f_d1x, p1, p2, x_resolution
            )
            # p_up_f_d1x = p_up_f_d1x_smooth

        if Xaxis == "dist":
            xaxis = distance
        elif Xaxis == "y":
            xaxis = xx
        else:
            xaxis = yy

        if smooth is not False:
            if smooth == True:
                p_up_f_d1x = _smooth_lowpass(xaxis, p_up_f_d1x)
            else:
                p_up_f_d1x = p_up_f_d1x_dict[smooth]

        # peak analysis
        MinMax_peaks, MinMax_peaks_p = _peaks_analysis(
            xaxis,
            p_up_f_d1x,
            fix_peak_nb=fix_peak_nb,
            method_peak=method_peak,
            proxy=None,
        )

        if np.array(MinMax_peaks).any():
            RI_minmax.append(np.hstack([[depth], xaxis[MinMax_peaks]]))
            RI_minmax_p.append(np.hstack([[depth], MinMax_peaks_p]))
        else:
            RI_minmax.append(np.hstack([[depth], []]))
            RI_minmax_p.append(np.hstack([[depth], []]))

        if showfig == True:
            if i == iplot:
                plt.subplot(3, 1, 3)
                plt.plot(xaxis, p_up_f_d1x, label="dx")
                for ind in range(len(MinMax_peaks)):
                    plt.scatter(
                        xaxis[MinMax_peaks[ind]],
                        p_up_f_d1x[MinMax_peaks[ind]],
                        color="r",
                    )
                    # plt.scatter([MinMax_peaks[ind]],0,color='g')
                plt.legend()
                plt.tight_layout()

    # R = [np.array(RI_minmax), np.array(RII_minmax), np.array(RIII_minmax)]
    dfI, dfII, dfIII = _ridges_2_df(RI_minmax, RII_minmax, RIII_minmax)
    dfIp, dfIIp, dfIIIp = _ridges_2_df(RI_minmax_p, RII_minmax_p, RIII_minmax_p)
    # R_fit = _build_ridge(RI_minmax,RII_minmax,RIII_minmax)

    if peakp_out == False:
        return dfI, dfII, dfIII, ax
    else:
        return dfI, dfII, dfIII, dfIp, dfIIp, dfIIIp, ax # , R, R_fit


def _peaks_analysis(
    x_axis, p_up_f, fix_peak_nb=None, method_peak="spline_roots", **kwargs
):
    """
    search for peaks in the signal using a given algorithm defined by method_peak-

    .. note:: fix_peak_nb constrainst help to differentiate between main peaks and secondary peaks often required to exclude

    Parameters:

    * x_axis : string
        The direction of the 2d profile
    * p_up_f : 1D-arrays
        Signal to analyse
    * fix_peak_nb : int
        Constrainst on number of peaks to detect
    * method_peak : string
        Constrainst on algoritm to use to detect the peaks

    **kwargs
        * showfig : True or False
            If True, will display the figure.
    Returns:

    * MinMax_peaks :
        numpy array containing min and max of peaks
    * p_up_f[MinMax_peaks]
        numpy array containing amplitude associated with the min and max of peaks

    """
    for key, value in kwargs.items():
        if key == "proxy":
            proxy = value
            if value == "width":
                pxy = 2
            else:
                pxy = 1

    if method_peak == "spline_roots":
        MinMax_peaks = _spline_roots(x_axis, p_up_f)

    else:

        if method_peak == "find_peaks":
            Max_peaks, _ = find_peaks(p_up_f, height=None)
            prominences = peak_prominences(p_up_f, Max_peaks)[0]
            results_half = peak_widths(p_up_f, Max_peaks, rel_height=0.5)[0]
            p_max = np.array([Max_peaks, prominences, results_half]).T

            if p_max.shape[0] > 2:
                p_max = p_max[p_max[:, pxy].argsort()[::-1]]

            # --- repeat for the min --------
            Min_peaks, _ = find_peaks(-p_up_f)
            # proxies to evaluate the peak
            prominences = peak_prominences(-p_up_f, Min_peaks)[0]
            results_half = peak_widths(-p_up_f, Min_peaks, rel_height=0.5)[0]

            p_min = np.array([Min_peaks, prominences, results_half]).T

            if p_min.shape[0] > 2:
                p_min = p_min[p_min[:, pxy].argsort()[::-1]]

            # MinMax_peaks= np.append(x_axis[Max_peaks],x_axis[Min_peaks])
            MinMax_peaks = np.append(Max_peaks, Min_peaks)

            if fix_peak_nb is not None:
                MinMax_peaks = _select_ridges_nb(fix_peak_nb, p_max, p_min)
            # MinMax_peaks_p
        elif method_peak == "peakdet":
            delta = 0.01
            Max_peaks, Min_peaks = peakdet(p_up_f, delta)

            # MinMax_peaks= np.append(x_axis[Max_peaks],x_axis[Min_peaks])
            MinMax_peaks = np.append(Max_peaks, Min_peaks)

    return MinMax_peaks, p_up_f[MinMax_peaks]


def _select_ridges_nb(fix_peak_nb, p_max, p_min):

    # print('_select_ridges_nb')
    # --- select a fixed number  --------
    if fix_peak_nb < p_max.shape[0]:
        Max_peaks_select = p_max[0:fix_peak_nb, 0].astype(int)
    else:
        Max_peaks_select = p_max[:, 0].astype(int)

    warn_peak = 0
    if fix_peak_nb < p_min.shape[0]:
        Min_peaks_select = p_min[0:fix_peak_nb, 0].astype(int)
    else:
        if p_min.shape[1] == 0:
            warn_peak = 1
        else:
            Min_peaks_select = p_min[:, 0].astype(int)

    if warn_peak == 0:
        MinMax_peaks = np.hstack([Min_peaks_select, Max_peaks_select])
    else:
        MinMax_peaks = np.hstack([Max_peaks_select])

    return MinMax_peaks


def _spline_roots(x_axis, p_up_f):

    spl = CubicSpline(x_axis, p_up_f)
    zeros_der = spl.derivative().roots()

    index_Max_peaks_der = []
    for mpder in zeros_der:
        index_Max_peaks_der.append(abs(x_axis - mpder).argmin())

    MinMax_peaks = np.array(zeros_der)

    return MinMax_peaks


def _connect_ridges_lines(df, max_distances, gap_thresh):
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


def filter_ridges(
    dfIf, dfIIf, dfIIIf, minDepth, maxDepth, minlength=3, rmvNaN=False, **kwargs
):
    """
    Filter non-rectiligne ridges (denoising)

    Parameters:
    dfI: dataframe
        contains ridges of type I
    dfII: dataframe
        contains ridges of type II
    dfIII: dataframe
        contains ridges of type III
    * minDepth
        Text here
    * maxDepth
    * kwargs
    heights: dataframe
        contains heights of ridges of type I


    Returns:

    * BB :
        Text here

    """


    height1f, height2f, height3f = [None, None, None]
    for key, value in kwargs.items():
        if key == "heights":
            height1f, height2f, height3f = value
            # print(height2)
    # -----------------------------------------------------------------------#
    # select a range of ridges within x limits
    for key, value in kwargs.items():
        if key == "xmin":
            minx = value

            for k in enumerate(
                dfIf.columns[1:]
            ):  # loop over ridges of the same familly
                id_2NaN = np.where(dfIf[k[1]] < minx)
                dfIf[k[1]].iloc[id_2NaN] = np.nan
                if height1f is not None:
                    height1f[k[1]].iloc[id_2NaN] = np.nan

            for k in enumerate(
                dfIIf.columns[1:]
            ):  # loop over ridges of the same familly
                id_2NaN = np.where(dfIIf[k[1]] < minx)
                dfIIf[k[1]].iloc[id_2NaN] = np.nan
                if height2f is not None:
                    height2f[k[1]].iloc[id_2NaN] = np.nan

            for k in enumerate(
                dfIIIf.columns[1:]
            ):  # loop over ridges of the same familly
                id_2NaN = np.where(dfIIIf[k[1]] < minx)
                dfIIIf[k[1]].iloc[id_2NaN] = np.nan
                if height3f is not None:
                    height3f[k[1]].iloc[id_2NaN] = np.nan

        if key == "xmax":
            maxx = value

            for k in enumerate(
                dfIf.columns[1:]
            ):  # loop over ridges of the same familly
                id_2NaN = np.where(dfIf[k[1]] > maxx)
                dfIf[k[1]].iloc[id_2NaN] = np.nan
                if height1f is not None:
                    height1f[k[1]].iloc[id_2NaN] = np.nan

            for k in enumerate(
                dfIIf.columns[1:]
            ):  # loop over ridges of the same familly
                id_2NaN = np.where(dfIIf[k[1]] > maxx)
                dfIIf[k[1]].iloc[id_2NaN] = np.nan
                if height2f is not None:
                    height2f[k[1]].iloc[id_2NaN] = np.nan

            for k in enumerate(
                dfIIIf.columns[1:]
            ):  # loop over ridges of the same familly
                id_2NaN = np.where(dfIIIf[k[1]] > maxx)
                dfIIIf[k[1]].iloc[id_2NaN] = np.nan
                if height3f is not None:
                    height3f[k[1]].iloc[id_2NaN] = np.nan

    # -----------------------------------------------------------------------#
    # remove lines NaN (produce when a peak defined only for some elevation levels)
    if rmvNaN == True:
        if dfIf.isnull().values.any():
            print("NaN or Inf detected - trying to remove")
            dfIf.dropna(axis=1, inplace=True)  # remove collumns
            dfIf = dfIf[~dfIf.isin([np.nan, np.inf, -np.inf]).any(1)]  # remove lines

            if height1f is not None:
                height1f.dropna(axis=1, inplace=True)  # remove collumns
                height1f = height1f[
                    ~height1f.isin([np.nan, np.inf, -np.inf]).any(1)
                ]  # remove lines

        if dfIIf.isnull().values.any():
            dfIIf.dropna(axis=1, inplace=True)  # remove collumns
            dfIIf = dfIIf[~dfIIf.isin([np.nan, np.inf, -np.inf]).any(1)]  # remove lines

            if height2f is not None:
                # print(height2)
                height2f.dropna(axis=1, inplace=True)  # remove collumns
                height2f = height2f[
                    ~height2f.isin([np.nan, np.inf, -np.inf]).any(1)
                ]  # remove lines
                # print(height2)

        if dfIIIf.isnull().values.any():
            dfIIIf.dropna(axis=1, inplace=True)  # remove collumns
            dfIIIf = dfIIIf[
                ~dfIIIf.isin([np.nan, np.inf, -np.inf]).any(1)
            ]  # remove lines

            if height3f is not None:
                height3f.dropna(axis=1, inplace=True)  # remove collumns
                height3f = height3f[
                    ~height3f.isin([np.nan, np.inf, -np.inf]).any(1)
                ]  # remove lines

    # -----------------------------------------------------------------------#
    # regional filtering between two elevations.
    # Particulary useful to remove noise for data close to the surface
    if dfIf is not None:
        dfIf = dfIf.loc[(dfIf["elevation"] > minDepth) & (dfIf["elevation"] < maxDepth)]
        if height1f is not None:
            height1f = height1f.loc[
                (height1f["elevation"] > minDepth) & (height1f["elevation"] < maxDepth)
            ]

    if dfIIf is not None:
        dfIIf = dfIIf.loc[
            (dfIIf["elevation"] > minDepth) & (dfIIf["elevation"] < maxDepth)
        ]
        if height2f is not None:
            height2f = height2f.loc[
                (height2f["elevation"] > minDepth) & (height2f["elevation"] < maxDepth)
            ]

    if dfIIIf is not None:
        dfIIIf = dfIIIf.loc[
            (dfIIIf["elevation"] > minDepth) & (dfIIIf["elevation"] < maxDepth)
        ]
        if height3f is not None:
            height3f = height3f.loc[
                (height3f["elevation"] > minDepth) & (height3f["elevation"] < maxDepth)
            ]

    # -----------------------------------------------------------------------#
    # check length of ridges (remove column if less than N points)
    if dfIf is not None:
        smallCol = dfIf.count()
        idrmv = np.where(smallCol < minlength)[0].tolist()
        dfIf = dfIf.drop(dfIf.columns[idrmv], axis=1)
        if height1f is not None:
            height1f = height1f.drop(height1f.columns[idrmv], axis=1)

    if dfIIf is not None:
        smallCol = dfIIf.count()
        idrmv = np.where(smallCol < minlength)[0].tolist()
        dfIIf = dfIIf.drop(dfIIf.columns[idrmv], axis=1)
        if height2f is not None:
            height2f = height2f.drop(height2f.columns[idrmv], axis=1)

    if dfIIIf is not None:
        smallCol = dfIIIf.count()
        idrmv = np.where(smallCol < minlength)[0].tolist()
        dfIIIf = dfIIIf.drop(dfIIIf.columns[idrmv], axis=1)
        if height3f is not None:
            height3f = height3f.drop(height3f.columns[idrmv], axis=1)

    if height1f is not None:
        return dfIf, dfIIf, dfIIIf, height1f, height2f, height3f
    else:
        return dfIf, dfIIf, dfIIIf


def fit_ridges(df, rmvOutliers=False):
    """
    Fit ridges and return points and fit equations to plot

    Parameters:

    * df
        dataframe including all tree types of ridges

    Returns:

    * BB :
        points and fit equations to plot

    """
    if len(df) == 1:
        df = [df]

    if len(df[0]) == 0:
        raise ValueError("No data to fit")

    df_Rfit = []
    for r_type in range(len(df)):  # loop over ridges type I, II, III
        fit_ridges_all = []
        lable = []
        cols = []

        for k in enumerate(
            df[r_type].columns[1:]
        ):  # loop over ridges of the same familly
            # if abs(np.mean(np.diff(df[r_type][k[1]])))>1: # check if ridge is vertical
            # print(df[r_type].columns[1:])
            sign = np.mean(np.diff(df[r_type][k[1]]))  # slope sign
            slope =  abs(min(df[r_type][k[1]])-max(df[r_type][k[1]]))/abs(min(df[r_type]["elevation"])-max(df[r_type]["elevation"]))
            # print(sign)
            if sign == 0:  # check if ridge is vertical
                # print( abs(np.mean(np.diff(df[r_type][k[1]]))))
                # if abs(np.mean(np.diff(df[r_type][k[1]])))>1: # check if ridge is vertical
                # print("vertical ridge type:" + str(r_type) + " / ridgenb:" + k[1])
                fit_name = "R" + str(r_type) + " Vert." + k[1]
                y_fit = np.linspace(
                    -max(df[r_type]["elevation"]) * 2, max(df[r_type]["elevation"]), 100
                )
                x_fit = df[r_type][k[1]].iloc[[0]].to_numpy() * np.ones(len(y_fit))
                # print(x_fit,y_fit)
            else:
                # print("oblique ridge type:" + str(r_type) + " / ridgenb:" + k[1])
                fit_name = "R" + str(r_type) + " Obl." + k[1]



                x_values = df[r_type][k[1]]

                # print(x_values)
                # print(x_max)
                # print(x_min)
                y_values = df[r_type]["elevation"]
                # print(x_values)
                # objective function
                def objective(x, a, b, c):
                	return a * x + b

                # fit curve
                popt, _ = curve_fit(objective, x_values, y_values)

                # define new input values
                x_new=x_values
                # unpack optima parameters for the objective function
                a, b, c = popt
                # use optimal parameters to calculate new values
                y_new = objective(x_new, a, b, c)

                # x_fit, y_fit, _ = _fit(
                #     df[r_type][k[1]],
                #     df[r_type]["elevation"],
                #     slope=sign,
                #     rmvOutliers=rmvOutliers,
                # )  # fit function

                x_fit, y_fit = x_new, y_new
                # print(a)

                if a < 0:
                    x_min = min(df[r_type][k[1]])
                    x_max = max(df[r_type][k[1]])  + abs(a)*10
                if a > 0:
                    x_min = min(df[r_type][k[1]])
                    x_max = max(df[r_type][k[1]]) - a*10


                x_new = np.linspace(x_min,x_max,100)
                # # use optimal parameters to calculate new values
                y_new = objective(x_new, a, b, c)

                x_fit, y_fit = x_new, y_new

                # x_fit = df[r_type][k[1]]
                # y_fit = slope*x_fit #- min(df[r_type][k[1]])
                # plt.plot(df[r_type][k[1]].to_numpy(), m*xx + c, 'r', label='Fitted line')
            # plt.figure()
            # plt.title(r_type)
            # plt.plot(x_fit, y_fit, '--')
            # plt.plot(df[r_type][k[1]],df[r_type]['elevation'], 'b*')
            # print(len(x_fit))

            fit_xy = np.array([x_fit, y_fit]).T
            lable_xy = np.array(["x", "y"])
            lable = np.array([fit_name])

            cols = pd.MultiIndex.from_product([lable, lable_xy])
            fit_tmp = pd.DataFrame(fit_xy, columns=cols)

            if k[0] == 0:
                fit_ridges_all = fit_tmp
            else:
                fit_ridges_all = pd.concat([fit_ridges_all, fit_tmp], axis=1)

        df_Rfit.append(fit_ridges_all)  # merge ridges from different fanilly

    return df_Rfit


def ridges_intersection_Z0(fit, ax=None, ridge_nb=None):
    """
    Find intersection of ridges (NOT YET IMPLEMENTED)

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


def scalFUN(df, EXTnb=[1], z0=0, **kwargs):
    """
    Analysis of ridges (NOT YET IMPLEMENTED)

    Parameters:

    * a
        Text here

    Returns:

    * BB :
        Text here

    """
    rmvOutliers = False
    for key, value in kwargs.items():
        rmvOutliers = value

    SI = []
    FIT = []
    PT = []
    if df.isnull().values.any():
        print("NaN or Inf detected - better to filter data first!")
        df.dropna(axis=1, inplace=True)  # remove collumns
        df = df[~df.isin([np.nan, np.inf, -np.inf]).any(1)]  # remove lines

        # Tau = np.gradient(np.log(up_f_Centralridge)) / np.gradient(np.log(z_r))

    # for i in enumerate(EXTnb):
    # for i in EXTnb:
    i = EXTnb

    if rmvOutliers == True:
        # _ , filtered_entries, quartileSet = _removeOutliers(df['EX_xpos'+str(i[1])])
        _, filtered_entries, quartileSet = _removeOutliers(df[i])
        # print("------")
        df = df.iloc[filtered_entries, :]

    # print(df['EX_xpos'+str(i[1])])
    # num = np.gradient(np.log(np.abs(df['EX_xpos'+str(i[1])])))

    # plt.figure()
    # plt.plot(np.log(np.abs(df[i])))

    num = np.gradient(np.log(np.abs(df[i])))
    den = np.gradient(np.log(df["elevation"]))
    # print('num')
    # print(num)

    # plt.figure()
    # # plt.scatter(np.log(df['elevation']),np.log(np.abs(df[i])))
    # plt.scatter(np.log(1/df['elevation']),np.log(np.abs(df[i])))
    # plt.xlabel('log (1/z)')
    # plt.ylabel('log (ridge amplitude)')

    # print(den)
    Tau = num / den
    # print(df[i])

    # dzz = np.log(df[1]['elevation'].iloc[1])-np.log(df[0]['elevation'].iloc[0])
    Tau2 = np.gradient(np.log(np.abs(df[i])), np.log(df["elevation"]))
    # Tau2 = np.gradient(np.log(np.abs(df[i])),4.2)
    # print('Tau')
    # print(Tau)
    # print('Tau2')
    # print(Tau2)

    q = 1.0 / df["elevation"]

    factor = (df["elevation"] - z0) / df["elevation"]
    # factor = df['elevation']/(df['elevation'] - z0)
    # factor = 1/(z0)
    # factor = 1
    Tau = Tau2 * factor

    points = np.array([q, Tau]).T
    x_fit, f, si = _fit(q, Tau, xmin=0, SI=True)
    fit = np.array([x_fit, f]).T

    FIT.append(fit)
    PT.append(points)
    SI.append(si)

    return np.array(PT), np.array(FIT), np.array(SI), EXTnb


def scalEULER(df, EXTnb=[1], z0=0):
    """
    Analysis (Euler deconvolution of ridges, Fedi 2019) (NOT YET IMPLEMENTED)

    Parameters:

    * a
        Text here

    Returns:

    * BB :
        Text here

    """

    return points, fit  # , SI


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

    if zmin==0: # avoid downward continuation
        zmin =+ 1e-3
    # This way, if z is not an array, it is now
    z = z * np.ones_like(x)
    # Remove the last z because I only want depths to the top of the layers
    depths = mesh.get_zs()[:-1]

    upw_f = []
    upw_f_dq = []
    # Offset by the data z because in the paper the data is at z=0
    for depth in depths - z[0]:

        # continued field calculation
        upw_fhi = transform.upcontinue(x, y, data, shape, depth)

        # qorder vertical derivate of the continued field
        upw_f_dqhi = transform.derivz(x, y, upw_fhi, shape, order=qorder)

        # print(np.mean(upw_fhi))
        # print(np.mean(upw_f_dqhi))


        upw_f.extend(upw_fhi)
        upw_f_dq.extend(upw_f_dqhi)

    label_prop = "upwc_q" + str(qorder)
    # mesh.addprop('upwc', np.array(upw_f))
    mesh.addprop(label_prop, np.array(upw_f_dq))
    return mesh, label_prop


def dEXP(x, y, z, data, shape, zmin, zmax, nlayers, qorder=0, SI=1):
    """
    DEXP model (Fedi, 2012). (NOT YET TESTED)

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
    mesh_dexp = _makemesh(x, y, shape, zmin, zmax, nlayers)
    # This way, if z is not an array, it is now
    z = z * np.ones_like(x)
    # Remove the last z because I only want depths to the top of the layers
    depths = mesh_dexp.get_zs()[:-1]
    weights = (np.abs(depths)) ** ((SI + qorder) / 2)
    csd = []
    # Offset by the data z because in the paper the data is at z=0
    for depth, weight in zip(depths - z[0], weights):

        # continued field calculation
        upw_f = transform.upcontinue(x, y, data, shape, depth)

        # qorder vertical derivate of the continued field
        upw_f_dq = transform.derivz(x, y, upw_f, shape, order=qorder)

        # the continued field weigted (=DEXP)
        upw_f_dq_w = upw_f_dq * weight
        csd.extend(upw_f_dq_w)

    label_prop = "dexp_q" + str(qorder)
    mesh_dexp.addprop(label_prop, np.array(csd))
    return mesh_dexp, label_prop


def dEXP_ratio(x, y, z, data, shape, zmin, zmax, nlayers, qorders=[1, 0], returnField=False):
    """
    DEXP ratio model (NOT YET validated)
    Abbas, M. A., and Fedi, M. (2014). Automatic DEXP
    imaging of potential fields independent of the structural index. Geophysical Journal
    International, 199 (3), 1625-1632.

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
    * qorders : 1D-array
        The order of the derivatives for the ratio calculation

    Returns:

    * mesh : :class:`fatiando.mesher.PrismMesh`
        The estimated physical property distribution set in a prism mesh (for
        easy 3D plotting)

    """
    mesh_dexp = _makemesh(x, y, shape, zmin, zmax, nlayers)
    # This way, if z is not an array, it is now
    z = z * np.ones_like(x)
    # Remove the last z because I only want depths to the top of the layers
    depths = mesh_dexp.get_zs()[:-1]
    weights = (np.abs(depths)) ** ((qorders[0] - qorders[1]) / 2)
    csd = []
    # Offset by the data z because in the paper the data is at z=0
    for depth, weight in zip(depths - z[0], weights):

        # continued field calculation
        upw_f = transform.upcontinue(x, y, data, shape, depth)

        # qorder vertical derivate of the continued field
        upw_f_dq_0 = transform.derivz(x, y, upw_f, shape, order=qorders[0])
        upw_f_dq_1 = transform.derivz(x, y, upw_f, shape, order=qorders[1])
        
        # upw_f_dq_0 = transform.derivy(x, y, upw_f, shape, order=qorders[0])
        # upw_f_dq_1 = transform.derivy(x, y, upw_f, shape, order=qorders[1])
        
        
        ratio = upw_f_dq_0 / upw_f_dq_1
        # the continued field weigted (=DEXP)
        upw_f_dq_w = weight * ratio
        # upw_f_dq_w =  ratio
        csd.extend(upw_f_dq_w)

    label_prop = "dexp_q" + str(qorders)
    mesh_dexp.addprop(label_prop, np.array(csd))
    
    if returnField:
        return mesh_dexp, label_prop, upw_f_dq_0, upw_f_dq_1, weight, upw_f_dq_w
    else:
        return mesh_dexp, label_prop




# def auto_dEXP():
def _jumpAnalysis(x):
    a = np.array(x)
    # Get forward and backward derivatives (their absolute value as we care about magnitude of differences not direction):
    d1 = np.r_[0, np.abs(a[1:] - a[:-1])]
    d2 = np.r_[np.abs(a[1:] - a[:-1]), 0]
    largest_jumps = np.argmax(d1)

    jump = []
    if d1[largest_jumps] > 20:
        jump.append(largest_jumps)

    return jump


def _removeOutliers(x, outlierConstant=1):
    a = np.array(x)
    # print(a)
    upper_quartile = np.percentile(a, 75)
    lower_quartile = np.percentile(a, 25)
    # print('upper_quartile: ' + str(upper_quartile))
    # print('lower_quartile: ' + str(lower_quartile))
    IQR = (upper_quartile - lower_quartile) * outlierConstant
    quartileSet = (lower_quartile - IQR, upper_quartile + IQR)
    # print("quartileSet: " + str(quartileSet))

    resultList = []
    iTrue = []
    for i, y in enumerate(a.tolist()):
        if y >= quartileSet[0] and y <= quartileSet[1]:
            resultList.append(y)
            iTrue.append(i)
    # print(iTrue)
    return resultList, iTrue, quartileSet


def _fit(x, y, **kwargs):
    """
    Curve least square fit.

    Parameters:

    *

    Returns:

    * Intersect y(0) (with the y-axis, useful for scalFUN function)

    """

    def f(x, A, B):  # this is your 'straight line' y=f(x)
        return A * x + B

    if np.count_nonzero(~np.isnan(x)) < 2:
        raise ValueError("Need at least 3 points to fit the data")

    else:
        try:

            for key, value in kwargs.items():
                if key == "rmvOutliers":
                    # print(x)
                    _, filtered_entries, quartileSet = _removeOutliers(x)
                    # print("------")
                    # print(filtered_entries)
                    # ------------------------------------ #
                    # print('remove outliers points')
                    df = pd.concat([x, y], axis=1).reindex(x.index)
                    # if np.mean(x)<0:
                    #     filtered_entries = (x < 0)
                    # else:
                    #     filtered_entries = (x > 0)
                    new_df = df.iloc[filtered_entries, :]

                    # new_df = df[filtered_entries]
                    x = new_df.iloc[:, 0]
                    y = new_df.iloc[:, 1]
                    # ------------------------------------ #

            popt, pcov = curve_fit(f, x, y)  # your data x, y to fit
            # popt, pcov = curve_fit(f,x,y,loss='soft_l1') # your data x, y to fit
            x_min = min(x)
            x_max = max(x)  # min/max values for x axis

            # x_min = quartileSet[0]
            # x_max = quartileSet[1]                               #min/max values for x axis
            # print('x_min: ' + str(x_min))
            # print('x_max: ' + str(x_max))

            for key, value in kwargs.items():
                if key == "slope":  # prolongation until the curve intersect
                    slope = value
                    if slope < 0:
                        x_min_reg = x_min + 1.2 * np.abs(x_max)
                        x_max_reg = x_max
                    if slope > 0:
                        x_min_reg = x_min - 2 * x_max
                        x_max_reg = x_max  # + 2*np.abs(max(df[0][k[1]]))

                    x_fit = np.linspace(
                        x_min_reg, x_max_reg, 100
                    )  # range of x values used for the fit function

                elif (
                    key == "SI"
                ):  # study of ridges intersection with the origin to infer the structural index
                    fit_SI = value
                    if fit_SI == True:
                        x_min_SI = 0
                        x_max_SI = max(x)  # min/max values for x axis
                        x_fit = np.linspace(
                            x_min_SI, x_max_SI, 100
                        )  # range of x values used for the fit function

            #     if key == 'xmin':
            #        x_min = value
            #     if key == 'xmax':
            #        x_max = value
            # print('x_min_reg: ' + str(x_min_reg))
            # print('x_max_reg: ' + str(x_max_reg))

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
    dfI = dfI.add_prefix("EX_xpos")
    dfI = dfI.rename(columns={"EX_xpos0": "elevation"})
    dfI.head(5)

    dfII = pd.DataFrame(RII_minmax)
    # df[0] = ['layer']
    dfII = dfII.add_prefix("EX_xpos")
    dfII = dfII.rename(columns={"EX_xpos0": "elevation"})
    dfII.head(5)

    dfIII = pd.DataFrame(RIII_minmax)
    # df[0] = ['layer']
    dfIII = dfIII.add_prefix("EX_xpos")
    dfIII = dfIII.rename(columns={"EX_xpos0": "elevation"})
    dfIII.head(5)

    return dfI, dfII, dfIII


def pad_edges(xp, yp, U, shape, pad_type=2):

    padtypes = [
        "0",
        "mean",
        "edge",
        "lintaper",
        "reflection",
        "oddreflection",
        "oddreflectiontaper",
    ]
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
    ax.pcolormesh(
        pad_y.reshape(padshape), pad_x.reshape(padshape), padded_data, cmap="RdBu_r"
    )
    ax.set_xlim(pad_y.min(), pad_y.max())
    ax.set_ylim(pad_x.min(), pad_x.max())

    shape = padded_data.shape
    U = padded_data.reshape(shape[0] * shape[1])
    xp = pad_x
    yp = pad_y

    return xp, yp, U, shape


def profile_extra(x, y, v, point1, point2, size, algorithm="cubic"):
    """
    Extract a profile between 2 points from spacial data.

    Uses interpolation to calculate the data values at the profile points.

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
    * algorithm : string
        Interpolation algorithm. Either ``'cubic'``, ``'nearest'``,
        ``'linear'`` (see scipy.interpolate.griddata).

    Returns:

    * [xp, yp, distances, vp] : 1d arrays
        ``xp`` and ``yp`` are the x, y coordinates of the points along the
        profile. ``distances`` are the distances of the profile points from
        ``point1``. ``vp`` are the data points along the profile.

    """

    return


def profile_noInter(x, y, v, point1, point2, size=None, **kwargs):
    # https://stackoverflow.com/questions/7878398/how-to-extract-an-arbitrary-line-of-values-from-a-numpy-array
    """
    Extract a profile between 2 points from spacial data.

    NO interpolation to calculate the data values at the profile points (find nearest point).

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

    Xaxis = []
    showfig = False
    for key, value in kwargs.items():
        if key == "Xaxis":
            Xaxis = value
        if key == "showfig":
            showfig = value

    x1, y1 = point1
    x2, y2 = point2
    maxdist = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

    if size is None:
        size = len(v)

    distances = np.linspace(0, maxdist, size)
    angle = np.arctan2(y2 - y1, x2 - x1)
    xp = x1 + distances * np.cos(angle)
    yp = y1 + distances * np.sin(angle)

    if Xaxis == "dist":
        xaxis = distances
    elif Xaxis == "y":
        xaxis = xp
    else:
        xaxis = yp

    nodes = np.array([x, y]).T
    points_p = np.array([xp, yp]).T
    # find nearest point

    # from progressbar import ProgressBar
    # pbar = ProgressBar()
    # vp = []
    # for p in pbar(points_p):
    #     ind = _closest_node(p, nodes)
    #     vp.append(v[ind])

    # from progressbar import ProgressBar
    # pbar = ProgressBar()
    vp = []
    for p in points_p:
        ind = _closest_node(p, nodes)
        vp.append(v[ind])

    vp_smooth_dict = _smooth_allfcts(xaxis, vp, showfig)

    return xp, yp, distances, vp, vp_smooth_dict


def _smooth_allfcts(xaxis, vp, showfig=False):

    # window_size, poly_order = 101, 3
    # vp_smooth = savgol_filter(vp, window_size, poly_order)
    # print(xaxis)
    # spl = UnivariateSpline(xaxis, vp, s=10)
    # plt.plot(xaxis, spl(xaxis), 'g', lw=3)
    # vp_smooth_spline = np.array(spl(xaxis))

    # xaxis = distances
    # xnew = np.linspace(min(xaxis),
    #                 max(xaxis),len(x))
    # vp_smooth = gridder.interp_at(x, y, v, xp, yp, algorithm='cubic', extrapolate=True)

    #%%
    # _smooth_lowpass

    #%%
    # if smooth == True:

    vp_smooth_0 = _smooth_lowpass(xaxis, vp)

    vp_smooth_1 = _smooth1d(xaxis, vp)
    vp_smooth_2 = _smooth1d_old(xaxis, vp)
    vp_smooth_3 = _smooth_lowpass(xaxis, vp_smooth_1)

    # f = interpolate.interp1d(xaxis, vp, fill_value='extrapolate',kind='cubic')
    # xnew = np.linspace(min(xaxis),
    #                     max(xaxis),len(xaxis))
    # vp_smooth_4 = f(xnew)
    from scipy.ndimage.filters import uniform_filter1d

    N = int(len(vp_smooth_2) / 20)
    vp_smooth_4 = uniform_filter1d(
        vp_smooth_2, size=N, mode="nearest", cval=min(vp_smooth_2)
    )

    if showfig == True:
        plt.figure()
        plt.plot(xaxis, vp, label="raw", marker="+")
        plt.plot(xaxis, vp_smooth_0, label="lowpass")
        plt.plot(xaxis, vp_smooth_1, label="hanning_window")
        plt.plot(xaxis, vp_smooth_3, label="hanning_window + lowpass")
        plt.plot(xaxis, vp_smooth_2, label="CubicSmoothingSpline")
        plt.plot(xaxis, vp_smooth_4, label="CubicSmoothingSpline + interp1d")
        plt.show()
        plt.legend()

    vp_smooth_dict = {
        "Lowpass": vp_smooth_0,
        "Hanning": vp_smooth_1,
        "Hanning+Lowpass": vp_smooth_3,
        "CubicSmoothingSpline": vp_smooth_2,
        "CubicSmoothingSpline + interp1d": vp_smooth_4,
    }

    return vp_smooth_dict


def _closest_node(node, nodes):
    nodes = np.asarray(nodes)
    dist_2 = np.sum((nodes - node) ** 2, axis=1)
    return np.argmin(dist_2)


# ----- peaks detections

from numpy import NaN, Inf, arange, isscalar, asarray, array


def peakdet(v, delta, x=None):
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
        sys.exit("Input vectors v and x must have same length")

    if not isscalar(delta):
        sys.exit("Input argument delta must be a scalar")

    if delta <= 0:
        sys.exit("Input argument delta must be positive")

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
            if this < mx - delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn + delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return np.array(maxtab), np.array(mintab)


def _smooth1d_old(x_axis, p_up):
    import csaps

    spl = csaps.CubicSmoothingSpline(x_axis, p_up, smooth=0.005)
    # spl = UnivariateSpline(x_axis, p_up, s=3)
    # spl.set_smoothing_factor(0.005)

    return np.array(spl(x_axis))


def _smooth1d(x_axis, p_up, window_len=6, window="hanning"):
    """smooth the data using a window with requested size."""

    # from scipy.signal import savgol_filter
    # from scipy.interpolate import interp1d

    # itp = interp1d(x,y, kind='linear')
    ws = int(len(p_up) / 3)
    if ws % 2 == 0:
        ws = ws + 1

    window_size, poly_order = ws, 2
    filtdata = savgol_filter(p_up, window_size, poly_order)

    return filtdata


def _smooth_lowpass(x_axis, p_up, **kwargs):

    cutoff = 0.015  # desired cutoff frequency of the filter, Hz ,      slightly higher than actual 1.2 Hz
    for key, value in kwargs.items():
        if key == "cutoff":
            cutoff = value

    # Filter requirements.
    fs = abs(1 / (x_axis[0] - x_axis[1]))  # sample rate, Hz
    nyq = 0.5 * fs  # Nyquist Frequency
    order = 1  # sin wave can be approx represented as quadratic
    # print('cutoff '+ str(cutoff))
    # print('fs '+ str(fs))
    # print('order '+ str(order))

    filtdata = butter_lowpass_filter(p_up, cutoff, nyq, order)

    return np.array(filtdata)


def butter_lowpass_filter(p_up, cutoff, nyq, order):
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients
    b, a = butter(order, normal_cutoff, btype="low", analog=False)
    y = filtfilt(b, a, p_up)
    return y


def smooth2d(x, y, U, sigma=10,plot=False):

    plt.figure()
    plt.subplot(1, 2, 1)
    plt.scatter(x, y, c=U, cmap="viridis", vmax=max(U), vmin=min(U))
    plt.colorbar()
    plt.axis("square")
    # plt.show()

    U2d = U.reshape(int(np.sqrt(U.shape)), int(np.sqrt(U.shape)))
    U2d_f = gaussian_filter(U2d, sigma=sigma)

    U_f = np.copy(U)
    U_f = U2d_f.reshape(U.shape)
    
    if plot:
        plt.subplot(1, 2, 2)
        plt.scatter(x, y, c=U_f, cmap="viridis", vmax=max(U), vmin=min(U))
        plt.colorbar()
        plt.axis("square")
    plt.show()

    return U_f
