"""
Sensitivity analysis of DEXP to anomaly width on Mise-Ã -la-masse 
-----------------------------------------------------------------

This code shows a step-by-step processing of potential field imaging aiming at giving an estimate of electrical sources positions and depth using the dEXP tranformation method.
dEXP method implementation from Fedi et al. 2012. 
Calculations used :mod:`dEXP`, while plotting use the :mod:`plot_dEXP` module.

Application on a anomaly of electrical resistivity.
The model data was created using geometric objects from :mod:`pygimli.meshtools`. The forward simulation of the data was done using :mod:`pygimli.ERTsimulate` module.


.. note::

    This is part of a larger project aiming at inverting current sources density (see more at: https://icsd-dev.readthedocs.io/en/latest/)


**References**

Uieda, L., V. C. Oliveira Jr, and V. C. F. Barbosa (2013), Modeling the Earth with Fatiando a Terra, Proceedings of the 12th Python in Science Conference, pp. 91 - 98.

Uieda, L. (2018). Verde: Processing and gridding spatial data using Green's functions. Journal of Open Source Software, 3(29), 957. doi:10.21105/joss.00957

Fedi, M., and M. Pilkington (2012), Understanding imaging methods for potential
field data, Geophysics, 77(1), G13, doi:10.1190/geo2011-0078.1

Rücker, C., Günther, T., Wagner, F.M., 2017. pyGIMLi: An open-source library for modelling and inversion in geophysics, Computers and Geosciences, 109, 106-123, doi: 10.1016/j.cageo.2017.07.011

"""

from fatiando.vis.mpl import square
from fatiando import gridder

# my own functions
import lib.dEXP as dEXP
import lib.plot_dEXP as pEXP
import lib.set_parameters as para

# exemples
import examples.malm.loadmalm.Load_sens_MALM as MALM

import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 15


#%%
# MALM DATA synthetic anomaly: analysis of sensitivity

MESH = []
LABEL = []
DF_F = []
DF_FIT = []
XXZZ = []
CTm = []


filenames = ['MSoilR1000.0AnoR1Z-13.75W5H2.5L5Noise0',
             'MSoilR1000.0AnoR1Z-13.75W25H2.5L5']

for fi in filenames:
    x_raw, y_raw, z_raw, U_raw, maxdepth, shape_raw, p1, p2, _ , ano_prop = MALM.load_MALM_sens3d(filename='./loadmalm/' +
                                                                fi + '.pkl')
    shape = (150,150)
    xp,yp,U = gridder.interp(x_raw,y_raw,U_raw,shape)
    
    
    parameters = para.set_par(shape=shape,max_elevation=abs(maxdepth))
    interp = True
    scaled = parameters[0]
    SI = parameters[1]
    zp, qorder, nlay = parameters[2:5]
    minAlt_ridge, maxAlt_ridge = parameters[5:7]
    x_axis = 'y'
    #%%
    # ridges analysis parameters
    nlay = 25
    max_elevation = 20
    minAlt_ridge = max_elevation*0.05
    maxAlt_ridge = max_elevation*0.65
    
    interp = True
    smooth = False 
    
    #%%
    # Anomalies properties
    # HDWL : height, Depth, Width (x), Lenght (y)
    x1, x2, z1, z2 = [max(x_raw)/2-ano_prop['HWDL'][1]/2,max(x_raw)/2 + ano_prop['HWDL'][1]/2,
                    ano_prop['HWDL'][2]+ ano_prop['HWDL'][0]/2,
                    ano_prop['HWDL'][2]- ano_prop['HWDL'][0]/2]
    xxzz = [x1, x2, z1, z2]
    CT = ano_prop['SoilR']/ano_prop['AnoR']
    
    #%% 
    # Plot the data 
    pEXP.plot_line(xp, yp, U,p1,p2, interp=interp,Xaxis=x_axis)
        
    #%% 
    # Pad the edges of grids (if necessary)
    # xp,yp,U, shape = dEXP.pad_edges(xp,yp,U,shape,pad_type=0) # reflexion=5
    # pEXP.plot_line(xp, yp,U,p1,p2, interp=interp)
        
        
    #%% 
    # Upward continuation of the field data
    
    mesh, label_prop = dEXP.upwc(xp, yp, zp, U, shape, 
                     zmin=0, zmax=max_elevation, nlayers=nlay, 
                     qorder=qorder)
    
    plt, cmap = pEXP.plot_xy(mesh, label=label_prop,Xaxis=x_axis)
    plt.colorbar(cmap, label=label_prop)
    
    
    #%%
    # Ridges identification
    # dEXP.ridges_minmax_plot(xp, yp, mesh, p1, p2,
    #                                       label=label_prop,
    #                                       fix_peak_nb=2,
    #                                       method_peak='find_peaks')  
    
    # or  find_peaks or peakdet or spline_roots
    dfI,dfII, dfIII = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,
                                          label=label_prop,
                                          fix_peak_nb=2,
                                          method_peak='find_peaks',
                                          Xaxis=x_axis,
                                          showfig=True) 
    
    df = dfI, dfII, dfIII

    #%% 
    # Plot ridges over continuated section
        
    # fig = plt.figure()
    # ax = plt.gca()
    # pEXP.plot_xy(mesh, label=label_prop, ax=ax) #, ldg=)
    # pEXP.plot_ridges_harmonic(dfI,dfII,dfIII,ax=ax)
    
    #%%
    # Filter ridges (regionally constrainsted)
    
    dfI_f,dfII_f, dfIII_f = dEXP.filter_ridges(dfI,dfII,dfIII,
                                                minDepth=minAlt_ridge,maxDepth=maxAlt_ridge,
                                                minlength=10,rmvNaN=True,
                                                Xaxis=x_axis)
    df_f = dfI_f, dfII_f, dfIII_f
    
    #%%
    # Fit ridges (regionally constrainsted)

    # df_fit = dEXP.fit_ridges(df_f, rmvOutliers=True) # fit ridges on filtered data
    df_fit = dEXP.fit_ridges(df_f, rmvOutliers=True) # fit ridges on filtered data

    #%% 
    # save data loop

    MESH.append(mesh)
    LABEL.append(label_prop)
    DF_F.append(df_f)
    DF_FIT.append(df_fit)
    XXZZ.append(xxzz)
    CTm.append(CT)




#%%

i = 0
dfI_f,dfII_f,dfIII_f = DF_F[i]
fig, ax1 = plt.subplots(figsize=(15,3))
x1, x2, z1, z2 = XXZZ[i]
square([x1, x2, z1, z2])
plt, cmap = pEXP.plot_xy(MESH[i], label=LABEL[i], ax=ax1, Xaxis='y',
          Vminmax=[0,10])
pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax1,label=True)
ax2 = pEXP.plot_ridges_sources(DF_FIT[i], ax=ax1, z_max_source=-max_elevation*1.2,
                      ridge_type=[0,1,2],ridge_nb=None)


#%% 

i = 1
dfI_f,dfII_f,dfIII_f = DF_F[i]
fig, ax1 = plt.subplots(figsize=(15,3))
x1, x2, z1, z2 = XXZZ[i]
square([x1, x2, z1, z2])
plt, cmap = pEXP.plot_xy(MESH[i], label=LABEL[i], ax=ax1, Xaxis='y',
          Vminmax=[0,10])
pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax1,label=True)
ax2 = pEXP.plot_ridges_sources(DF_FIT[i], ax=ax1, z_max_source=-max_elevation*1.2,
                      ridge_type=[0,1,2],ridge_nb=None)
