"""
Sensitivity analysis of DEXP to noise level  on Mise-Ã -la-masse 
----------------------------------------------------------------

This code shows a step-by-step processing of potential field imaging aiming at giving an estimate of electrical sources positions and depth using the dEXP tranformation method.
dEXP method implementation from Fedi et al. 2012. 
Calculations used :mod:`dEXP`, while plotting use the :mod:`plot_dEXP` module.

Application on a anomaly of electrical resistivity.
The model data was created using geometric objects from :mod:`pygimli.meshtools`. The forward simulation of the data was done using :mod:`pygimli.ERTsimulate` module.
The noise is added to the data. Here 1% plus 1µV. Note, we force a specific noise seed as we want reproducable results for testing purposes.

.. note::

    This is part of a larger project aiming at inverting current sources density (see more at: https://icsd-dev.readthedocs.io/en/latest/)


**References**

Uieda, L., V. C. Oliveira Jr, and V. C. F. Barbosa (2013), Modeling the Earth with Fatiando a Terra, Proceedings of the 12th Python in Science Conference, pp. 91 - 98.

Uieda, L. (2018). Verde: Processing and gridding spatial data using Green's functions. Journal of Open Source Software, 3(29), 957. doi:10.21105/joss.00957

Fedi, M., and M. Pilkington (2012), Understanding imaging methods for potential
field data, Geophysics, 77(1), G13, doi:10.1190/geo2011-0078.1

Rucker, C., Gunther, T., Wagner, F.M., 2017. pyGIMLi: An open-source library for modelling and inversion in geophysics, Computers and Geosciences, 109, 106-123, doi: 10.1016/j.cageo.2017.07.011

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


# filenames = ['MSoilR1000.0AnoR1Z-13.75W5H2.5L5Noise0',
#              'MSoilR1000.0AnoR1Z-13.75W5H2.5L5Noise0.01', # 10% of noise
#              'MSoilR1000.0AnoR1Z-13.75W5H2.5L5Noise0.02', # 20% of noise
#              'MSoilR1000.0AnoR1Z-13.75W5H2.5L5Noise0.02', # use smooth fct
#              'MSoilR1000.0AnoR1Z-13.75W5H2.5L5Noise0.02'] # use derivative order 1




filenames = ['MSoilR10AnoR10Z-13.75W5H2.5L5S0Noise0',
             'MSoilR100AnoR1Z-13.75W5H2.5L5S0Noise10', # 10% of noise
             'MSoilR100AnoR1Z-13.75W5H2.5L5S0Noise20', # 20% of noise
             'MSoilR100AnoR10Z-13.75W5H2.5L5S0Noise0'] # use derivative order 1

for i, fi in enumerate(filenames):
    x_raw, y_raw, z_raw, U_raw, maxdepth, shape_raw, p1, p2, SimName, ano_prop = MALM.load_MALM_sens3d(filename='./loadmalm/' +
                                                                fi + '.pkl')
    shape = (150,150)
    # print(U_raw[1:10])
    xp,yp,U = gridder.interp(x_raw,y_raw,U_raw,shape)
    
    
    parameters = para.set_par(shape=shape,max_elevation=abs(maxdepth))
    interp = True
    scaled = parameters[0]
    SI = parameters[1]
    zp, qorder, nlay = parameters[2:5]
    minAlt_ridge, maxAlt_ridge = parameters[5:7]
    fix_peak_nb = 2
    #%%
    # ridges analysis parameters
    nlay = 25
    max_elevation = 20
    minAlt_ridge = max_elevation*0.25
    maxAlt_ridge = max_elevation*0.75
    
    interp = True
    smooth = False 
    x_axis = 'y'

    #%%
    # Anomalies properties
    # HDWL : height, Depth, Width (x), Lenght (y)
    x1, x2, z1, z2 = [max(x_raw)/2-ano_prop['HWD'][1]/2,max(x_raw)/2 + ano_prop['HWD'][1]/2,
                    ano_prop['HWD'][2]+ ano_prop['HWD'][0]/2,
                    ano_prop['HWD'][2]- ano_prop['HWD'][0]/2]
    xxzz = [x1, x2, z1, z2]
    CT = ano_prop['SoilR']/ano_prop['AnoR']
    
    #%% 
    # some constrainsts
    if i==2:
    # smooth the data 
        U = dEXP.smooth2d(xp, yp, U, sigma=1)
        # smooth = True 

    if i==3:
    # inscrease derivative order 
        qorder = qorder + 1
        # fix_peak_nb = 4
        U = dEXP.smooth2d(xp, yp, U, sigma=1)
    #%% 
    # Plot the data 
    # pEXP.plot_line(xp, yp, U,p1,p2, interp=interp)
    
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
    dfI,dfII, dfIII, ax = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,
                                          label=label_prop,
                                          fix_peak_nb=fix_peak_nb,
                                          method_peak='find_peaks',
                                          smooth=smooth,
                                          Xaxis=x_axis,
                                          showfig=True
                                          )  
        
    #%%
    # Filter ridges (regionally constrainsted)
    
    dfI_f,dfII_f, dfIII_f = dEXP.filter_ridges(dfI,dfII,dfIII,
                                                minDepth=minAlt_ridge,maxDepth=maxAlt_ridge,
                                                minlength=7,rmvNaN=True,
                                                xmin=100, xmax=300,
                                                Xaxis=x_axis
                                                )
    df_f = dfI_f, dfII_f, dfIII_f     
        
    #%%
    # fit 
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
fig, ax = plt.subplots(figsize=(15,3))
# ax = plt.gca()
pEXP.plot_xy(MESH[i], label=LABEL[i], ax=ax) #, ldg=)
dfI_f,dfII_f,dfIII_f = DF_F[i]
pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax,label=False,legend=False)   
pEXP.plot_ridges_sources(DF_FIT[i], ax=ax, z_max_source=-max_elevation*1.2,
                          ridge_type=[0,1,2],ridge_nb=None)
# plt.xlim([100,300])
x1, x2, z1, z2 = XXZZ[i]
square([x1, x2, z1, z2])
plt.annotate(CTm[i],[(x1 + x2)/2, (z1+z2)/2])
plt.title('0% of noise')

#%% 

i = 1
fig, ax = plt.subplots(figsize=(15,3))
pEXP.plot_xy(MESH[i], label=LABEL[i], ax=ax) #, ldg=)
dfI_f,dfII_f,dfIII_f = DF_F[i]
pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax,label=False)   
pEXP.plot_ridges_sources(DF_FIT[i], ax=ax, z_max_source=-max_elevation*1.2,
                          ridge_type=[0,1,2],ridge_nb=None)
x1, x2, z1, z2 = XXZZ[i]
square([x1, x2, z1, z2])
plt.annotate(CTm[i],[(x1 + x2)/2, (z1+z2)/2])
plt.title('10% of noise')
plt.savefig('10_noise_q0.png', dpi=450)

i = 2
fig, ax = plt.subplots(figsize=(15,3))
pEXP.plot_xy(MESH[i], label=LABEL[i], ax=ax) #, ldg=)
dfI_f,dfII_f,dfIII_f = DF_F[i]
pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax,label=False)   
pEXP.plot_ridges_sources(DF_FIT[i], ax=ax, z_max_source=-max_elevation*1.2,
                          ridge_type=[0,1,2],ridge_nb=None)
x1, x2, z1, z2 = XXZZ[i]
square([x1, x2, z1, z2])
plt.annotate(CTm[i],[(x1 + x2)/2, (z1+z2)/2])
plt.title('10% of noise - smoothed')
plt.savefig('10_noise_q0_smooth.png', dpi=450)


i = 3
fig, ax = plt.subplots(figsize=(15,3))
pEXP.plot_xy(MESH[i], label=LABEL[i], ax=ax) #, ldg=)
dfI_f,dfII_f,dfIII_f = DF_F[i]
pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax,label=False)   
pEXP.plot_ridges_sources(DF_FIT[i], ax=ax, z_max_source=-max_elevation*1.2,
                          ridge_type=[0,1,2],ridge_nb=None)
x1, x2, z1, z2 = XXZZ[i]
square([x1, x2, z1, z2])
plt.annotate(CTm[i],[(x1 + x2)/2, (z1+z2)/2])
plt.title('10% of noise - q-order=1')

plt.tight_layout(pad=0.2, h_pad=-8.2)
plt.savefig('10_noise_q1.png', dpi=450)

# # Loop on source depth
# fig, axs = plt.subplots(1,len(filenames))

# for i in range(len(filenames)):
    
#     pEXP.plot_xy(MESH[i], label=LABEL[i], ax=axs[i]) #, ldg=)
#     dfI_f,dfII_f,dfIII_f = DF_F[i]
#     pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=axs[i],label=False, legend=False)    
#     pEXP.plot_ridges_sources(DF_FIT[i], ax=axs[i], z_max_source=-max_elevation*1.2,
#                                   ridge_type=[0,1,2],ridge_nb=None)
    # x1, x2, z1, z2 = XXZZ[i]
    # square([x1, x2, z1, z2])
    # plt.annotate(CTm[i],[(x1 + x2)/2, -(z1+z2)/2])
