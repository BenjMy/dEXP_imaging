"""
Sensitivity analysis of DEXP to anomaly depth on Mise-a-la-masse 
----------------------------------------------------------------

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

Rucker, C., Gunther, T., Wagner, F.M., 2017. pyGIMLi: An open-source library for modelling and inversion in geophysics, Computers and Geosciences, 109, 106-123, doi: 10.1016/j.cageo.2017.07.011

"""
#%%

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

import numpy as np

#%%
# MALM DATA synthetic anomaly: analysis of sensitivity

MESH = []
LABEL = []
MESHratio = []
LABELratio = []
DF_F = []
DF_FIT = []
XXZZ = []
CTm = []
UPW_f_dq_0 = []
UPW_f_dq_1 = []

filenames = ['MSoilR1000.0AnoR1Z-3.75L5h2.5',
              'MSoilR1000.0AnoR1Z-13.75L5h2.5',
              'MSoilR1000.0AnoR1Z-23.75L5h2.5']



# import pickle
# infile = open('./loadmalm/' + filenames[0] + '.pkl','rb')
# data = pickle.load(infile)
# infile.close()
    
    
# import pickle
# filenames = ['MSoilR1000.0AnoR1Z-3.75L5h2.5']
x_axis='y'
for fi in filenames:
    x_raw, y_raw, z_raw, U_raw, maxdepth, shape_raw, p1, p2, SimName, ano_prop = MALM.load_MALM_sens3d(filename='./loadmalm/' +
                                                                fi + '.pkl')

    shape = (200,200)
    xp,yp,U = gridder.interp(x_raw,y_raw,U_raw,shape)
    
    
    parameters = para.set_par(shape=shape,max_elevation=abs(maxdepth))
    interp = True
    scaled = parameters[0]
    SI = parameters[1]
    zp, qorder, nlay = parameters[2:5]
    minAlt_ridge, maxAlt_ridge = parameters[5:7]
    
    #%%
    # ridges analysis parameters
    nlay = 25
    max_elevation = 30
    minAlt_ridge = max_elevation*0.05
    maxAlt_ridge = max_elevation*0.65
    
    interp = True
    smooth = False 
    
    #%%
    # Anomalies properties
    # HDWL : height, Depth, Width (x), Lenght (y)

    x1, x2, z1, z2 = [max(x_raw)/2-ano_prop['HWD'][1]/2,max(x_raw)/2 + ano_prop['HWD'][1]/2,
                    ano_prop['HWD'][2]+ ano_prop['HWD'][0]/2,
                    ano_prop['HWD'][2]- ano_prop['HWD'][0]/2]
    xxzz = [x1, x2, z1, z2]
    CT = ano_prop['SoilR']/ano_prop['AnoR']
      
    #%% 

    qratio = [2,1]
    mesh_ratio, label_ratio, upw_f_dq_0, upw_f_dq_1, weight, upw_f_dq_w  = dEXP.dEXP_ratio(xp, yp, zp, U, shape, 
                                                                              zmin=0, zmax=max_elevation, nlayers=nlay, 
                                                                              qorders=qratio,returnField=True)

    #%% 
    # save data loop

    MESHratio.append(mesh_ratio)
    UPW_f_dq_0.append(upw_f_dq_0)
    UPW_f_dq_1.append(upw_f_dq_1)
    MESHratio.append(mesh_ratio)
    LABELratio.append(label_ratio)
    XXZZ.append(xxzz)
    CTm.append(CT)


    ratio = UPW_f_dq_0[0] / UPW_f_dq_1[0]
    fig, ax = plt.subplots(2,2, sharex=True, sharey=True)
    pEXP.plot_field(xp, yp, U, shape, ax=ax[0,0])
    ax[0,0].set_title('Field u')
    pEXP.plot_field(xp, yp, UPW_f_dq_0[0], shape, ax=ax[0,1])
    ax[0,1].set_title('Field derivate' + str(qratio[0]))
    pEXP.plot_field(xp, yp, UPW_f_dq_1[0], shape, ax=ax[1,0])
    ax[1,0].set_title('Field derivate' + str(qratio[1]))
    pEXP.plot_field(xp, yp, upw_f_dq_w, shape, ax=ax[1,1])
    ax[1,1].set_title('Scaled Field ratio')
    selec_x = list(np.arange(70,200-70))

#%% 
# Plot the results
#
# .. important::
# 
#     True depth are respectively 3.75, 13.75, 23.75
# 


scl = 0
i = 0

#%% 

fig = plt.figure()
ax = plt.gca()
plt, cmap = pEXP.plot_xy(MESHratio[i],scaled=scl, label=LABELratio[i],
              markerMax=True,qratio=str(qratio),aspect_equal=True,
              ax=ax, Xaxis=x_axis,p1p2=np.array([p1,p2]),Vminmax=[-0.4,0.4],
              regional_cut = selec_x) 
# plt.colorbar(cmap)
x1, x2, z1, z2 = XXZZ[i]
square([x1, x2, -z1, -z2])
plt.annotate(CTm[i],[(x1 + x2)/2, -(z1+z2)/2])
# plt.tight_layout()

#%% 
i = 1

fig = plt.figure()
ax = plt.gca()
plt, cmap = pEXP.plot_xy(MESHratio[i],scaled=scl, label=LABELratio[i],
              markerMax=True,qratio=str(qratio),aspect_equal=True,
              Vminmax=[-0.4,0.4],
              ax=ax, Xaxis=x_axis,p1p2=np.array([p1,p2]),
              regional_cut = selec_x) 
# plt.colorbar(cmap)
x1, x2, z1, z2 = XXZZ[i]
square([x1, x2, -z1, -z2])
plt.annotate(CTm[i],[(x1 + x2)/2, -(z1+z2)/2])
# plt.tight_layout()

#%% 
i = 2

fig = plt.figure()
ax = plt.gca()
plt, cmap = pEXP.plot_xy(MESHratio[i], scaled=scl, label=LABELratio[i],
              markerMax=True,qratio=str(qratio),aspect_equal=True,
              Vminmax=[-0.4,0.4],
              ax=ax, Xaxis=x_axis,p1p2=np.array([p1,p2]),
              regional_cut = selec_x) 
# plt.colorbar(cmap)
# plt.annotate('True depth is', xy, args, kwargs)()
x1, x2, z1, z2 = XXZZ[i]
square([x1, x2, -z1, -z2])
plt.annotate(CTm[i],[(x1 + x2)/2, -(z1+z2)/2])
# plt.title('23.75m depth anomaly')
# plt.tight_layout()

