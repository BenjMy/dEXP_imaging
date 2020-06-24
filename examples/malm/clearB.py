"""
MALM DEXP preprocessing
-----------------------

This code shows how to remove the influence of the return electrode B and correct as much as it is possible the field before dEXP analysis.
Calculations used :mod:`utils_dEXP`, while plotting use the :mod:`plot_dEXP` module.

Application on a anomaly of electrical resistivity.
The model data was created using geometric objects from :mod:`pygimli.meshtools`. The forward simulation of the data was done using :mod:`pygimli.ERTsimulate` module.


.. note::

    This is part of a larger project aiming at inverting current sources density (see more at: https://icsd-dev.readthedocs.io/en/latest/)


**References**

Rücker, C., Günther, T., Wagner, F.M., 2017. pyGIMLi: An open-source library for modelling and inversion in geophysics, Computers and Geosciences, 109, 106-123, doi: 10.1016/j.cageo.2017.07.011

----
"""

import utils_dEXP as uEXP
from fatiando import gridder
import matplotlib.pyplot as plt

import examples.malm.loadmalm.Load_sens_MALM as MALM

#%%
# Import MALM data total field U_{T} = U{a} - U{b}
# The total voltage is the superposition of the injection +I on electrode A (masse) and injection -I on the return (and remote) electrode B

#%%
# Load data
xp, yp, z, uA, maxdepth, shape, p1, p2, SimName, model_prop = MALM.load_MALM_sens3d(filename='./loadmalm/' +
                                                            'MSoilR1000.0AnoR1Z-23.75L5h2.5.pkl')
Bpos= model_prop['EA_EB_EN'][1]
uT_elecs= model_prop['u_elecs']
xyzs = model_prop['elecs']

#%%
# remove B return electrode effects using :mod:`uEXP.cor_field_B` using a resistivity of 100 Ohm.m

U_cor = uEXP.cor_field_B(xyzs[:-3,0],xyzs[:-3,1],xyzs[:-3,2],uT_elecs,Bpos,rho=100)

#%%
# Compare the voltage seen by the surface electrode using 2d plot over a line
p1_elecs= [min(xyzs[:-3,0]),(max(xyzs[:-3,0])+min(xyzs[:-3,0]))/2]
p2_elecs= [max(xyzs[:-3,0]),(max(xyzs[:-3,0])+min(xyzs[:-3,0]))/2]

xx, yy, distance, profile = gridder.profile(xyzs[:-3,0],xyzs[:-3,1], uT_elecs, p1_elecs, p2_elecs, 1000)
plt.figure()
# plt.title(strname + '_data' + str(ZZ), fontsize=15)
plt.plot(xx, profile, '.k')

xx, yy, distance, profile = gridder.profile(xyzs[:-3,0],xyzs[:-3,1], U_cor, p1_elecs, p2_elecs, 1000)
plt.plot(xx, profile, '.r')
plt.grid()

#%%
# Compare with numerical solution
uT = model_prop['uTz0_grid']
uT_field_cor = uEXP.cor_field_B(xp,yp,z,uT,Bpos,rho=100)

#%%
# Compare plot over a line
xx, yy, distance, puT = gridder.profile(xp,yp, uT, p1, p2, 1000)

plt.figure()
plt.subplot(2,1,1)
plt.plot(xx, puT, '*b', label='U_{T}')

xx, yy, distance, puA = gridder.profile(xp,yp, uA, p1, p2, 1000)
plt.plot(xx, puA, '.k', label='U_{A}')

xx, yy, distance, puT_field_cor = gridder.profile(xp,yp, uT_field_cor, p1, p2, 1000)
plt.plot(xx, puT_field_cor, '.r', label='U_{cor} = U_{T} - U_{B}')
plt.grid()
plt.legend()

plt.subplot(2,1,2)
diff1 = puT - puA
diff2 = puT - puT_field_cor
plt.plot(xx, diff1, '.b', label='U_{T} - U_{A}')
plt.plot(xx, diff2, '.g', label='U_{T} - U_{cor}')

plt.legend()


