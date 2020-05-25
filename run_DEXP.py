import os

import dEXP as dEXP
import numpy as np
import matplotlib.pyplot as plt


SI = 2
scaled=0
# -------------------------------  Model parameters
ZZ = -3.75
maxdepth=20
Rsoil = 1000
x1, x2, y1, y2, z1, z2 = -5,5,-5,5,ZZ-2.5/2,ZZ+2.5/2
#square([y1, y2, x1, x2])

filename = 'MSoilR' + str(Rsoil) + 'AnoR1Z'+ str(ZZ) +'L5h2.5'
MainPath='E:/Padova/Simulation/MALM_SensitivityAnalysis/Sens3d/' + filename + '/Data/'
#MainPath='E:/Padova/Simulation/MALM_SensitivityAnalysis/Sens3d/' + filename + '/Data/'
os.chdir(MainPath)

# ------------------------------- Load data
shape = (50, 50)
x,y,z,gz=np.loadtxt('3d_SENS_FAT.txt', unpack=True)

# ------------------------------- remove B return electrode effects 
B = [-104.16666666666667, -104.16666666666667, 0]
gz_cor = dEXP.cor_field_B(x,y,z,gz,B,rho=100)


# ------------------------------- Plot the data
plt.figure()
plt.axis('scaled')
mpl.contourf(yp, xp, gz, shape, 30)
plt.colorbar()
plt.xlabel('East (km)')
plt.ylabel('North (km)')
mpl.m2km()
plt.show()
#plt.title(strname +'ztop' + str(za) +'_zbot'+ str(zb) + '_data', fontsize=20)
#plt.savefig(pathFig+strname + '_ExempleFig_z' + str(za) + str(zb) + '_data' + '.png')
#x1, x2, y1, y2, z1, z2 = np.array(model[0].get_bounds())
square([y1, y2, x1, x2])

