"""
GravMag: 3D imaging on synthetic gravity data
(simple example)
"""

from fatiando import gridder, mesher, utils
from fatiando.gravmag import prism, imaging
from fatiando.vis import mpl #, myv
from fatiando.vis.mpl import square
import numpy as np
import matplotlib.pyplot as plt
import os

pathData='./'+ 'Data/'
# pathFig='./' + strname + '/'


# -------------------------------  Model
za= 1000
zb= 1500
dens = 1200
simname = 'za_' + str(za) + 'zb_' + str(zb) + 'dens_' + str(dens) 
model = [mesher.Prism(-1000, 1000, -1000, 1000, za, zb, {'density': dens})]
#model = [mesher.Prism(-4000, 0, -4000, -2000, za, zb, {'density': 1200}),
#         mesher.Prism(-1000, 1000, -1000, 1000, 5, 7000, {'density': -800})]

shape = (30, 30)
#xp, yp, zp = gridder.scatter((-6000, 6000, -6000, 6000), shape[0]*shape[1], z=0)
xp, yp, zp = gridder.regular((-6000, 6000, -6000, 6000), shape, z=0)

# gz = utils.contaminate(prism.gz(xp, yp, zp, model), 0.1)
gz = prism.gz(xp, yp, zp, model)

x1, x2, y1, y2, z1, z2 = np.array(model[0].get_bounds())/1000


import pickle
d = { "xyzg" : [xp, yp, zp, gz], "shape" : shape , "model" : model}
afile = open(simname + '.pkl', 'wb')
pickle.dump(d, afile)
afile.close()

import pickle

#reload object from file

    
file2 = open(simname + '.pkl', 'rb')
u = pickle._Unpickler(file2)
u.encoding = 'latin1'
test = u.load()

file2.close()

# new_d['model']



# ------------------------------- Plot the data
plt.figure()
plt.axis('scaled')
mpl.contourf(yp, xp, gz, shape, 30)
plt.colorbar()
plt.xlabel('x (m)')
plt.ylabel('y (m)')
mpl.m2km()
plt.show()
plt.tight_layout()
# plt.title(strname + '_Z_' + str(ZZ) + '_data')
# plt.savefig(pathFig + strname + '_Z_' + str(ZZ) + '_data' + '.png')
square([y1, y2, x1, x2])


# model
# data2write=np.zeros([len(distance),4])
# data2write[:,0:3]=xp, yp, zp
# data2write[:,4]= gz
