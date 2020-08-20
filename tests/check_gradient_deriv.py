# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 13:54:17 2020

@author: Benjamin
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 17:55:57 2020

@author: Benjamin
"""

import matplotlib.pyplot as plt
import numpy as np
import dEXP as dEXP
from dEXP import _fit
import plot_dEXP as pEXP
import set_parameters as para
from fatiando.gravmag import prism, imaging, transform
from fatiando.vis.mpl import square

# exemples
import examples.magnetic.fwdmag.fwd_mag_sphere as magfwd
import examples.gravimetry.loadgrav.grav_models as grav
import utils_dEXP as uEXP

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 15


#%% 
# load model previously generated using Fatiando a terra package
# asymetric square: za3000_zb3500_l500_ofs3000_dens1200
# asymetric rectangle: mod_asy_rect.pkl or mod_asy_rect_z1.2km or mod_asy_rect_z3.2km
# symetric:  za_1000zb_1500dens_1200 # or sym_square_z3.2km
# a/symetric with a rectangle space around sym_square_z3.2km_rectArea / mod_asy_rect_z3.2km_rectArea
# za3000_zb3500_l500_ofs0_dens1200

dataname = 'za_1000zb_1500dens_1200'
data_struct = grav.load_grav_fatiando(name='../examples/gravimetry/loadgrav/' + dataname)
xp,yp,zp,U = data_struct['xyzg']
shape = data_struct['shape']
model = data_struct['model']
dens  = data_struct['density']
# scaled, SI, zp, qorder, nlay, minAlt_ridge, maxAlt_ridge = para.set_par(shape=shape,max_elevation=max_elevation)

x1, x2, y1, y2, z1, z2 = np.array(model[0].get_bounds())

# p1 =[min(yp),0]
# p2 =[max(yp),0]

max_elevation=z2*1.2
scaled, SI, zp, qorder, nlay, minAlt_ridge, maxAlt_ridge = para.set_par(shape=shape,max_elevation=max_elevation)
interp = True
qorder = 1


for slicedir in enumerate('x'):
    if slicedir[1]=='x':
        print(slicedir[1])
        p1 =[(x1+x2)/2,min(yp)]
        p2 =[(x1+x2)/2,max(yp)]
        x_axis = slicedir[1] # direction of p1p2 profil
        SI = 1.5
    else:
        p1 =[min(xp),(y1+y2)/2]
        p2 =[max(xp),(y1+y2)/2]
        x_axis = slicedir[1] # direction of p1p2 profil
        SI = 2

    #%%
    # Plot field data over a 2d line crossing the anomalies
    # pEXP.plot_line(xp, yp, U,p1,p2, interp=False, Xaxis='x')
    
    pEXP.plot_line(xp, yp, U,p1,p2, interp=True, Xaxis=x_axis)
    square([x1, x2, y1, y2])

    pEXP.plot_field(xp,yp,U, shape)
    square([x1, x2, y1, y2])
    
# Gaussian function to derivate 
    #%% ------------------------------- Plot the derivatives
    # http://campar.in.tum.de/Chair/HaukeHeibelGaussianDerivatives
    
    xderiv = transform.derivx(xp, yp, U, shape,order=qorder)
    yderiv = transform.derivy(xp, yp, U, shape,order=qorder)
    zderiv = transform.derivz(xp, yp, U, shape,order=qorder)
    
    # # plt.plot(xderiv)
    # # plt.plot(yderiv)
    
    # # interp = True
    pEXP.plot_line(xp, yp, xderiv ,p1,p2,title='xderiv',savefig=False, interp=interp, Xaxis=x_axis)
    pEXP.plot_line(xp, yp, yderiv ,p1,p2,title='yderiv',savefig=False, interp=interp, Xaxis=x_axis)
    
    # # p1_perp,p2_perp = uEXP.perp_p1p2(p1,p2, offset=0)
    # # pEXP.plot_line(xp, yp, yderiv ,p1_perp,p2_perp,title='yderiv',savefig=False, interp=interp, Xaxis=x_axis)    
    pEXP.plot_line(xp, yp, zderiv ,p1,p2,title='zderiv',savefig=False, interp=interp, Xaxis=x_axis)
    
    #%% ------------------------------- Plot the derivatives

    xx, yy, distance, profile = pEXP.plot_line(xp, yp, U ,p1,p2,title='U',savefig=False, interp=interp, Xaxis=x_axis)
    plt.plot(xx,profile)
    
    # Tau = np.gradient(profile,xx)
    Tau = np.gradient(np.log(profile),np.log(xx))
    # Tau = np.gradient(profile)
    plt.plot(xx,Tau)
