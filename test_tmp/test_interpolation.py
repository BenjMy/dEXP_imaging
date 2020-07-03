# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 17:51:17 2020

@author: Benjamin
"""
import numpy as np
from scipy import interpolate

x = np.arange(0,10)
y = np.exp(-x/3.0)
f = interpolate.interp1d(x, y, fill_value='extrapolate')

p = 10
xnew = np.linspace(min(x)-p,max(x)+p,len(x)*2)
# len(x)
# len(xnew)

plt.scatter(x,y)
plt.plot(xnew,f(xnew))


# vp = scipy.interpolate.griddata((x, y), v, (xp, yp),
#                                 method=algorithm).ravel()
# if extrapolate and algorithm != 'nearest' and np.any(np.isnan(vp)):
#     fill_nans(x, y, v, xp, yp, vp)