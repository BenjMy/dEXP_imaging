# -*- coding: utf-8 -*-
"""
Created on Thu May 28 08:01:49 2020

@author: Benjamin
"""

import pickle

def load_MALM_sens3d(filename=None):

    file = open(filename, 'rb')
    u = pickle._Unpickler(file)
    u.encoding = 'latin1'
    data = u.load()

    SimName='M' + 'SoilR' + str(data['SoilR']) + 'AnoR' + str(data['AnoR']) + 'Z' + str(data['HWD'][0]) + 'W' + str(data['HWD'][1]) +  'D' + str(data['HWD'][2])

    maxdepth = data['HWD'][2] * 1.5
    shape = data['shape']
    p1 = data['p12'][0]
    p2 = data['p12'][1]
    xyzu = data['XYU']
    xp, yp, z, U  = xyzu[:,0],  xyzu[:,1], xyzu[:,2],  xyzu[:,3]
    
    return xp, yp, z, U, maxdepth, shape, p1, p2, SimName, data

