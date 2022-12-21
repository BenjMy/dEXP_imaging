<<<<<<< Updated upstream
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 08:01:49 2020

@author: Benjamin
"""

import pickle
import numpy as np

def load_MALM_sens3d(filename=None):

    file = open(filename, 'rb')
    u = pickle._Unpickler(file)
    u.encoding = 'latin1'
    data = u.load()
    

    # SimName='M' + 'SoilR' + str(data['SoilR']) + 'AnoR' + str(data['AnoR']) + 'Z' + str(data['HWD'][0]) + 'W' + str(data['HWD'][1]) +  'D' + str(data['HWD'][2])
    SimName = None

    # different version of the pickle file
    if "HWDL" in data:
        maxdepth = data['HWDL'][2] * 1.5
    elif "HWD" in data:
        maxdepth = data['HWD'][2] * 1.5
    elif "HWDLs" in data:
        maxdepth = data['HWDLs'][2] * 1.5

    shape = data['shape']
    p1 = data['p12'][0]
    p2 = data['p12'][1]
    
    if "HWDL" in data:
        xyzu = data['HWDL']
        xp, yp, z, U  = xyzu[:,0],  xyzu[:,1], xyzu[:,2],  xyzu[:,3]

        # return xp, yp, z, U, maxdepth, shape, p1, p2, SimName, data
    elif "HWD" in data:
        xyzu = data['uAz0_grid']
        xp, yp, z, U  = xyzu[:,0],  xyzu[:,1], xyzu[:,2],  xyzu[:,3]
        
    elif "HWDLs" in data:
        xyzu = data['XYU']
        xp, yp, z, U  = xyzu[:,0],  xyzu[:,1], xyzu[:,2],  xyzu[:,3]

        # return maxdepth, shape, p1, p2, SimName, data
    
    # add gaussian noise
    # ----------------
    relative_error_tmp = filename.split("Noise",1)[1]
    relative_error = int(relative_error_tmp.split(".",1)[0])/100
    noise_floor = 0.0
    
    print(relative_error)
    std = np.sqrt((relative_error * np.abs(U)) ** 2 + noise_floor ** 2)
    noise = std * np.random.randn(*U.shape)
    Unew = U + noise
            
            
    # noiseU = np.random.normal(np.mean(U),np.std(U)*noise/100, U.shape)*U 
    # noiseU = np.random.normal(np.mean(U),np.std(U)*noise/100, U.shape)*U 
    # Unew = U + np.random.randn(1)*U*noise/100
    
    
    # Unew = U + noiseU
    # np.mean(Unew)
    # np.mean(U)

    
    return xp, yp, z, Unew, maxdepth, shape, p1, p2, SimName, data


=======
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 08:01:49 2020

@author: Benjamin
"""

import pickle

def load_MALM_sens3d(filename=None):

    #file = open(filename, 'rb')
    #u = pickle._Unpickler(file)
    #u.encoding = 'latin1'
    #data = u.load()

    infile = open(filename,'rb')
    #data = pickle.load(infile,encoding='latin1')
    data = pickle.load(infile)
    infile.close()
    SimName = None

    # different version of the pickle file
    if "HWDL" in data:
        maxdepth = data['HWDL'][2] * 1.5
    elif "HWD" in data:
        maxdepth = data['HWD'][2] * 1.5
    elif "HWDLs" in data:
        maxdepth = data['HWDLs'][2] * 1.5

    shape = data['shape']
    p1 = data['p12'][0]
    p2 = data['p12'][1]
    
    
    if "HWDL" in data:
        xyzu = data['XYU']
        xp, yp, z, U  = xyzu[:,0],  xyzu[:,1], xyzu[:,2],  xyzu[:,3]

        # return xp, yp, z, U, maxdepth, shape, p1, p2, SimName, data
    elif "HWD" in data:
        try:
            xyzu = data['uAz0_grid']
        except:
            try:
                xyzu = data['XYU']
            except: print('not a valid input file')
                
        xp, yp, z, U  = xyzu[:,0],  xyzu[:,1], xyzu[:,2],  xyzu[:,3]
        
    elif "HWDLs" in data:
        xyzu = data['XYU']
        xp, yp, z, U  = xyzu[:,0],  xyzu[:,1], xyzu[:,2],  xyzu[:,3]

        # return maxdepth, shape, p1, p2, SimName, data
    return xp, yp, z, U, maxdepth, shape, p1, p2, SimName, data


>>>>>>> Stashed changes
