# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 20:09:26 2020

@author: Benjamin
"""
import numpy as np 
import matplotlib.pyplot as plt

def _point_against_line(xp,yp,p1,p2):
# https://stackoverflow.com/questions/45766534/finding-cross-product-to-find-points-above-below-a-line-in-matplotlib

    isabove = lambda p, p1,p2: np.cross(p-p1, p2-p1) < 0
    
    # p1 = np.array([1,1])
    # p2 = np.array([4,3])
    
    fig, (ax,ax2) = plt.subplots(ncols=2, sharex=True, sharey=True)

    # p = np.random.rand(10,2)*5
    p = np.vstack([xp,yp]).T
    print(p)
    
    ax2.plot([p1[0],p2[0]],[p1[1],p2[1]], marker="o", color="k")
    ax2.scatter(p[:,0],p[:,1], c=isabove(p,p1,p2), cmap="bwr", vmin=0, vmax=1)
    
    
    # ax.set_xlim(0,6)
    # ax.set_ylim(0,6)
    
    plt.show()

    return

def mirror(xp,yp, data, p1, p2, mirror='above'):
    
    # put zero over/under a certain line
    _point_against_line(xp,yp, p1, p2)
    
    # replicate values to mirror
    
    # put zero to 
    # data_sym = np.zeros()
    
    # for row in range(data.shape[0]):
    #     for col in range(data.shape[1]):
    #         data_sym = data()
            
    return


a = np.array([[1,2,3,4],[1,2,3,4],[1,2,3,4],[1,2,3,4]])
a[0,3] = a[3,0]