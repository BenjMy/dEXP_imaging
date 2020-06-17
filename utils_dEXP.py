# -*- coding: utf-8 -*-
"""
Miscellaneous utility functions.

@author: Benjamin
"""
import numpy as np 
import matplotlib.pyplot as plt

#%% Functions SPECIFIC FOR MISE-A-LA-MASSE data preparation to dEXP transformation 



#%% Geometrical correction/filtering

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


def perp_p1p2(p1,p2, offset=0):
    
    midX=(p1[0]+p2[0])/2
    midY=(p1[1]+p2[1])/2 
    
    plt.scatter(midX,midY,c='red')

    # midX=(p1[0]+offset+p2[0]+offset)/2
    # midY=(p1[1]+offset+p2[1]+offset)/2 
    
    # plt.scatter(midX,midY,c='green')


    new_p2 = [midX-p2[1]+p1[1], midY + p2[0]-p1[0]]
    new_p1 = [midX+p2[1]-p1[1], midY - p2[0]+p1[0]]

    p12x=[p1[0],p2[0]]
    p12y=[p1[1],p2[1]]
    
    p12x_new=[new_p1[0],new_p2[0]]
    p12y_new=[new_p1[1],new_p2[1]]
    
    plt.plot(p12x,p12y,c='red')
    plt.plot(p12x_new,p12y_new,c='green')
    plt.axis('square')
    
    
    return new_p1, new_p2