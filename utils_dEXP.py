# -*- coding: utf-8 -*-
"""
Miscellaneous utility functions.

@author: Benjamin
"""
import numpy as np 
import matplotlib.pyplot as plt
import math 
from numpy import inf

#%% Functions SPECIFIC FOR MISE-A-LA-MASSE data preparation to dEXP transformation 

def cor_field_B(x,y,z,u,B,rho=100):
    """
    Calculates the potential field (electric) produced by a current injection in B (return electrode) for a
    given homogeneous electrical resistivity rho
    """
    I = 1 # injected current (A)
    num = rho*I
    dist = np.sqrt((x-B[0])**2 + (y-B[1])**2 + (z-B[2])**2) # distance between B and potential electrodes
    den = 2*math.pi*dist
    u_B = num/den
    
    print(u_B)
    u_B[u_B == inf] = 0
    
    u_cor = u - u_B # correct total measured potential from influence of B

    plt.figure() # figsize=(20,10)
    plt.subplot(2,2,1)
    plt.tricontourf(x, y, dist, 50, cmap='RdBu_r')
    cbar = plt.colorbar(orientation='vertical', aspect=50)
    cbar.set_label('Distance from B (m)')
    # plt.scatter(B[0],B[1], marker='v', color='red')
    # plt.tight_layout()
    plt.axis('square')

    plt.subplot(2,2,2)
    plt.tricontourf(x, y, u_B, 50, cmap='RdBu_r')
    cbar = plt.colorbar(orientation='vertical', aspect=50)
    cbar.set_label('$u_{B}$ (V)')
    # plt.scatter(B[0],B[1], marker='v', color='red')
    # plt.tight_layout()
    plt.axis('square')
    
    plt.subplot(2,2,3)
    plt.tricontourf(x, y, u, 50, cmap='RdBu_r')
    cbar = plt.colorbar(orientation='vertical', aspect=50)
    cbar.set_label('$u = U_{T} $(V)')
    # plt.scatter(B[0],B[1], marker='v', color='red')
    # plt.tight_layout()
    plt.axis('square')
    
    plt.subplot(2,2,4)
    plt.tricontourf(x, y, u_cor, 50, cmap='RdBu_r')
    cbar = plt.colorbar(orientation='vertical', aspect=50)
    cbar.set_label('$u = U_{T} - u_{B}$ (V)')
    # plt.scatter(B[0],B[1], marker='v', color='red')
    # plt.tight_layout()
    plt.axis('square')

    return u_cor   


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