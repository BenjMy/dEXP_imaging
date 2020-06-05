# -*- coding: utf-8 -*-
"""
Some functions to plot results from DEXP transformation
"""

import matplotlib.pyplot as plt
import numpy as np 
from fatiando import gridder, mesher, utils

from collections import OrderedDict
from scipy.optimize import curve_fit

# my own functions
import dEXP as dEXP

def plot_z(mesh):

    # ------------------------------- Plot the data
    image = mesh.props['density'].reshape(mesh.shape)
    # if scaled == 1:
    #     scale = 0.1*np.abs([image.min(), image.max()]).max()
    #     mins, maxs = [-scale, scale]
    # else:
    #     mins, maxs = [image.min(),image.max()]
    
    fig, ax = plt.subplots(nrows=1, ncols=5,figsize=(25,8))
    ## Get a horizontal section at the z
    for ss in range(5):
        ll = ss+1
        section = mesh.get_layer(ll*4+ll-1)
        zz=mesh.get_zs()
        ish= ax[ss].imshow(np.array(mesher.extract('density',section)).reshape(shape[0], shape[1]),
                   extent=mesh.bounds[0:4],origin='lower', cmap="viridis", vmin=mins, vmax=maxs)
        # square([y1, y2, x1, x2])
        ax[ss].set_title('z={:.2} m'.format(zz[ll*4+ll]), fontsize=20)
        ax[ss].set_aspect('equal')
        plt.show()
    #    plt.colorbar(ish, ax=ax[ss], pad=0)
        ax[ss].set_xlabel('x (m)', fontsize=20)
        ax[ss].set_ylabel('y (m)', fontsize=20)
    
    plt.tight_layout()
    # plt.suptitle(strname + '_Z_' + str(ZZ), fontsize=40)
    # plt.savefig(pathFig +strname + '_Z_' + str(ZZ) + '.png')

    return

def plot_xy(mesh,scaled=0,label=None, ax=None, markerMax=False):

    if label not in mesh.props:
        raise ValueError("mesh doesn't have a '%s' property." % (label))

    image = mesh.props[label].reshape(mesh.shape)

    if scaled == 1:
        scale = 0.1*np.abs([image.min(), image.max()]).max()
        mins, maxs = [-scale, scale]
    else:
        mins, maxs = [image.min(),image.max()]
    
    if ax == None:
        fig = plt.subplots()
        ax = plt.gca()

    
    x = mesh.get_xs()
    y = mesh.get_ys()
    z = mesh.get_zs()
    
    #ax = plt.subplot(1, 2, 1)
    ax.set_title('Model slice at y={} m'.format(y[len(y)//2]))
    ax.set_aspect('equal')
    
    cmap = ax.pcolormesh(x, z, image[:, :, mesh.shape[1]//2])
    # square([x1, x2, z1, z2])
    
    if markerMax == True:
        # search for the max
        ind = np.unravel_index(np.argmax(image[:, :, mesh.shape[1]//2], axis=None), 
                               image[:, :, mesh.shape[1]//2].shape)
        image[:, :, mesh.shape[1]//2][ind]
        zmax = z[ind[0]]
        #ymax = y[ind[1]]
        xmax = -x[ind[1]]
        
        ax.scatter(xmax,zmax, s=70, c='r', marker='v')
    



    ax.set_ylim(z.max(), z.min())
    ax.set_xlim(x.min(), x.max())
    ax.set_xlabel('x (m)')
    ax.set_ylabel('z (m)')
    ax.set_title(label)
    
    if 'upwc' in label:
        plt.gca().invert_yaxis()

    #ax = plt.subplot(1, 2, 2)
    #ax.set_title('Model slice at x={} m'.format(x[len(x)//2]))
    #ax.set_aspect('equal')
    #ax.pcolormesh(y, z, image[:, mesh.shape[2]//2, :], cmap="cubehelix",
    #              vmin=mins, vmax=maxs)
    #square([y1, y2, z1, z2])
    #ax.set_ylim(z.max(), z.min())
    #ax.set_xlim(y.min(), y.max())
    #ax.set_xlabel('y (km)')
    #ax.set_ylabel('z (km)')
    
    #plt.tight_layout()
    plt.show()
    
    #plt.tight_layout()
    # plt.suptitle(strname + '_xy_' + str(ZZ), fontsize=15)
    # plt.savefig(pathFig +strname + '_xy_' + str(ZZ) + '.png')

    return plt, cmap

# def plot_line_mesh(mesh, lnb= 0, p1,p2,ax=None):
    
#     # Extract a profile between points 1 and 2
#     xx, yy, distance, profile = gridder.profile(x, y, data, p1, p2, 1000)
    
#     # Plot the profile and the original map data
#     plt.figure()
#     ax = ax or plt.gca()
#     ax = ax or plt.gca()
#     plt.subplot(2, 1, 1)
#     # plt.title(strname + '_data' + str(ZZ), fontsize=15)
#     plt.plot(distance, profile, '.k')
#     plt.xlim(distance.min(), distance.max())
#     plt.grid()
#     plt.subplot(2, 1, 2)
#     plt.title("Original data")
#     plt.plot(xx, yy, '-k', label='Profile', linewidth=2)
#     scale = np.abs([data.min(), data.max()]).max()
#     plt.tricontourf(x, y, data, 50, cmap='RdBu_r', vmin=-scale, vmax=scale)
#     plt.colorbar(orientation='horizontal', aspect=50)
#     plt.legend(loc='lower right')
#     plt.tight_layout()
#     plt.show()
#     #plt.suptitle(strname + '_ztop' + str(za) +'_zbot'+ str(zb), fontsize=15)
#     # plt.savefig(pathFig+ strname + '_data' + str(ZZ) + '.png')

#     return ax

def plot_line(x,y,data,p1,p2,ax=None,interp=True,**kwargs):
    
    # Extract a profile between points 1 and 2
    if interp == False:
        xx, yy, distance, profile = dEXP.profile_noInter(x, y, data, p1, p2, 1000)
    else:
        xx, yy, distance, profile = gridder.profile(x, y, data, p1, p2, 1000)
    
    # Plot the profile and the original map data
    plt.figure()
    ax = ax or plt.gca()
    plt.subplot(2, 1, 1)
    # plt.title(strname + '_data' + str(ZZ), fontsize=15)
    plt.plot(xx, profile, '.k')
    plt.xlim(xx.min(), xx.max())
    plt.grid()
    plt.subplot(2, 1, 2)
    plt.title("Original data")
    plt.plot(xx, yy, '-k', label='Profile', linewidth=2)
    # scale = np.abs([data.min(), data.max()]).max()
    plt.tricontourf(x, y, data, 50, cmap='RdBu_r', vmin=min(data), vmax=max(data))
    plt.colorbar(orientation='horizontal', aspect=50)
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.axis('square')

    for key, value in kwargs.items():
        # print("{0} = {1}".format(key, value))

        if key == 'title':
            plt.title(value, fontsize=15)        
        if key == 'suptitle':
            plt.suptitle(value, fontsize=15)
        if key == 'savefig':
            if value == True:
            # plt.savefig(pathFig+ strname + '_data' + str(ZZ) + '.png')
                plt.savefig('fig2d.png', r=400)
        
    plt.show()

    return ax


def plot_ridges_harmonic(RI=None,RII=None,RIII=None,ax=None,label=False,**kwargs):
    """
    Plot ridges in the harmonic domain

    Parameters:

    * a
        Text here

    Returns:

    * BB : 
        Text here

    """
    # depths = R[:,:][i][1]

    if ax == None:
        fig = plt.subplots()
        ax = plt.gca()
    
    bbox = {'fc': '0.8', 'pad': 0}

    if RI is not None:
        for i, cl in enumerate([datacol for datacol in RI.columns if datacol != 'elevation']):
            RI.plot(x=cl, y='elevation', kind="scatter", ax=ax,label='Ridge I',c='r')
            if label==True:
                ax.text(RI[cl].iloc[0], max(RI['elevation']), 'RI'+str(i),rotation=90, bbox= bbox)
    if RII is not None:
        for i, cl in enumerate([datacol for datacol in RII.columns if datacol != 'elevation']):
            RII.plot(x=cl, y='elevation', kind="scatter", ax=ax,label='Ridge II',c='b')
            if label==True:
                ax.text(RII[cl].iloc[0], max(RII['elevation']), 'RII'+str(i),rotation=90, bbox= bbox)
    if RIII is not None:
        for i, cl in enumerate([datacol for datacol in RIII.columns if datacol != 'elevation']):
            RIII.plot(x=cl, y='elevation', kind="scatter", ax=ax,label='Ridge III',c='g')
            if label==True:
                ax.text(RIII[cl].iloc[0], max(RIII['elevation']), 'RIII'+str(i),rotation=90, bbox= bbox)
# plt.text(60, .025, r'$\mu=100,\ \sigma=15$')

            
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())


    return ax


def plot_ridges_sources(df_fit, ax=None, ridge_type=None, ridge_nb=None, z_max_source=None):
    """
    Plot ridges in the source domain and observe intersection point

    Parameters:
    * df_fit
        dataframe fit
    Returns:

    * BB : 
        Text here

    """
    if ax == None:
        fig = plt.subplots()
        ax = plt.gca()
    plt.rcParams['font.size'] = 15

    if ridge_type is None:
        ridge_type = np.arange(0,len(df_fit))

    ridge_nb_n = []
    if ridge_nb is None:
        for r_type in enumerate(ridge_type):
            ridge_nb_n_tmp = np.arange(0,df_fit[r_type[1]].shape[1])
            ridge_nb_n.append(ridge_nb_n_tmp)
    else:
        ridge_nb_n = ridge_nb
    
    # print(ridge_nb_n)
            
    for r_type in enumerate(ridge_type):
        
        for r_nb in enumerate(ridge_nb_n[r_type[0]]):
            print(r_type[1],r_nb[1])
            
            name_col = df_fit[r_type[1]].columns[r_nb[1]][0]
            # print(name_col)
            # print(r_type[1])
            ax.plot(df_fit[r_type[1]][name_col]['x'],
                    df_fit[r_type[1]][name_col]['y'], 'g--')


            # plt.plot(df_fit[r_type[1]][r_nb[1]])


            # ax.plot(fit[r_nb[0]][:,0], fit[r_nb[0]][:,1], 'g--')
            # # ax.scatter(points[i[0]][:,0], points[i[0]][:,1],marker='*')
            
            ymin, ymax = ax.get_ylim()
    
            if z_max_source is None:
                ax.set_ylim([-ymax*3,ymax])   
            else:
                ax.set_ylim([z_max_source,ymax])   
    
            ax.set_xlabel('x (m)', size=20)
            # ax.set_ylabel('depth (m)')
            # plt.title(r'$\frac{\partial log(f)}{\partial log(z)}$', size=20)
            plt.grid()
            # plt.legend()
            # the text bounding box
            plt.ylabel(" ")
            bbox = {'fc': '0.8', 'pad': 0}
            ax.text(-0.15, 0.3, 'depth', transform=ax.transAxes, fontsize=14, rotation=90, bbox= bbox)
            ax.text(-0.15, 0.8, 'altitude', transform=ax.transAxes, fontsize=14, rotation=90, bbox= bbox)

    return ax
    

def plot_scalFUN(points, fit, ax=None, z0=None):
    """
    Plot scalfun function analysis

    Parameters:

    * a
        Text here

    Returns:

    * BB : 
        Text here

    """

    if ax == None:
        fig = plt.subplots()
        ax = plt.gca()

    if ~isinstance(fit, list):
    #    fit = np.array([fit])
    #    points = np.array([points])
        z0 = np.array([z0])

    # # if len(z0)<2:
    # #     z0 = [z0]
    # fit = [fit]
    
    for z in enumerate(z0):
        ax.plot(fit[z[0]][:,0], fit[z[0]][:,1], '--',
                 label='fit_z0=' + str(z[1]), color='red')
        ax.scatter(points[z[0]][:,0], points[z[0]][:,1],marker='*')
        print(points[z[0]][:,0])
        ax.set_xlim([0,max(points[z[0]][:,0])])
    ax.set_ylim([-5,5])        
    ax.set_xlabel('q (m)', size=20)
    ax.set_ylabel('$\\tau_{f}$', size=20)
    # plt.title(r'$\frac{\partial log(f)}{\partial log(z)}$', size=20)
    plt.grid()
    plt.legend()
    
    return ax


# def plotDEXP(x,y,depths,list_up, list_indmax):
    
#     for uu in enumerate(list_up):
#         xx, yy, distance, p_up_f = gridder.profile(x, y, uu[1], p1, p2, 1000)

#     X, Y = np.meshgrid(xx, depths)
    
#     # Plot the profile and the original map data
#     plt.figure()
#     plt.subplot(3, 3, 1)
#     plt.plot(xx, profile, '.k')
#     plt.xlim(xx.min(), xx.max())
#     plt.xlabel('position (m)')
#     plt.ylabel('Field u')
#     plt.grid()
#     #
#     plt.subplot(3, 3, 2)
#     d1 = transform.derivz(xp, yp, gz, shape,order=1)
#     xx, yy, distance, p_d1 = gridder.profile(xp, yp, d1, p1, p2, 1000)
#     plt.plot(xx, p_d1, '.k')
#     plt.xlim(xx.min(), xx.max())
#     plt.xlabel('position (m)')
#     plt.ylabel('1st vertical derivative')
#     plt.grid()
    
#     plt.subplot(3, 3, 3)
#     d2 = transform.derivz(xp, yp, gz, shape,order=2)
#     xx, yy, distance, p_d2 = gridder.profile(xp, yp, d2, p1, p2, 1000)
#     plt.plot(xx, p_d2, '.k')
#     plt.xlim(xx.min(), xx.max())
#     plt.xlabel('position (m)')
#     plt.ylabel('2nd vertical derivative')
#     plt.grid()
    
#     plt.subplot(3, 3, 4)
#     plt.contourf(X, Y, up_f_sec)
#     plt.xlabel('position (m)')
#     plt.ylabel('Field u continuated (altitude)')
#     plt.grid()
    
    
#     plt.subplot(3, 3, 5)
#     plt.contourf(X, Y, up_f_d1_sec)
#     plt.xlabel('position (m)')
#     plt.ylabel('\delta u_c / \delta z')
#     plt.grid()
    
#     plt.subplot(3, 3, 6)
#     plt.contourf(X, Y, up_f_d2_sec)
#     plt.xlabel('position (m)')
#     plt.ylabel('\delta^2 u_c / \delta^2 z')
#     plt.grid()
    
#     plt.subplot(3, 3, 7)
#     plt.contourf(X, Y, up_f_w_sec)
#     plt.scatter(X[list_indmax[3]],Y[list_indmax[3]], s=70, c='w', marker='v')
#     plt.xlabel('position (m)')
#     plt.ylabel('dEXP(u_c)')
#     plt.grid()
#     x1, x2, y1, y2, z1, z2 = np.array(model[0].get_bounds())
#     square([x1, x2,z1, z2])
#     plt.gca().invert_yaxis()
    
#     plt.subplot(3, 3, 8)
#     plt.contourf(X, Y, up_f_d1_w_sec)
#     plt.scatter(X[list_indmax[4]],Y[list_indmax[4]], s=70, c='w', marker='v')
#     plt.xlabel('position (m)')
#     plt.ylabel('dEXP(\delta u_c / \delta z)')
#     plt.grid()
#     x1, x2, y1, y2, z1, z2 = np.array(model[0].get_bounds())
#     square([x1, x2,z1, z2])
#     plt.gca().invert_yaxis()
    
#     plt.subplot(3, 3, 9)
#     plt.contourf(X, Y, up_f_d1_w_sec)
#     plt.scatter(X[list_indmax[5]],Y[list_indmax[5]], s=70, c='w', marker='v')
#     plt.xlabel('position (m)')
#     plt.ylabel('dEXP(\delta u_c / \delta z)')
#     plt.grid()
#     x1, x2, y1, y2, z1, z2 = np.array(model[0].get_bounds())
#     square([x1, x2,z1, z2])
#     if len(model)>1:
#         x1, x2, y1, y2, z1, z2 = np.array(model[1].get_bounds())
#         square([x1, x2,z1, z2])
#     plt.gca().invert_yaxis()
