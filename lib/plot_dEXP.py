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
import lib.dEXP as dEXP
from mpl_axes_aligner import align


def plot_field(x, y, field, shape, **kwargs):

    Vminmax = None
    for key, value in kwargs.items():
        if key == "Vminmax":
            Vminmax = value
            mins = value[0]
            maxs = value[1]

    # Make maps of all fields calculated
    fig = plt.figure()
    ax = plt.gca()
    plt.rcParams["font.size"] = 10
    X, Y = x.reshape(shape), y.reshape(shape)

    if Vminmax is None:
        scale = np.abs([field.min(), field.max()]).max()
        mins = -scale
        maxs = scale

    # ax.set_title(field_name)
    plot = ax.pcolormesh(
        X, Y, field.reshape(shape), cmap="RdBu_r", vmin=mins, vmax=maxs, rasterized=True
    )
    cbar = plt.colorbar(plot, ax=ax)
    cbar.set_label("voltage (V)")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    plt.tight_layout(pad=0.5)
    # plt.show()

    return ax, plt


def plot_z(mesh):
    """
    Get horizontal sections at z levels and plot it

    Parameters:

    * mesh
        Fatiando mesh type

    """
    image = mesh.props["density"].reshape(mesh.shape)
    # if scaled == 1:
    #     scale = 0.1*np.abs([image.min(), image.max()]).max()
    #     mins, maxs = [-scale, scale]
    # else:
    #     mins, maxs = [image.min(),image.max()]

    fig, ax = plt.subplots(nrows=1, ncols=5, figsize=(25, 8))
    ## Get a horizontal section at the z
    for ss in range(5):
        ll = ss + 1
        section = mesh.get_layer(ll * 4 + ll - 1)
        zz = mesh.get_zs()
        ish = ax[ss].imshow(
            np.array(mesher.extract("density", section)).reshape(shape[0], shape[1]),
            extent=mesh.bounds[0:4],
            origin="lower",
            cmap="viridis",
            vmin=mins,
            vmax=maxs,
        )
        # square([y1, y2, x1, x2])
        ax[ss].set_title("z={:.2} m".format(zz[ll * 4 + ll]), fontsize=20)
        ax[ss].set_aspect("equal")
        plt.show()
        #    plt.colorbar(ish, ax=ax[ss], pad=0)
        ax[ss].set_xlabel("x (m)", fontsize=20)
        ax[ss].set_ylabel("y (m)", fontsize=20)

    plt.tight_layout()
    # plt.suptitle(strname + '_Z_' + str(ZZ), fontsize=40)
    # plt.savefig(pathFig +strname + '_Z_' + str(ZZ) + '.png')

    return


def plot_xy(mesh, scaled=0, label=None, ax=None, markerMax=False, **kwargs):
    """
    Get vertical xy section of the mesh at max/2 and plot it

    Parameters:

    * mesh
        Fatiando mesh type
    * scaled
        Txt here
    * label
        Txt here
    * ax
        Txt here
    * markerMax
        Txt here
    """

    x = mesh.get_xs()
    y = mesh.get_ys()
    z = mesh.get_zs()

    Xaxis = "x"
    p1p2 = None
    dEXP_bool = False
    SI = None
    q_ratio = None
    Vminmax = None
    aspect_equal = False
    for key, value in kwargs.items():
        if key == "Xaxis":
            Xaxis = value
        if key == "p1p2":
            p1p2 = value
        if key == "markerMax":
            markerMax = value
        if key == "SI":
            SI = value
        if key == "qratio":
            q_ratio = value
        if key == "Vminmax":
            Vminmax = value
            vmin = Vminmax[0]
            vmax = Vminmax[1]
        if key == "aspect_equal":
            aspect_equal = True

    if label not in mesh.props:
        raise ValueError("mesh doesn't have a '%s' property." % (label))

    image = mesh.props[label].reshape(mesh.shape)

    if Vminmax is not None:
        mins, maxs = [vmin, vmax]
    else:
        if scaled == 1:
            scale = 1e-5 * np.abs([image.min(), image.max()]).max()
            mins, maxs = [-scale, scale]
        else:
            mins, maxs = [image.min(), image.max()]
            # mins, maxs = 0,0.01

    if ax == None:
        fig = plt.subplots(figsize=(15, 3))
        ax = plt.gca()

    if p1p2 is not None:
        p_xaxis = []
        for i in p1p2[0]:
            if i in p1p2[1]:
                p_xaxis.append(i)
            else:
                # print(p_xaxis)
                # print(p1p2)
                print("need to rotate first?")

        if len(p_xaxis) > 1:
            print("len must be <2")

        if Xaxis is "x":
            # p_xaxis= p1p2[0][0]
            # print("paxis=" + str(p_xaxis))
            id_p_xaxis = (np.abs(x - p_xaxis)).argmin()
            ax.set_title(
                "slice at x={} m".format(int(y[id_p_xaxis])), fontsize=10, loc="left"
            )
            cmap = ax.pcolormesh(y, z, image[:, id_p_xaxis, :], vmin=mins, vmax=maxs)
            # cmap = ax.pcolormesh(x, z, image[:, :, id_p_xaxis])
            # ax.set_xlim(y.min(), y.max())
            ax.set_xlabel("y (m)")
        else:
            # p_xaxis= p1p2[0][1]
            id_p_xaxis = (np.abs(y - p_xaxis)).argmin()
            ax.set_title("slice at y={} m".format(int(x[id_p_xaxis])), fontsize=10)
            cmap = ax.pcolormesh(x, z, image[:, :, id_p_xaxis], vmin=mins, vmax=maxs)
            # cmap = ax.pcolormesh(y, z, image[:, id_p_xaxis, :])
            # ax.set_xlim(x.min(), x.max())
            ax.set_xlabel("x (m)")
    else:
        cmap = ax.pcolormesh(
            x, z, image[:, :, mesh.shape[1] // 2], vmin=mins, vmax=maxs
        )
        ax.set_xlim(x.min(), x.max())
        ax.set_xlabel("x (m)")
    # ax.set_aspect('equal')
    # plt.show()

    # if markerMax == True:
    #     # search for the max
    #     ind = np.unravel_index(np.argmax(image[:, :, id_p_xaxis], axis=None),
    #                            image[:, :, mesh.shape[1]//2].shape)
    #     image[:, :, mesh.shape[1]//2][ind]
    #     zmax = z[ind[0]]
    #     #ymax = y[ind[1]]
    #     xmax = -x[ind[1]]

    #     ax.scatter(xmax,zmax, s=70, c='r', marker='v')

    ax.set_ylim(z.max(), z.min())

    ax.set_ylabel("z (m)")
    ax.set_ylabel("depth\n(m)")
    # ax.set_title(label)

    if "upwc" in label:
        plt.gca().invert_yaxis()

    if markerMax == True:
        if p1p2 is not None:
            if Xaxis is "x":
                x_axis = x
                id_p_xaxis = (np.abs(x_axis - p_xaxis)).argmin()
                ind = np.unravel_index(
                    np.argmax(image[:, id_p_xaxis, :], axis=None),
                    image[:, id_p_xaxis, :].shape,
                )
            else:
                # search for the max
                x_axis = y
                id_p_xaxis = (np.abs(x_axis - p_xaxis)).argmin()
                ind = np.unravel_index(
                    np.argmax(image[:, :, id_p_xaxis], axis=None),
                    image[:, :, id_p_xaxis].shape,
                )
            z_exp = z[ind[0]]
            if Xaxis is "x":
                x_axis_exp = y[ind[1]]
            elif Xaxis is "y":
                x_axis_exp = x[ind[1]]

            # x_axis_exp = x_axis[ind[1]]
            print("Markermax_z=" + str(z_exp))
            print("Markermax_x=" + str(x_axis_exp))
            ax.scatter(x_axis_exp, z_exp, s=70, c="r", marker=".")

        if SI is not None:
            ax.legend(["SI=" + str(np.round(SI, 1))])
        if q_ratio is not None:
            ax.legend(
                ["q-ratio=" + str(q_ratio)]
            )  # bbox_to_anchor=(1.05, 0.05), loc='lower left')

    # ax = plt.subplot(1, 2, 2)
    # ax.set_title('Model slice at x={} m'.format(x[len(x)//2]))
    if aspect_equal:
        ax.set_aspect("equal")
    # ax.pcolormesh(y, z, image[:, mesh.shape[2]//2, :], cmap="cubehelix",
    #              vmin=mins, vmax=maxs)
    # square([y1, y2, z1, z2])
    # ax.set_ylim(z.max(), z.min())
    # ax.set_xlim(y.min(), y.max())
    # ax.set_xlabel('y (km)')
    # ax.set_ylabel('z (km)')

    # plt.tight_layout()
    # plt.show()

    # plt.tight_layout()
    # plt.suptitle(strname + '_xy_' + str(ZZ), fontsize=15)
    # plt.savefig(pathFig +strname + '_xy_' + str(ZZ) + '.png')
    # ax.set_aspect(5)

    return plt, cmap


def slice_mesh(x, y, mesh, label, p1, p2, ax=None, interp=True, **kwargs):
    """
    Get vertical xy section of the mesh at max/2 and plot it

    Parameters:

    * mesh
        Fatiando mesh type
    * scaled
        Txt here
    * label
        Txt here
        Txt here
    * markerMax
        Txt here
    """

    if ax == None:
        fig = plt.subplots()
        ax = plt.gca()

    Xaxis = []
    for key, value in kwargs.items():
        if key == "Xaxis":
            Xaxis = value

    image = mesh.props[label]
    toslice = np.reshape(image, [mesh.shape[0], mesh.shape[1] * mesh.shape[2]])
    depths = mesh.get_zs()[:-1]
    x_resolution = 1000

    Z = []
    X = []
    Y = []
    DIST = []
    IMG = []
    z = 0 * np.ones_like(x)
    interp = True
    for i, depth in enumerate(depths - z[0]):  # Loop for RII extremas
        u_layer = toslice[i, :]  # analysing extrema layers by layers from top to bottom

        if interp == True:
            xx, yy, distance, u_layer_p1p2 = gridder.profile(
                x, y, u_layer, p1, p2, x_resolution
            )
        else:
            xx, yy, distance, u_layer_p1p2, u_layer_p1p2_smooth = profile_noInter(
                x, y, u_layer, p1, p2, x_resolution
            )

        zz = depth * np.ones(len(u_layer_p1p2))
        Z.append(zz)
        X.append(xx)
        Y.append(yy)
        DIST.append(distance)
        IMG.append(u_layer_p1p2)

    if "dist" in Xaxis:
        xaxis = DIST
    else:
        xaxis = X

    cmap = ax.pcolormesh(np.array(xaxis), np.array(Z), np.array(IMG))

    # if 'upwc' in label:
    #     plt.gca().invert_yaxis()

    # ax.set_ylim(np.max(Z), np.min(Z))
    # ax.set_xlim(np.max(xaxis), np.min(xaxis))
    ax.set_xlabel("x (m)")
    ax.set_ylabel("z (m)")
    ax.set_title(label)
    ax.set_aspect("equal")

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


def plot_line(x, y, data, p1, p2, ax=None, interp=True, **kwargs):
    """
    

    Parameters
    ----------
    x : TYPE
        x coords of the 2d map.
    y : TYPE
        y coords of the 2d map..
    data : TYPE
        2d map dataset.
    p1 : list
        initial point coordinate (x,y) of the extracted profile.
    p2 : list
        final point coordinate (x,y) of the extracted profile.
    **kwargs 

    Returns
    -------
    xx : np.array
        interpolated x.
    yy : np.array
        interpolated y.
    distance : np.array
        profile distance.
    profile : TYPE
        data values along the profile.
    """

    Xaxis = []
    Vminmax = []
    smooth_type = False
    showfig = False
    vmin = min(data)
    vmax = max(data)
    limx = None
    limy = None
    x_resolution = None

    for key, value in kwargs.items():
        if key == "Xaxis":
            Xaxis = value
        if key == "smooth":
            smooth_type = value
        if key == "showfig":
            showfig = value
        if key == "Vminmax":
            Vminmax = value
            vmin = Vminmax[0]
            vmax = Vminmax[1]
        if key == "limx":
            limx = value
        if key == "limy":
            limy = value
        if key == "x_resolution":
            x_resolution = value
    # Extract a profile between points 1 and 2
    if interp == False:
        xx, yy, distance, profile, vp_smooth_dict = dEXP.profile_noInter(
            x, y, data, p1, p2, size=x_resolution, showfig=showfig
        )
        # xx, yy, distance, profile, profile_smooth = dEXP.profile_extra(x, y, data, p1, p2, 1000)
    else:
        xx, yy, distance, profile = gridder.profile(x, y, data, p1, p2, 1000, "nearest")

    if Xaxis == "dist":
        xaxis = distance
    elif Xaxis == "y":
        xaxis = xx
    else:
        xaxis = yy

    if smooth_type is not False:
        if smooth_type == True:
            profile = vp_smooth_dict["Lowpass"]
        else:
            profile = vp_smooth_dict[smooth_type]

    # Plot the profile and the original map data
    plt.figure()
    # ax = ax or plt.gca()
    plt.subplot(1, 2, 1)
    # plt.title(strname + '_data' + str(ZZ), fontsize=15)
    plt.plot(xaxis, profile, ".k")
    # plt.xlim(xaxis.min(), xaxis.max())
    plt.grid()

    if Xaxis == "dist":
        plt.xlabel("distance (m)")
    elif Xaxis == "y":
        plt.xlabel("x (m)")
    else:
        plt.xlabel("y (m)")

    plt.ylabel("voltage (V)")

    plt.subplot(1, 2, 2)
    plt.title("Original data")
    plt.plot(xx, yy, "-k", label="Profile", linewidth=2)
    # scale = np.abs([data.min(), data.max()]).max()
    # plt.tricontourf(x, y, data, 50, cmap='RdBu_r', vmin=min(data), vmax=max(data))
    plt.scatter(x, y, c=data, cmap="RdBu_r", vmin=vmin, vmax=vmax)
    plt.colorbar(orientation="vertical", aspect=50)
    plt.legend(loc="lower right")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")

    plt.axis("square")

    if limx is not None:
        plt.xlim(limx[0], limx[1])
    if limx is not None:
        plt.ylim(limy[0], limy[1])

        # plt.tight_layout()

    for key, value in kwargs.items():
        # print("{0} = {1}".format(key, value))

        if key == "title":
            plt.title(value, fontsize=15)
        if key == "suptitle":
            plt.suptitle(value, fontsize=15)
        if key == "savefig":
            if value == True:
                # plt.savefig(pathFig+ strname + '_data' + str(ZZ) + '.png')
                plt.savefig("fig2d.png", r=400)

    plt.show()

    return xx, yy, distance, profile, ax, plt


def plot_ridges_harmonic(
    RI=None, RII=None, RIII=None, ax=None, label=False, legend=True, **kwargs
):
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

    bbox = {"fc": "0.8", "pad": 0}

    if RI is not None:
        for i, cl in enumerate(
            [datacol for datacol in RI.columns if datacol != "elevation"]
        ):
            RI.plot(x=cl, y="elevation", kind="scatter", ax=ax, label="Ridge I", c="r")
            if label == True:
                ax.text(
                    RI[cl].iloc[0],
                    max(RI["elevation"]),
                    "RI" + str(i),
                    rotation=90,
                    bbox=bbox,
                    fontsize=8,
                )
    if RII is not None:
        for i, cl in enumerate(
            [datacol for datacol in RII.columns if datacol != "elevation"]
        ):
            RII.plot(
                x=cl, y="elevation", kind="scatter", ax=ax, label="Ridge II", c="b"
            )
            if label == True:
                ax.text(
                    RII[cl].iloc[0],
                    max(RII["elevation"]),
                    "RII" + str(i),
                    rotation=90,
                    bbox=bbox,
                    fontsize=8,
                )
    if RIII is not None:
        for i, cl in enumerate(
            [datacol for datacol in RIII.columns if datacol != "elevation"]
        ):
            RIII.plot(
                x=cl, y="elevation", kind="scatter", ax=ax, label="Ridge III", c="g"
            )
            if label == True:
                ax.text(
                    RIII[cl].iloc[0],
                    max(RIII["elevation"]),
                    "RIII" + str(i),
                    rotation=90,
                    bbox=bbox,
                    fontsize=8,
                )
    # plt.text(60, .025, r'$\mu=100,\ \sigma=15$')

    if legend == True:
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys(), fontsize=11)
    else:
        ax.get_legend().remove()

    ax.set_ylabel("Elevation (m)", loc="top")

    return ax


def plot_ridges_sources(
    df_fit, ax=None, ridge_type=None, ridge_nb=None, z_max_source=None, **kwargs
):
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
        ax2 = plt.gca()
    else:
        # print(ax)
        ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis

    plt.rcParams["font.size"] = 15

    if ridge_type is None:
        ridge_type = np.arange(0, len(df_fit))

    ridge_nb_n = []
    if ridge_nb is None:
        for r_type in enumerate(ridge_type):
            if len(np.array(df_fit[r_type[1]])) > 1:
                ridge_nb_n_tmp = np.arange(0, df_fit[r_type[1]].shape[1])
                ridge_nb_n.append(ridge_nb_n_tmp)
            else:
                ridge_nb_n.append(None)
    else:
        ridge_nb_n = ridge_nb
    # print(ridge_nb_n)

    for r_type in enumerate(ridge_type):

        if ridge_nb_n[r_type[1]] is not None:
            # print(ridge_nb_n[r_type[1]])
            for r_nb in enumerate(ridge_nb_n[r_type[0]]):
                # print(r_type[1],r_nb[1])

                name_col = df_fit[r_type[1]].columns[r_nb[1]][0]
                # print(name_col)
                # print(r_type[1])
                ax2.plot(
                    df_fit[r_type[1]][name_col]["x"],
                    df_fit[r_type[1]][name_col]["y"],
                    "k--",
                )

                # print(df_fit[1]['R1 Vert.EX_xpos1']['x'])
                # # plt.plot(df_fit[r_type[1]][r_nb[1]])
                # ax.plot(df_fit[1]['R1 Vert.EX_xpos1']['x'],
                #         df_fit[1]['R1 Vert.EX_xpos1']['y'], 'g--')

                # ax.plot(fit[r_nb[0]][:,0], fit[r_nb[0]][:,1], 'g--')
                # # ax.scatter(points[i[0]][:,0], points[i[0]][:,1],marker='*')


        ymin, ymax = ax.get_ylim()
        # print(z_max_source)

        if z_max_source is None:
            ax2.set_ylim([-ymax * 3, ymax])
        else:
            ax2.set_ylim([z_max_source, ymax])

        for key, value in kwargs.items():
            if key == "xmin":
                x_min = value
                ax2.set_xlim([x_min, None])
            if key == "xmax":
                x_max = value
                ax2.set_xlim([None, x_max])


        color = "tab:red"
        ax2.set_ylabel(
            "depth\n(m)", color=color, loc="bottom"
        )  # we already handled the x-label with ax1
        ax2.set_xlabel("y (m)")  # we already handled the x-label with ax1
        # ax2.xaxis.set_label_coords(0,-3000)
        ax2.tick_params(axis="y", labelcolor=color)
        # plt.title(r'$\frac{\partial log(f)}{\partial log(z)}$', size=20)
        plt.grid()
        # plt.legend()
        # the text bounding box
        # plt.ylabel("")

        # plt.axis('square')

        # bbox = {'fc': '0.8', 'pad': 0}
        # ax.text(-0.15, 0.3, 'depth', transform=ax.transAxes, fontsize=14, rotation=90, bbox= bbox)
        # ax.text(-0.15, 0.8, 'altitude', transform=ax.transAxes, fontsize=14, rotation=90, bbox= bbox)

        ax2.spines["right"].set_position(("axes", 0.6))
        ax2.spines['right'].set_color('red')

        labels_ax1 = ax.get_yticks() 
        labels_ax1= labels_ax1[labels_ax1>0]
        
        labels_ax2 = ax2.get_yticks() 
        labels_ax2= labels_ax2[labels_ax2<0]
        
        ax.set_yticks(labels_ax1)
        ax2.set_yticks(labels_ax2)
        ax2.spines["right"].set_visible(False)

        # Adjust the plotting range of two y axes
        org1 = 0.0  # Origin of first axis
        org2 = 0.0  # Origin of second axis
        pos = 0.5  # Position the two origins are aligned
        align.yaxes(ax, org1, ax2, org2, pos)


    return ax2


def plot_scalFUN(points, fit, ax=None, z0=None, label=None):
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
        ax.plot(
            fit[z[0]][:, 0],
            fit[z[0]][:, 1],
            "--",
            label="fit_z0=" + str(z[1]),
            color="red",
        )
        ax.scatter(
            points[z[0]][:, 0],
            points[z[0]][:, 1],
            marker="*",
            label="Ridge id:" + str(label),
        )
        # print(points[z[0]][:,0])
        ax.set_xlim([0, max(points[z[0]][:, 0])])
    ax.set_ylim([-5, 5])
    ax.set_xlabel("q (m) = 1/elevation", size=20)
    ax.set_ylabel("$\\tau_{f}$ (Structural index)", size=20)
    # plt.title(r'$\frac{\partial log(f)}{\partial log(z)}$', size=20)
    ax.grid()
    ax.legend()

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
