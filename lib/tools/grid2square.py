# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 09:57:38 2020

@author: Benjamin
"""


import pyvista as pv
import os


def squaremat(r=20):
    # MainPath= 'E:/Padova/Experiments/GC_2019_Landfill_Porto_MALM'
    # os.chdir(MainPath)
    # grid = pv.read('Ano_field_EA.vtk')
    # grid.set_active_scalars("_Marker")
    # cpos = grid.plot()

    # bodies = grid.split_bodies()
    # for i, body in enumerate(bodies):
    #     print('Body {} volume: {:.3f}'.format(i, body.volume))

    # position of the anomaly
    ya = 5.036227e6
    xa = (284272.0 + 284250.00) / 2

    # landfill geometry
    # utm coordinates
    coords_liner = [
        (284046.43, 5036328.39),
        (284132.24, 5036277.32),
        (284146, 5036297),
        (284281.07, 5036214.66),
        (284317.46, 5036281.81),
        (284245, 5036313),
        (284097.55, 5036411.08),
    ]
    coords_liner = np.asarray(coords_liner)

    slope = (coords_liner[2, 1] - coords_liner[3, 1]) / (
        coords_liner[2, 0] - coords_liner[3, 0]
    )  # (yb-ya)/(xb-xa)
    x1 = xa + r * np.cos(slope)
    y1 = ya + r * np.sin(slope)

    x2 = xa - r * np.cos(slope)
    y2 = ya - r * np.sin(slope)

    plt.figure(figsize=(20, 10))
    ax = plt.subplot()
    plt.plot(coords_liner[:, 0], coords_liner[:, 1], "*-")
    for i in range(len(coords_liner[:, 0])):
        ax.annotate(str(i), (coords_liner[i, 0], coords_liner[i, 1]))

    plt.scatter(xa, ya, c="red")
    plt.scatter(x1, y1, c="blue")
    plt.scatter(x2, y2, c="blue")
    plt.axis("square")


return x1, y1, x2, y2
