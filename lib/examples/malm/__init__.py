import numpy as np

import pygimli as pg
import pygimli.meshtools as mt
from pygimli.viewer.pv import drawSensors

elecs = np.loadtxt('elecs_ref.txt')
sensors = np.zeros((np.shape(elecs)[0],3))
sensors[:,0:2] = elecs


plc0 = mt.createCube(size=[400, 400, 9], pos=[min(sensors[:,0]), min(sensors[:,1]), -9/2], marker=1, boundaryMarker=1)
plc1 = mt.createCube(size=[400, 400, 7], pos=[min(sensors[:,0]), min(sensors[:,1]),-9-7/2], marker=2, boundaryMarker=1)
plc2 = mt.createCube(size=[400, 400, 4], pos=[min(sensors[:,0]), min(sensors[:,1]),-9-7-2], marker=1, boundaryMarker=1)

plc= plc0 + plc1 + plc2
mesh = mt.createMesh(plc)

mesh.exportVTK('mesh.vtk')
mesh.exportVTK('elecs.vtk', sensors)


plotter, _ = pg.show(mesh, mesh.cellMarkers(), alpha=0.6, hold=True, notebook=True)
drawSensors(plotter, sensors, diam=5, color='yellow')
plotter.show()


bounds = [2,4.5, 2,4.5, 1,3]
clipped = dataset.clip_box(bounds)

p = pv.Plotter()
p.add_mesh(dataset, style='wireframe', color='blue', label='Input')
p.add_mesh(clipped, label='Clipped')
p.add_legend()
p.show()


n_sensors = 8
sensors = np.zeros((n_sensors, 3))
sensors[0, 0] = 15
sensors[0, 1] = -10
sensors[1:, 0] = -15
sensors[1:, 1] = np.linspace(-15, 15, n_sensors - 1)


