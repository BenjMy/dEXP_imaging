"""
Imaging methods for potential fields.

Implements the DEXP method described in Fedi and Pilkington (2012).

.. note::

	This is part of a larger project aiming at inverting current sources density (see more at: https://icsd-dev.readthedocs.io/en/latest/)


**References**

Fedi, M., and M. Pilkington (2012), Understanding imaging methods for potential
field data, Geophysics, 77(1), G13, doi:10.1190/geo2011-0078.1

----
"""

import math
import numpy as np
import matplotlib.pyplot as plt


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

	u_cor = u - u_B # correct total measured potential from influence of B

	plt.figure()
	plt.subplot(1,3,1)
	plt.tricontourf(x, y, dist, 50, cmap='RdBu_r')
	cbar = plt.colorbar(orientation='vertical', aspect=50)
	cbar.set_label('Distance from B (m)')
	plt.tight_layout()
	plt.axis('square')

	plt.subplot(1,3,2)
	plt.tricontourf(x, y, u_B, 50, cmap='RdBu_r')
	cbar = plt.colorbar(orientation='vertical', aspect=50)
	cbar.set_label('$u_{B}$')
	plt.tight_layout()
	plt.axis('square')

	plt.subplot(1,3,3)
	plt.tricontourf(x, y, u_cor, 50, cmap='RdBu_r')
	cbar = plt.colorbar(orientation='vertical', aspect=50)
	cbar.set_label('$u = U_{T} - u_{B}$ (V)')
	plt.tight_layout()
	plt.axis('square')

	return u_cor

def ridges():
	"""
	Text here

	Parameters:

	* a
		Text here

	Returns:

	* BB : 
		Text here - Panda dataframe containing RI, RII and RII

	"""
	return

def geom_Z0():
	"""
	Text here

	Parameters:

	* a
		Text here

	Returns:

	* BB : 
		Text here

	"""
	return

def scalFUN():
	"""
	Text here

	Parameters:

	* a
		Text here

	Returns:

	* BB : 
		Text here

	"""
	return Tau, q, SI

def dEXP(x, y, z, data, shape, zmin, zmax, nlayers, qorder=0, SI=1):
	"""
	DEXP model (Fedi, 2012).

	Calculates a physical property distribution given potential field data on a
	**regular grid**. Uses depth weights.
	Parameters:

	* x, y : 1D-arrays
		The x and y coordinates of the grid points
	* z : float or 1D-array
		The z coordinate of the grid points
	* data : 1D-array
		The potential field at the grid points
	* shape : tuple = (ny, nx)
		The shape of the grid
	* zmin, zmax : float
		The top and bottom, respectively, of the region where the physical
		property distribution is calculated
	* nlayers : int
		The number of layers used to divide the region where the physical
		property distribution is calculated
	* qorder : float
		The order of the vertical derivative
	* SI : float
		The structural index

	Returns:

	* mesh : :class:`fatiando.mesher.PrismMesh`
		The estimated physical property distribution set in a prism mesh (for
		easy 3D plotting)

	"""
	mesh = _makemesh(x, y, shape, zmin, zmax, nlayers)
	# This way, if z is not an array, it is now
	z = z * numpy.ones_like(x)
	freq, dataft = _getdataft(x, y, data, shape)
	# Remove the last z because I only want depths to the top of the layers
	depths = mesh.get_zs()[:-1]
	weights = (np.abs(depths)) ** ((SI+qorder)/2)
	density = []
	# Offset by the data z because in the paper the data is at z=0
	for depth, weight in zip(depths - z[0], weights):

		# continued field calculation
		upw_f= transform.upcontinue(xp, yp, gz, shape, depth)

		# qorder vertical derivate of the continued field
		upw_f_dq = transform.derivz(xp, yp, upw_f, shape,order=qorder)
		
		# the continued field weigted (=DEXP)
		upw_f_dq_w= upw_f_dq*weight

	mesh.addprop('density', numpy.array(upw_f_dq_w))
	return mesh

# def auto_dEXP():


