Getting Started
===============

.. The getting-started should show some primary use cases in more detail. The reader will follow a step-by-step procedure to set-up a working prototype

pyDEXP aims to process Mise-à-la-masse (MALM) datasets for a variety of applications. pyDEXP has been initially developed for plant root imaging.

The simpliest processing can be achieved with the python API.
You'll first need to import the pyDEXP package::

	import lib.dEXP as dEXP
	import lib.plot_dEXP as pEXP

Then after loading you data make basics operation such as upward continuation (Uieda et al., 2013, 2018)::

	mesh, label_prop = dEXP.upwc(xp, yp, zp, U, shape,
		         zmin=0, zmax=max_elevation, nlayers=nlay,
		         qorder=qorder)
                 
Searching for ridges (Fedi et al., 2012)::

	dfI,dfII, dfIII = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,
		                              label=label_prop,
		                              fix_peak_nb=2,


and finally plot the results::

	fig = plt.figure()
	ax = plt.gca()
	pEXP.plot_xy(mesh, label=label_prop, ax=ax)
	pEXP.plot_ridges_harmonic(dfI,dfII,dfIII,ax=ax)
		                              
		                                      
More examples are available in the Example gallery section.



**References**

Uieda, L., V. C. Oliveira Jr, and V. C. F. Barbosa (2013), Modeling the Earth with Fatiando a Terra, Proceedings of the 12th Python in Science Conference, pp. 91 - 98.

Uieda, L. (2018). Verde: Processing and gridding spatial data using Green's functions. Journal of Open Source Software, 3(29), 957. doi:10.21105/joss.00957

Fedi, M., and M. Pilkington (2012), Understanding imaging methods for potential
field data, Geophysics, 77(1), G13, doi:10.1190/geo2011-0078.1



