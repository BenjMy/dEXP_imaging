"""
Magnetic field data analysis using pyDEXP: effect of derivative order q
-----------------------------------------------------------------------

This code shows a step-by-step processing of potential field imaging aiming at giving an estimate of magnetic sources positions and depth using the dEXP tranformation method.
dEXP method implementation from Fedi et al. 2012. 
Calculations used :mod:`dEXP`, while plotting use the :mod:`plot_dEXP` module.

The model data was created using geometric objects from :mod:`fatiando.mesher`. The forward simulation of the data was done using :mod:`fatiando.gravmag` module.

Sources locations:
    - S_{A} = [10e3,10e3,2e3] # xyz coordinates
    - S_{B} = [25e3,10e3,1e3]

Sources properties: 
    - radius = 1.5e3
    - inc = 50
    - dec = -30

.. note::

    This is part of a larger project aiming at inverting current sources density (see more at: https://icsd-dev.readthedocs.io/en/latest/)

**References**

Uieda, L., V. C. Oliveira Jr, and V. C. F. Barbosa (2013), Modeling the Earth with Fatiando a Terra, Proceedings of the 12th Python in Science Conference, pp. 91 - 98.

Uieda, L. (2018). Verde: Processing and gridding spatial data using Green's functions. Journal of Open Source Software, 3(29), 957. doi:10.21105/joss.00957

Fedi, M., and M. Pilkington (2012), Understanding imaging methods for potential
field data, Geophysics, 77(1), G13, doi:10.1190/geo2011-0078.1

"""
import matplotlib.pyplot as plt
import numpy as np
import dEXP as dEXP
from dEXP import _fit
import plot_dEXP as pEXP
import set_parameters as para
import examples.magnetic.fwdmag.fwd_mag_sphere as magfwd


#%%
# Write code showing the effect of the derivative numbers on the nb of ridges and on the solution accuracy