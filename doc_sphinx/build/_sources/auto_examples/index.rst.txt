:orphan:



.. _sphx_glr_auto_examples:

Introduction
============

Below is a gallery of examples for different physics: **gravimetry**, **magnetic** and **geoelectrical methods (SP and active ERT)**. For each we applied the DEXP transformation with a step-by-step processing description in order to estimate of the sources positions. The DEXP method is the implementation of Fedi et al. 2012. Calculations used :mod:`dEXP`, while plotting use the :mod:`plot_dEXP` module.



    
    
**References**

- Uieda, L., V. C. Oliveira Jr, and V. C. F. Barbosa (2013), Modeling the Earth with Fatiando a Terra, Proceedings of the 12th Python in Science Conference, pp. 91 - 98.

- Uieda, L. (2018). Verde: Processing and gridding spatial data using Green's functions. Journal of Open Source Software, 3(29), 957. doi:10.21105/joss.00957

- Fedi, M., and M. Pilkington (2012), Understanding imaging methods for potential
field data, Geophysics, 77(1), G13, doi:10.1190/geo2011-0078.1



.. raw:: html

    <div class="sphx-glr-clear"></div>



.. _sphx_glr_auto_examples_gravimetry:

Gravimetric potential field data
================================


.. note::

	The gravimetric model data was created using geometric objects from fatiando.mesher. The 		forward simulation of the data was done using fatiando.gravmag module.

	Sources properties: 
	    
	    
* estimate of the **depth of the density** anomaly using the dEXP tranformation method 



.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This code shows a step-by-step processing of potential field imaging aiming at giving an estima...">

.. only:: html

 .. figure:: /auto_examples/gravimetry/images/thumb/sphx_glr_run_grav_thumb.png
     :alt: Example of gravimetric potential field data analysis using pyDEXP

     :ref:`sphx_glr_auto_examples_gravimetry_run_grav.py`

.. raw:: html

    </div>


.. toctree::
   :hidden:

   /auto_examples/gravimetry/run_grav
.. raw:: html

    <div class="sphx-glr-clear"></div>



.. _sphx_glr_auto_examples_magnetic:

Magnetic potential field data
=============================

.. note::

	The magnetic model data was created using geometric objects from fatiando.mesher. The 		forward simulation of the data was done using fatiando.gravmag module.

	Sources properties: 
	    * radius = 1.5e3
	    * inc = 50
	    * dec = -30


* Identify 2 depths of sources produces by 2 distincts magnetic sources using the geometrical method


    



.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This code shows a step-by-step processing of potential field imaging aiming at giving an estima...">

.. only:: html

 .. figure:: /auto_examples/magnetic/images/thumb/sphx_glr_run_mag_2sources_DEXP_thumb.png
     :alt: Magnetic field data analysis using pyDEXP: a 2-sources case

     :ref:`sphx_glr_auto_examples_magnetic_run_mag_2sources_DEXP.py`

.. raw:: html

    </div>


.. toctree::
   :hidden:

   /auto_examples/magnetic/run_mag_2sources_DEXP
.. raw:: html

    <div class="sphx-glr-clear"></div>



.. _sphx_glr_auto_examples_malm:

Mise-à-la-masse analysis using the dEXP theory
==============================================

The theory of the dEXP theory is applied here for the first time (to our knowledges) for active geoelectrical methods. The theory behind this choice and required adjustements are explained in the documentation section.

[preprocessing](./examples/malm/clearB.py)

[preprocessing](clearB.py)

[preprocessing](examples/malm/clearB)

We show here: 

* The [preprocessing](examples/malm/clearB.py) of MALM for dEXP;
* a sensitivity analysis to **anomaly depth** (3d) of the dEXP theory apply to Mise-à-la-Masse data;
* effect of **noise level**;
* **in prep**: a case study where we used a Mise-à-la-Masse survey to identify a leakage in a landfill;
* **in prep**: effect of borehole data and downward continuation.

    



.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This code shows a step-by-step processing of potential field imaging aiming at giving an estima...">

.. only:: html

 .. figure:: /auto_examples/malm/images/thumb/sphx_glr_run_sens_depth_MALM_dexpratio_thumb.png
     :alt: Sensitivity analysis of DEXP to anomaly depth on Mise-a-la-masse

     :ref:`sphx_glr_auto_examples_malm_run_sens_depth_MALM_dexpratio.py`

.. raw:: html

    </div>


.. toctree::
   :hidden:

   /auto_examples/malm/run_sens_depth_MALM_dexpratio

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This code shows a step-by-step processing of potential field imaging aiming at giving an estima...">

.. only:: html

 .. figure:: /auto_examples/malm/images/thumb/sphx_glr_run_sens_depth_MALM_dexpratio_q12_thumb.png
     :alt: Sensitivity analysis of DEXP to anomaly depth on Mise-a-la-masse

     :ref:`sphx_glr_auto_examples_malm_run_sens_depth_MALM_dexpratio_q12.py`

.. raw:: html

    </div>


.. toctree::
   :hidden:

   /auto_examples/malm/run_sens_depth_MALM_dexpratio_q12

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This code shows a step-by-step processing of potential field imaging aiming at giving an estima...">

.. only:: html

 .. figure:: /auto_examples/malm/images/thumb/sphx_glr_run_sens_sigma_MALM_dexpratio_thumb.png
     :alt: Sensitivity analysis of DEXP to depth on Mise-a-la-masse

     :ref:`sphx_glr_auto_examples_malm_run_sens_sigma_MALM_dexpratio.py`

.. raw:: html

    </div>


.. toctree::
   :hidden:

   /auto_examples/malm/run_sens_sigma_MALM_dexpratio

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This code shows a step-by-step processing of potential field imaging aiming at giving an estima...">

.. only:: html

 .. figure:: /auto_examples/malm/images/thumb/sphx_glr_run_sens_MALM_thumb.png
     :alt: Sensitivity analysis of DEXP to contrast of resistivity on Mise-a-la-masse

     :ref:`sphx_glr_auto_examples_malm_run_sens_MALM.py`

.. raw:: html

    </div>


.. toctree::
   :hidden:

   /auto_examples/malm/run_sens_MALM

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This code shows a step-by-step processing of potential field imaging aiming at giving an estima...">

.. only:: html

 .. figure:: /auto_examples/malm/images/thumb/sphx_glr_run_sens_noise_MALM_thumb.png
     :alt: Sensitivity analysis of DEXP to noise level

     :ref:`sphx_glr_auto_examples_malm_run_sens_noise_MALM.py`

.. raw:: html

    </div>


.. toctree::
   :hidden:

   /auto_examples/malm/run_sens_noise_MALM

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="    Application on a anomaly of electrical resistivity (1000 Omega.m^{-1} of contrast).     The...">

.. only:: html

 .. figure:: /auto_examples/malm/images/thumb/sphx_glr_clearB_thumb.png
     :alt: Preprocessing of MALM for dEXP

     :ref:`sphx_glr_auto_examples_malm_clearB.py`

.. raw:: html

    </div>


.. toctree::
   :hidden:

   /auto_examples/malm/clearB
.. raw:: html

    <div class="sphx-glr-clear"></div>



.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
