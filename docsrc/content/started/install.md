Which Python?
-------------

For the moment, pydEXP is only tested on Python 3.7. First, make sure you have all dependencies (see below) installed. 



Clone the gitlab repository::

    git clone https://github.com/BenjMy/dEXP_imaging.git

From the same `dEXP_imaging` directory you can import the module from source using python. 


Dependencies
------------

pydEXP requires the following dependencies for running:

* numpy
* scipy
* matplotlib
* pandas
* mpl_axes_aligner

In order to make upward continuation and derivate the field, pyDEXP uses:
* `fatiando <https://legacy.fatiando.org/>`_


Optionnal Third party packages
------------------------------

In order to forward model the geolectrical data, pyDEXP uses:

* Pygimli
* Resipy

.. Testing the install
.. -------------------

You can test the installation running one of the exemple.
