[![DOI](https://zenodo.org/badge/266739726.svg)](https://zenodo.org/badge/latestdoi/266739726)

# pydEXP

pyDEXP is a open-source python package aiming at processing potential field data using the dEXP theory formulated by Fedi et al., 2012. The package largely benefits from the imaging methods module of the Fatiando a terra open-source python library.

**References**

* Uieda, L., V. C. Oliveira Jr, and V. C. F. Barbosa (2013), Modeling the Earth with Fatiando a Terra, Proceedings of the 12th Python in Science Conference, pp. 91 - 98.
* Fedi, M., and M. Pilkington (2012), Understanding imaging methods for potential field data, Geophysics, 77(1), G13, doi:10.1190/geo2011-0078.1

### Documentation

Many application examples can be found in the online documentation.
https://dexp-imaging.readthedocs.io/en/latest/

cond
## Install


Create a new environment:
```shell
conda env create -f environment.yml
```

Activate this environment:
```shell
conda activate pydexp
```

Install pydexp from setup.py:
```shell
python setup.py develop
```


Install pygimli from setup.py:
```shell
git clone https://github.com/gimli-org/gimli
cd gimli
conda develop .
```





Update pydexp:
```shell
conda env update --file environment.yml --prune
```


### Dependencies
* pygimli (for the ERT and MALM forward modelling)

Pinned to the version 1.1.0: conda install -c gimli -c conda-forge pygimli=1.1.0

* fatiando 0.5 (for derivation of the potential field and upward continuation)

Pinned to the version 0.5 (git clone the repo for dev purpose)

Others dependencies are detailed in the requirements.txt file.

### Cite

One article stemming from this work is in preparation for [JGR](https://agupubs.onlinelibrary.wiley.com/journal/21699356?utm_source=google&utm_medium=paidsearch&utm_campaign=R3MR425&utm_content=EarthSpaceEnvirSci&gclid=EAIaIQobChMI7rHW38Oj-QIVxITVCh1ygABKEAAYASAAEgLZYPD_BwE). All data and codes to reproduce results are located in the folder notebook_JGR
