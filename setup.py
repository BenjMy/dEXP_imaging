#!/usr/bin/env python
# import os
# import glob

from setuptools import setup, find_packages

version_long = '0.1.0.dev0'

if __name__ == '__main__':
    setup(
        name='dep_imaging',
        version=version_long,
        description='Potential field inversion of MALM data',
        long_description=open('Readme.md', 'r').read(),
        long_description_content_type="text/markdown",
        author='Benjamin Mary',
        author_email='benjamin.mary@unipd.it',
        license='MIT',
        #packages=find_packages("lib"),
        # package_dir={'': 'lib'},
        package_dir={   package_1:"lib", 
                        package_2:"examples",
                        package_3:"fatiando"}

        # install_requires=[
        #     'dicttoxml',
        #     'jupyter',
        #     'ipywidgets',
        #     'PyQT5',
        # ],
        # classifiers=(
        #     "Programming Language :: Python :: 3",
        #     "License :: OSI Approved :: MIT License",
        #     "Operating System :: OS Independent",
        # ),
    )
