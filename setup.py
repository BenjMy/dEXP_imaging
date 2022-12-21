import setuptools
from setuptools import setup, find_packages

with open("Readme.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="pyCATHY",
    version="0.0.1",
    author="B. Mary",
    author_email="bmary@lbl.gov",
    description="Potential field inversion of MALM data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/BenjMy/dEXP_imaging",
    packages=setuptools.find_packages("lib","fatiando"),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License " + "v3 (LGPLv3)",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
    install_requires=["numpy","future","matplotlib","scipy","pandas",
    "mpl_axes_aligner","csaps"]
) 


<<<<<<< Updated upstream
if __name__ == '__main__':
    setup(
        name='dexp_imaging',
        version=version_long,
        description='Potential field inversion of MALM data',
        long_description=open('Readme.md', 'r').read(),
        long_description_content_type="text/markdown",
        author='Benjamin Mary',
        author_email='benjamin.mary@unipd.it',
        license='MIT',
        packages=find_packages("lib","dEXP"),
        #package_dir={'': 'lib'},
        install_requires=['csaps','future']
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
=======
>>>>>>> Stashed changes
