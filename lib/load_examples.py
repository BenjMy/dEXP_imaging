"""
dEXP for potential field source depth estimation
"""
from urllib.request import urlretrieve
import os
import tempfile

DEXP_ExampleDataPath='BenjMy/dEXP_imaging/master/lib/examples/'
# Example data repository
exampleDataRepository = ''.join((
    'https://raw.githubusercontent.com/',  # RAW files
    DEXP_ExampleDataPath,  # Organization and repository
    #'master/'  # Branch
))
import pygimli as pg

def getExampleFile(path, load=False, verbose=False):
    """Download and return a filename to the example repository.
    Parameters
    ----------
    path: str
        Path to the remote repo
    Returns
    -------
    filename: str
        Filename to the data content
    """
    url = exampleDataRepository + path
    fileName = os.path.join(tempfile.gettempdir(), DEXP_ExampleDataPath, path)

    if os.path.exists(fileName):
        if verbose:
            os.makedirs(os.path.dirname(fileName), exist_ok=True)
        tmp = urlretrieve(url, fileName)
            
    return tmp[0]
