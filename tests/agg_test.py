import numpy as np

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"

def aggfunc(x):
    print(x)
    return np.array([np.nanmin(x), np.nansum(x)])