"""
Module to analyze the connectivity of a cytoskeletal network given a dataset that
contains a list of cross-linkers and what filaments they are attached to.
"""

import numpy as np
import pandas as pd
from scipy.interpolate import Rbf
from scipy.stats import binned_statistic_2d
import os
import fnmatch

def read_afines(filename, nframes, dt):
    '''
    Read outputs form 
    '''
    with open(filename) as f:
        header = f.readline()
        nparticles = int(header.split()[-1])

    xlink = pd.read_csv(filename, comment='t', delimiter='\t',
                        header=None, na_values=-1)
    xlink.columns = ['x0', 'y0', 'x1', 'y1', 'fidx0', 'fidx1', 'lidx0', 'lidx1']
    xlink['t'] = np.arange(0, nframes).repeat(nparticles) * dt

    return xlink
