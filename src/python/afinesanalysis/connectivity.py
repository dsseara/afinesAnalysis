"""
Module to analyze the connectivity of a cytoskeletal network given a dataset that
contains a list of cross-linkers and what filaments they are attached to.
"""

import numpy as np
import pandas as pd
from scipy.interpolate import Rbf
from scipy.stats import binned_statistic_2d
import os
import networkx as nx

def read_afines(filename, nframes, dt):
    '''
    Reads in motor.txt file from afines output

    Parameters
    ----------
    filename : str
        name of file to load. Filename should be amotors.txt or pmotors.txt
    nframes : int
        number of time points for the data
    dt : float
        time between time points

    Results
    -------
    motors : DataFrame
        pandas dataframe that contains the motor positions and head states
    '''
    with open(filename) as f:
        header = f.readline()
        nparticles = int(header.split()[-1])

    motors = pd.read_csv(filename, comment='t', delimiter='\t',
                        header=None, na_values=-1)
    motors.columns = ['x0', 'y0', 'x1', 'y1', 'fidx0', 'fidx1', 'lidx0', 'lidx1']
    motors['t'] = np.arange(0, nframes).repeat(nparticles) * dt

    return motors

def globalConnectivity()

def clusterFilaments():
    '''
    Get ids of filaments that are connected together

    Parameters
    ----------
     param1 : type
        description

    Results
    -------
     output1 : type
        description

    See also
    --------

    Example
    -------

    '''
