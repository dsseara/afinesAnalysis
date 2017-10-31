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


def globalConnectivity(xlinker, nfilaments):
    '''
    Measure bulk connectivity of cytoskeletal network

    Parameters
    ----------
    xlinker : DataFrame
        pandas DataFrame with the xlinker information from reading
        either amotors.txt or pmotors.txt
    nfilaments : array_like
        number of filaments in the simulation. Same size as nconnected

    Results
    -------
    connectivity : array_like
        average connectivity measured for each element of nconnected and
        nfilaments as: 2*(nconnected)/nfilaments

    See also
    --------

    Example
    -------

    '''
    nconnected = (xlinker.dropna(subset=['fidx0', 'fidx1'])
                  .groupby('t').count().values[:, 0])
    connectivity = 2 * nconnected / nfilaments
    return connectivity


def clusterFilaments(xlinker, filaments):
    '''
    Get ids of filaments that are connected together

    Treat each doubly-bound xlinker as an edge on a graph where
    each filament is a node on that graph. Then find all the
    connected components of the graph, and label each filament (node)
    by which connected component it belongs to.

    Parameters
    ----------
    xlinker : DataFrame
        pandas DataFrame with the xlinker information from reading either
        amotors.txt or pmotors.txt
    filaments : array_like
        array of number of filaments in each frame

    Results
    -------
    xlinker : DataFrame
        same Dataframe as input, but with extra column that labels

    See also
    --------

    Example
    -------

    '''
    filaments['clusterID'] = (-1 * np.ones(len(filaments))).astype(int)

    for time in np.unique(filaments.t):
        g = nx.Graph()
        edges = (xlinker[xlinker.t == time]
                 .dropna(subset=['fidx0', 'fidx1'])[['fidx0', 'fidx1']].values)
        g.add_edges_from(edges)
        ccs = nx.connected_components(g)
        for index, cc in enumerate(ccs):
            filaments.loc[(filaments.t == time).values &
                          (filaments.filamentID.isin(cc)).values,
                          'clusterID'] = index

    # for (time, xldata), (time, filadata) in zip(xlinker.groupby('t'),
    #                                             filaments.groupby('t')):
    #     # Fill in new column with -1 if not connected to any other filament
    #     g = nx.Graph()

    #     edges = (xldata.dropna(subset=['fidx0', 'fidx1'])[['fidx0', 'fidx1']]
    #              .values)
    #     g.add_edges_from(edges)
    #     # Get connected components
    #     ccs = nx.connected_components(g)
    #     for index, cc in enumerate(ccs):
    #         filadata.loc[filadata.filamentID.isin(cc), 'clusterID'] = index
