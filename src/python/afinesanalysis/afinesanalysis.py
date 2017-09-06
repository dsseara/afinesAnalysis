"""
A module to read in and analyze data output from AFiNES

Methods paper that describe AFiNES:

http://www.cell.com/biophysj/abstract/S0006-3495(17)30622-7

A Versatile Framework for Simulating the Dynamic Mechanical Structure of
Cytoskeletal Networks
By:
Simon L. Freedman, Shiladitya Banerjee, Glen M. Hocky, Aaron R. Dinner
Biophysical Journal, 2017

Daniel Seara
Created 09/01/2017
"""
import numpy as np
import os


def readData(directory, *args):
    """
    Reads output .txt files form AFiNES simulations
    INPUT
    -----
    directory = directory where the files are saved, i.e. foo/bar/txt_stack/
    args      = which particles are analyzed. Options are "actins",
                "amotors", "links", or "pmotors". Default is "actins". Can
                be more than one, in which case more than one array is output

    OUTPUT
    ------
    xyt = 3D numpy array that contains positional data in the form
          [time, bead, dim]. The x-positions of all beads in the 4th
          frame is: xyt[4, :, 1]


    TODO
    ----
    Allow export of different information for different particles. i.e. ID of
    links that motors or cross-linkers are attached to.
    """

    fname = os.path.join(directory, args[0])

    with open(fname) as f:
        header = f.readline()
        nparticles = int(header.split()[-1])

    data = np.loadtxt(fname, comments='t', usecols=(0, 1, 3))
    nframes = int(np.ceil(int(data.shape[0] / nparticles)))
    xyt = np.array(np.split(data[:int(nparticles * np.floor(nframes))],
                   np.floor(nframes)))

    return xyt


def plotData(xyid):
    """
    Takes output from readData and plots the output
    """
    return 0
