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
import matplotlib.pyplot as plt
from scipy.interpolate import Rbf
from scipy.stats import binned_statistic_2d
import os
import fnmatch


def readConfigs(directory=None):
    """
    Reads .cfg file and puts outputs as a dictionary

    Parameters
    ----------
    directory: string, optional
        Directory where the .cfg file is saved. If None, then use current
        directory.

    Returns
    -------
    configs: dict
        dictionary with keys as the terms on the left hand side of
        the '=' on each line of the .cfg file in directory, and values
        as the terms on the right hand side. If the rhs is a number,
        it is returned as a float
    """

    if directory is None:
        directory = '.'

    configs = {}
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, '*.cfg'):
            with open(os.path.join(directory, file)) as f:
                for line in f:
                    line = line.split()[0].split('=')
                    try:
                        configs[line[0]] = float(line[1])
                    except:
                        if line[1] == 'true':
                            configs[line[0]] = True
                        elif line[1] == 'false':
                            configs[line[0]] = False
                        else:
                            configs[line[0]] = line[1]

    return configs


def readData(directory, filename='actins.txt'):
    """Reads output .txt files from AFiNES simulations

    Parameters
    ----------
    directory : string
        directory where the files are saved, i.e. foo/bar/txt_stack/
    filename : string (optional)
        which particles are analyzed. Options are "actins.txt",
        "amotors.txt", "links.txt", or "pmotors.txt".
        Default is "actins".

    Returns
    ------
    xyt = array_like
        3D array of position of all particles over time. In shape
        [time, bead, xyID]. x(y)-position of the nth bead in the mth frame
        is: xyt[m, n, 0(1)].
        IF ACTINS: filament id of nth bead at mth frame is: xyt[m, n, 2]

    TODO
    ----
    - Allow input other than actins
    - Allow export of different information for different particles. i.e. ID of
    links that motors or cross-linkers are attached to.
    """

    fname = os.path.join(directory, filename)

    with open(fname) as f:
        header = f.readline()
        nparticles = int(header.split()[-1])

    data = np.loadtxt(fname, comments='t')
    nframes = int(np.ceil(int(data.shape[0] / nparticles)))
    xyt = np.array(np.split(data[:int(nparticles * np.floor(nframes))],
                   np.floor(nframes)))

    return xyt


def interpolateVelocity(txy, dt=10, domainSize=50, nbins=10, minpts=10, dr=1,
                        rbfFunc='gaussian', rbfEps=5, savepath=False):
    """
    Interpolates velocity field of particles to a grid

    Velocity interpolation is done using radial basis functions.

    Parameters
    ----------
    txy : array_like
        particle x and y positions of time, in shape (frame, bead, (x, y))
    dt : scalar, optional
        number of frames between which to measure velocity. Defaults to 10
    domainSize: scalar, optional
        size of domain used, in same units as the positions given. Used to make
        periodic boundary conditions. Defaults to 50
    nbins : scalar, optional
        number of bins in each direction to use in smoothing data.
        Defaults to 10.
    minpts : scalar, optional
        minimum number of beads in a bin to calculate average. Defaults to 10
    dr : scalar, optional
        meshgrid size. Defaults to 1
    rbfFunc : string, optional
        functional form to use for radial basis function interplation. If
        'gaussian' (default), use rbfEps as size of gaussian
    rbfEps : scalar, optional
        standard deviation of gaussian used in rbf if rbfFunc='gaussian'
    savepath : str, optional
        Directory to save outputs. If savepath=None, then use current directory

    Returns
    ------
    xx : array_like
        2D interpolation grid x positions
    yy : array_like
        2D interpolation grid y positions
    uut : array_like
        2+1-D x components of velocity at each grid point
    vvt : array_like
        2+1-D y-components of velocity at each grid point

    See Also
    --------
    scipy.interpolate.Rbf
    """

    # To be used in binning data, avoid areas with low number of beads
    def thresh_mean(arr):
        if len(arr) < minpts:
            return np.nan
        else:
            return arr.mean()

    # Tile positions to simulate periodic boundary conditions
    rxyt = txy + np.array([domainSize, 0])
    lxyt = txy + np.array([-domainSize, 0])
    uxyt = txy + np.array([0, domainSize])
    dxyt = txy + np.array([0, -domainSize])
    urxyt = txy + np.array([domainSize, domainSize])
    ulxyt = txy + np.array([-domainSize, domainSize])
    drxyt = txy + np.array([domainSize, -domainSize])
    dlxyt = txy + np.array([-domainSize, -domainSize])

    txyTiled = np.concatenate([txy, rxyt, lxyt, uxyt, dxyt, urxyt, ulxyt,
                               drxyt, dlxyt], axis=1)
    # Get velocity while enforcing periodicity
    dtxy = txy[dt:] - txy[:-dt]
    dtxyp = (dtxy + domainSize / 2) % domainSize - domainSize / 2

    # Need to repeat 9 times for the 3x3 tiling of xy positions done above
    tuv = np.concatenate([dtxyp / dt] * 9, axis=1)

    edgesx = np.linspace(-domainSize, domainSize, 2 * nbins + 1)
    edgesy = np.linspace(-domainSize, domainSize, 2 * nbins + 1)
    print('imported data')

    # Smooth data by averaging locally and only averaging with enough particles
    mxt = np.stack([binned_statistic_2d(txyTiled[t, :, 0],
                                        txyTiled[t, :, 1],
                                        txyTiled[t, :, 0],
                                        statistic=thresh_mean,
                                        bins=[edgesx, edgesy])[0].ravel()
                    for t in range(len(tuv))])
    myt = np.stack([binned_statistic_2d(txyTiled[t, :, 0],
                                        txyTiled[t, :, 1],
                                        txyTiled[t, :, 1],
                                        statistic=thresh_mean,
                                        bins=[edgesx, edgesy])[0].ravel()
                    for t in range(len(tuv))])
    mut = np.stack([binned_statistic_2d(txyTiled[t, :, 0],
                                        txyTiled[t, :, 1],
                                        tuv[t, :, 0],
                                        statistic=thresh_mean,
                                        bins=[edgesx, edgesy])[0].ravel()
                    for t in range(len(tuv))])
    mvt = np.stack([binned_statistic_2d(txyTiled[t, :, 0],
                                        txyTiled[t, :, 1],
                                        tuv[t, :, 1],
                                        statistic=thresh_mean,
                                        bins=[edgesx, edgesy])[0].ravel()
                    for t in range(len(tuv))])

    print('calculated bin statistic')

    mxyt = np.stack([mxt, myt], axis=2)
    muvt = np.stack([mut, mvt], axis=2)

    # Get interpolated grid positions
    xi = np.arange(-0.5 * domainSize, 0.5 * domainSize + dr, dr)
    yi = np.arange(-0.5 * domainSize, 0.5 * domainSize + dr, dr)
    xx, yy = np.meshgrid(xi, yi)
    uut = np.zeros((len(mxyt), len(xi), len(yi)))
    vvt = np.zeros((len(mxyt), len(xi), len(yi)))

    for t in range(len(mxyt)):
        idxs = np.logical_and(np.isfinite(mxyt[t, :, 0]),
                              np.isfinite(mxyt[t, :, 1]))
        if len(mxyt[t, idxs]) > 0:
            rbfu = Rbf(mxyt[t, idxs, 0], mxyt[t, idxs, 1], muvt[t, idxs, 0],
                       function=rbfFunc, epsilon=rbfEps)
            rbfv = Rbf(mxyt[t, idxs, 0], mxyt[t, idxs, 1], muvt[t, idxs, 1],
                       function=rbfFunc, epsilon=rbfEps)

            uu = rbfu(xx, yy)
            vv = rbfv(xx, yy)

            # xywuwve = np.stack([mxyt[t, idxs, 0].flatten(),
            #                     mxyt[t, idxs, 1].flatten(),
            #                     wu.flatten(),
            #                     wv.flatten(),
            #                     rbfEps * np.ones(wu.shape[0])], axis=1)

            if savepath:
                if savepath is not None:
                    np.savetxt(os.path.join(savepath, 'velInterp_t{0}.txt'.format(t)),
                               [uu.ravel(), vv.ravel()], '%10.5f')
                    np.savetxt(os.path.join(savepath, 'gridInterp_t{0}.txt'.format(t)),
                               [xx.ravel(), yy.ravel()], '%10.5f')
                else:
                    np.savetxt('velInterp_t{0}.txt'.format(t),
                               [uu.ravel(), vv.ravel()])
                    np.savetxt('gridInterp_t{0}.txt'.format(t),
                               [xx.ravel(), vv.ravel()])

            uut[t, ...] = uu
            vvt[t, ...] = vv
            print('t = {0}'.format(t), end='\r')
        else:
            print("No finite data points to use")

    return [xx, yy, uut, vvt]


def massWeightedVelocityDivergence(txy, uut, vvt, dt=10, domainSize=50, dr=1,
                                   savestuff=False):
    """ Calculate mass weighted velocity divergence

    """

    xi = np.arange(-0.5 * domainSize - dr, 0.5 * domainSize + dr, dr)
    yi = np.arange(-0.5 * domainSize - dr, 0.5 * domainSize + dr, dr)
    divV = np.zeros(len(txy))
    for frame, pos in enumerate(uut):
        posx = txy[frame, :, 0]
        posy = txy[frame, :, 1]
        binnedMass, xedges, yedges = np.histogram2d(posx, posy,
                                                    bins=(xi, yi))
        binnedMass = binnedMass.T
        uu = uut[frame, ...]
        vv = vvt[frame, ...]
        divergence = np.gradient(uu)[1] + np.gradient(vv)[0]
        weightedDivV = divergence * binnedMass
        divV[frame] = np.sum(weightedDivV)

    return divV
