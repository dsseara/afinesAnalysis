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
import pandas as pd
from scipy.interpolate import Rbf
from scipy.stats import binned_statistic_2d
import os
import glob
import warnings


def readConfigs(filename=None):
    """
    Reads .cfg file and puts outputs as a dictionary

    Parameters
    ----------
    filename: string, optional
        .cfg file to read. If None, then use .cfg file in current
        directory

    Returns
    -------
    configs: dict
        keys and values are parameter names and parameter inputs in filename.
        If parameter input is a number, it is returned as a float
    """

    if filename is None:
        if len(glob.glob('./*cfg')) > 1:
            raise ValueError('More than 1 cfg file in current directory.')
        elif len(glob.glob('./*cfg')) == 0:
            raise ValueError('0 cfg files in current directory')
        else:
            filename = glob.glob('./*cfg')[0]

    if not os.path.isfile(filename):
        raise ValueError('Could not find ' + filename)

    configs = {}
    with open(filename) as f:
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


def readData(filename, configs, dataframe=True):
    """Reads output .txt files from AFiNES simulations

    Parameters
    ----------
    filename : string
        name of file to load. Optional are "actins.txt",
        "amotors.txt", "links.txt", or "pmotors.txt".
        Defaults to "actins".
    dataframe : bool, optional
        If True, read the data as a pandas DataFrame. Else, read
        in as numpy array. Defaults to True.

    Returns
    ------
    data : DataFrame or array_like
        parsed from filename

    TODO
    ----
    """

    txtFile = filename.split(os.path.sep)[-1]
    if txtFile.split('.')[0] not in ['actins', 'amotors', 'pmotors', 'links']:
        raise ValueError(txtFile + ' not recognized. Options are actins.txt, \
                         amotors.txt, links.txt, or pmotors.txt')

    # nframes = configs['nframes']
    dt_frame = (configs['tfinal'] - configs['tinit']) / nframes

    with open(filename) as f:
        header = f.readline()
        nparticles = int(header.split()[-1])

    if nparticles > 0:
        if dataframe:
            data = pd.read_csv(filename, header=None, comment='t', delimiter='\t')
            nframes = len(data.iloc[:, 0]) / nparticles
            if txtFile in 'actins.txt':
                data.columns = ['x', 'y', 'link_length', 'fid']
            elif txtFile in ['amotors.txt', 'pmotors.txt']:
                data.columns = ['x0', 'y0', 'x1', 'y1',
                                'fidx0', 'fidx1', 'lidx0', 'lidx1']
            data['t'] = np.arange(0, nframes).repeat(nparticles) * dt_frame
        else:
            data = np.loadtxt(filename, comments='t')
            data = np.array(np.split(data[:int(nparticles * np.floor(nframes))],
                            np.floor(nframes)))
    else:
        warnings.warn(filename + ' is empty. Creating empty DataFrame in its place')
        data = pd.DataFrame()

    return data


def interpolateVelocity(data, configs, dt=10, nbins=10, minpts=10, dr=1,
                        rbfFunc='gaussian', rbfEps=5, savepath=False):
    """
    Interpolates velocity field of particles to a grid

    Velocity interpolation is done using radial basis functions.

    Parameters
    ----------
    data : DataFrame or array_like
        output from readData. If numpy array, shape (nframes, nbeads, 2).
    configs : dict
        dictionary of configurations used to run AFINES. See readConfigs()
    dt : scalar, optional
        number of frames between which to measure velocity. Defaults to 10
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
    afinesanalysis.afinesanalys.readConfigs()
    """

    def _thresh_mean(arr):
        '''
        Helper function for interpolateVelocity to calculate mean of array if
        threshold number of points available

        Parameters
        ----------
        arr : array_like
            1D array-like object
        Returns
        -------
        thresh_mean : scalar
            mean of arr iff len(arr) >= minimum number of points
        '''
        if len(arr) < minpts:
            return np.nan
        else:
            return arr.mean()

    # Load data depending on what type it is
    if type(data) == np.ndarray:
        txy = data[..., :2]
    elif type(data) == pd.core.frame.DataFrame:
        txy = data.iloc[:, :2].values
        txy = txy.reshape((int(configs['nframes']),
                           int(configs['npolymer'] * configs['nmonomer']),
                           2))

    # Only considers square domains for now
    domainSize = configs['xrange']

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

    # Smooth data by averaging locally only in regions with enough particles
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
                    np.savetxt(os.path.join(savepath,
                               'velInterp_t{0}.txt'.format(t)),
                               [uu.ravel(), vv.ravel()], '%10.5f')
                    np.savetxt(os.path.join(savepath,
                               'gridInterp_t{0}.txt'.format(t)),
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
