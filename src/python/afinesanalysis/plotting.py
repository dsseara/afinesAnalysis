"""
Module to plot outputs of afines
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import pandas as pd
import os


def animateQuiver(xx, yy, uut, vvt, taxis=0):
    '''
    creates an animated quiver plot of a vector field (u(x,y), v(x,y))

    Parameters
    ----------
    xx : array_like
        2D array of grid points
    yy : array_like
        2D array of grid points
    uut : array_like
        x-component of vector field over time
    vvt : array_like
        y-component of vector field over time
    taxis : scalar, optional
        axis along which time is measured. Defaults to 0

    Returns
    -------
    anim : animation object
        animation object must be returned to display animation when
        calling animation.FuncAnimation from a function.

    '''

    # Make sure time is first dimension
    uut = np.rollaxis(uut, taxis)
    vvt = np.rollaxis(vvt, taxis)

    fig, ax = plt.subplots()
    Q = ax.quiver(xx, yy, uut[0, ...], vvt[0, ...], pivot='mid')
    anim = animation.FuncAnimation(fig, _update_quiver, frames=uut.shape[0],
                                   fargs=(Q, uut, vvt), blit=False)

    return anim


def _update_quiver(frame, Q, uut, vvt):
    '''
    Helper function for animateQuiver
    '''
    Q.set_UVC(uut[frame, ...], vvt[frame, ...])


def filaments(filas, configs, dt=1, dfilament=1, savepath=False):
    """Takes output from readData and plots the output as a series of .pngs
    Only works for actin data

    Parameters
    ----------
    data : pandas DataFrame
        Output of readData containing filament data.
    configs : dict
        dictionary of afines configuration from readConfigs
    dt : scalar
        plot every dt-th frame
    dfilament : scalar
        plot every dfilament-th filament
    savepath : bool or string, optional
        path to save series of .pngs. Creates subfolder savepath in current
        dir. Creates subfolder imgSeq in current dir if True. If False, does
        not save. Defaults to False

    Returns
    -------
    saves series of .pngs of simulation data

    See also
    --------
    afinesanalysis.readConfigs()
    """

    if savepath is True:
        savepath = os.path.join(os.curdir, 'imgSeq')

    if savepath:
        if not os.path.exists(os.path.join(os.curdir, savepath)):
            os.mkdir(os.path.join(os.curdir, savepath))
        savepath = os.path.join(os.curdir, savepath)

    rangex = configs['xrange']
    rangey = configs['yrange']

    for tind, time in enumerate(pd.unique(filas.t)[::dt]):
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        data = filas.loc[filas.t == time]
        for idindex, ID in enumerate(pd.unique(data.fid)[::dfilament]):
            actin = data[data.fid == ID][['x', 'y']].values
            dists = np.linalg.norm(np.diff(actin, 1, 0), axis=1)
            bools = np.reshape(dists > np.mean((rangex, rangey)) * 0.9,
                               [dists.size, 1])
            mask = np.vstack([np.hstack([bools, bools]), [False, False]])
            masked_actin = np.ma.MaskedArray(actin, mask)
            ax.plot(masked_actin[:, 0], masked_actin[:, 1], 'k', linewidth=0.5)
        # End loop over all filaments
        ax.set_aspect('equal')
        ax.set_xlim([-rangex / 2, rangex / 2])
        ax.set_ylim([-rangey / 2, rangey / 2])
        plt.tight_layout()

        if savepath:
            fig.savefig(os.path.join(savepath, 't{n:0.1f}.png'.format(n=tind)))
            plt.close(fig)
    # End loop over all times
# End filaments


def all(filamentData, pmotorData, amotorData, configs,
        dt=1, dfilament=1, dpmotors=1, damotors=1, savepath=False):
    """Takes output from readData and plots the output as a series of .pngs
    Plots filaments, crosslinkers (pmotors) and motors (amotors)

    Parameters
    ----------
    filamentData : pandas DataFrame
        Output of readData containing filament data.
    pmotorData : pandas DataFrame
        Output of readData containing pmotor data.
    amotorData : pandas DataFrame
        Output of readData containing amotor data.
    configs : dict
        dictionary of afines configuration from readConfigs
    dt : scalar
        plot every dt-th frame. Defaults to 1
    dfilament : scalar
        plot every dfilament-th filament. Defaults to 1
    dpmotor : scalar
        plot every dpmotor-th pmotor. Defaults to 1
    damotor : scalar
        plot every damotor-th amotor. Defaults to 1
    savepath : {bool, string}, optional
        path to save series of .pngs. Creates subfolder savePath/imgSeq. If,
        False, does not save. If True, uses current directory.
        Defaults to False.

    Returns
    -------
    saves series of .pngs of simulation data

    See also
    --------
    readConfigs()
    """
    if savepath is True:
        savepath = os.path.join(os.curdir, 'imgSeq')

    if savepath:
        if not os.path.exists(os.path.join(os.curdir, savepath)):
            os.mkdir(os.path.join(os.curdir, savepath))
        savepath = os.path.join(os.curdir, savepath)

    rangex = configs['xrange']
    rangey = configs['yrange']

    for tind, time in enumerate(pd.unique(filamentData.t)[::dt]):
        fig, ax = plt.subplots()

        data = filamentData.loc[filamentData.t == time]
        for idindex, ID in enumerate(pd.unique(data.fid)[::dfilament]):
            actin = data[data.fid == ID][['x', 'y']].values
            dists = np.linalg.norm(np.diff(actin, 1, 0), axis=1)
            bools = np.reshape(dists > np.mean((rangex, rangey)) * 0.9,
                               [dists.size, 1])
            mask = np.vstack([np.hstack([bools, bools]), [False, False]])
            masked_actin = np.ma.MaskedArray(actin, mask)
            ax.plot(masked_actin[:, 0], masked_actin[:, 1], 'k', linewidth=0.5)
        # End loop over all filaments

        # plot motors and cross-linkers. Note that the x/y position of the
        # second head is really x0/y0 + x1/y1, not just x1/y1
        if not amotorData.empty:
            amotorxs = amotorData.loc[amotorData.t == time][['x0', 'x1']].values
            amotorxs[:, 1] = amotorxs.sum(axis=1)
            # amotorxs = amotorxs[amotorxs[:, 1] < (rangex / 2), :]

            amotorys = amotorData.loc[amotorData.t == time][['y0', 'y1']].values
            amotorys[:, 1] = amotorys.sum(axis=1)
            # amotorys = amotorys[amotorys[:, 1] < (rangey / 2), :]

            ax.plot(amotorxs.T, amotorys.T, '.-', color='c', MarkerSize=2)  #, alpha=0.2)
        if not pmotorData.empty:
            pmotorxs = pmotorData.loc[pmotorData.t == time][['x0', 'x1']].values
            pmotorxs[:, 1] = pmotorxs.sum(axis=1)
            # pmotorxs = pmotorxs[pmotorxs[:, 1] < (rangex / 2), :]

            pmotorys = pmotorData.loc[pmotorData.t == time][['y0', 'y1']].values
            pmotorys[:, 1] = pmotorys.sum(axis=1)
            # pmotorys = pmotorys[pmotorys[:, 1] < (rangey / 2), :]

            ax.plot(pmotorxs.T, pmotorys.T, '.-', color='r', MarkerSize=2)  #, alpha=0.2)

        ax.set_aspect('equal')
        ax.set_xlim([-rangex / 2, rangex / 2])
        ax.set_ylim([-rangey / 2, rangey / 2])
        plt.tight_layout()

        if savepath:
            fig.savefig(os.path.join(savepath, 't{n:0.1f}.png'.format(n=time)))
            plt.close(fig)
    # End loop over all times
# End filaments
