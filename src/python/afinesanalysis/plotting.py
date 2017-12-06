"""
Module to plot outputs of afines
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
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


def filaments(data, configs, dt=1, dfilament=1, savepath=False):
for index, time in enumerate(pd.unique(filaments.t)):
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    data = filaments.loc[filaments.t==time]
    colors = sns.color_palette('viridis', len(pd.unique(data.clusterID)))
    #print(colors)
    for ID, data2 in data.groupby('filamentID'):
        cid = pd.unique(data2.clusterID)[0]
        #print(cid)
        actin = data2[['x','y']].values
        dists = np.linalg.norm(np.diff(actin, 1, 0), axis=1)
        bools = np.reshape(dists > np.mean((domain_xrange, domain_yrange)) * 0.9,
                           [dists.size, 1])
        mask = np.vstack([np.hstack([bools, bools]), [False, False]])
        masked_actin = np.ma.MaskedArray(actin, mask)
        if cid == -1:
            ax.plot(masked_actin[:, 0], masked_actin[:, 1], color='k')
        else:
            ax.plot(masked_actin[:, 0], masked_actin[:, 1], color=colors[cid])
    fig.savefig('coloredCluster{0}.png'.format(index))
    plt.close('all')


def plotFilaments(xyid, configs, dt=1, dfilament=1, savepath=False):
>>>>>>> origin/master
    """Takes output from readData and plots the output as a series of .pngs
    Only works for actin data

    Parameters
    ----------
    data : array_like
        Output of readData, shape (nframes, nbeads, 4). Last index refers to
        x position, y position, link_length, and filament_id, respectively
    configs : dict
        dictionary of afines configuration from readConfigs
    dt : scalar
        plot every dt-th frame
    dfilament : scalar
        plot every dfilament-th filament
    savepath : bool or string, optional
        path to save series of .pngs. Creates subfolder savePath/imgSeq. If not
        specified, does not save. If True, uses current directory

    Returns
    -------
    saves series of .pngs of simulation data

    See also
    --------
    afinesanalysis.readConfigs()
    """
    if savepath is True:
        savepath = os.curdir

    if savepath:
        if not os.path.exists(os.path.join(savepath, 'imgSeq')):
            os.mkdir(os.path.join(savepath, 'imgSeq'))

    rangex = configs['xrange']
    rangey = configs['yrange']

    for frame, pos in enumerate(data[::dt, ...]):
        fig, ax = plt.subplots()
        ax.set_xlim(-np.floor(rangex / 2), np.ceil(rangex / 2))
        ax.set_ylim(-np.floor(rangey / 2), np.ceil(rangey / 2))

        for actinID in np.unique(pos[..., -1])[::dfilament]:
            actin = pos[pos[..., -1] == actinID][:, :2]
            dists = np.linalg.norm(np.diff(actin, 1, 0), axis=1)
            bools = np.reshape(dists > np.mean((rangex, rangey)) * 0.9,
                               [dists.size, 1])
            mask = np.vstack([np.hstack([bools, bools]), [False, False]])
            masked_actin = np.ma.MaskedArray(actin, mask)
            # if any(np.abs(np.diff(actin))) > domainrange / 2:

            # dx
            ax.plot(masked_actin[:, 0], masked_actin[:, 1], 'm')

        ax.set_aspect('equal')
        if savepath:
            fig.savefig(os.path.join(savepath, 'imgSeq',
                                     'frame{number}.png'.format(number=frame)))
            plt.close(fig)


def motors(data, configs, dt=1, dfilament=1, savepath=False):
    """Takes output from readData and plots the output as a series of .pngs
    Only works for actin data

    Parameters
    ----------
    data : array_like
        Output of readData, shape (nframes, nbeads, 4). Last index refers to
        x position, y position, link_length, and filament_id, respectively
    configs : dict
        dictionary of afines configuration from readConfigs
    dt : scalar
        plot every dt-th frame
    dfilament : scalar
        plot every dfilament-th filament
    savepath : bool or string, optional
        path to save series of .pngs. Creates subfolder savePath/imgSeq. If not
        specified, does not save. If True, uses current directory

    Returns
    -------
    saves series of .pngs of simulation data

    See also
    --------
    readConfigs()
    """
