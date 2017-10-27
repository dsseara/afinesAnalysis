"""
Module to plot outputs of afines
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation


def animateQuiver(xx, yy, uut, vvt):
    '''
    creates an animated quiver plot of velocities
    '''
    fig, ax = plt.subplots()
    Q = ax.quiver(xx, yy, uut[0, ...], vvt[0, ...], pivot='mid')
    anim = animation.FuncAnimation(fig, _update_quiver, frames=uut.shape[0],
                                   fargs=(Q, uut, vvt), blit=False)

    plt.show()


def plotFilaments(xyid, configs, dt=1, dfilament=1, savepath=False):
    """Takes output from readData and plots the output as a series of .pngs
    Only works for actin data

    Parameters
    ----------
    xyid : array_like
        3D array output of readData
    configs : dict
        dictionary of afines configuration from readConfigs
    dt : scalar
        plot every dt-th frame
    dfilament : scalar
        plot every dfilament-th filament
    savepath : string (optional)
        path to save series of .pngs. Creates subfolder savePath/imgSeq. If not
        specified, does not save


    Returns
    -------
    saves series of .pngs of simulation data

    See also
    --------
    readConfigs()
    """
    if savepath:
        if not os.path.exists(os.path.join(savepath, 'imgSeq')):
            os.mkdir(os.path.join(savepath, 'imgSeq'))

    domain_xrange = configs['xrange']
    domain_yrange = configs['yrange']

    for frame, pos in enumerate(xyid[::dt, ...]):
        fig, ax = plt.subplots()
        ax.set_xlim(-np.floor(domain_xrange / 2), np.ceil(domain_xrange / 2))
        ax.set_ylim(-np.floor(domain_yrange / 2), np.ceil(domain_yrange / 2))

        for actinID in np.unique(pos[..., 2])[::dfilament]:
            actin = pos[pos[..., 2] == actinID][:, :2]
            dists = np.linalg.norm(np.diff(actin, 1, 0), axis=1)
            bools = np.reshape(dists > np.mean((domain_xrange, domain_yrange)) * 0.9,
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


def _update_quiver(frame, Q, uut, vvt):
    '''
    Helper function for animateQuiver
    '''
    u = uut[frame, ...]
    v = vvt[frame, ...]
    Q.set_UVC(u,v)
    return Q,
