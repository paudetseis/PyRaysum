# Copyright 2020 Pascal Audet

# This file is part of PyRaysum.

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

'''

Functions to plot single station displacement or receiver function seismograms.

'''

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import hilbert
from obspy import Stream, Trace


def stack_all(stream, pws=False):
    """
    Stacks all traces in ``Stream`` objects.

    Args:
        stream (:class:`~obspy.core.Stream`)
            Contains traces to stack
        pws (bool):
            Enables Phase-Weighted Stacking

    Returns:
        (tuple): tuple containing:

            stack (:class:`~obspy.core.Trace`):
                Stacked trace

    """

    # Copy stats from stream
    str_stats = stream[0].stats

    # Initialize arrays
    tmp = np.zeros(len(stream[0].data))
    weight = np.zeros(len(stream[0].data), dtype=complex)

    # Stack all traces
    for tr in stream:
        tmp += tr.data
        hilb = hilbert(tr.data)
        phase = np.arctan2(hilb.imag, hilb.real)
        weight += np.exp(1j*phase)

    # Normalize
    tmp = tmp/np.float(len(stream))

    # Phase-weighting
    if pws:
        weight = weight/np.float(len(stream))
        weight = np.real(abs(weight))
    else:
        weight = np.ones(len(stream[0].data))

    # Put back into traces
    stack = Trace(data=weight*tmp, header=str_stats)

    return stack


def rf_wiggles(rflist, btyp='baz', wvtype='P', pws=False, tmin=-5., tmax=20,
               scale=None, save=False, ftitle='Figure_rf_wiggle',
               fmt='png'):
    """
    Plots receiver function seismograms sorted by back-azimuth or slowness.

    Args:
        rflist (list):
            list of :class:`~obspy.core.Stream` objects containing receiver
            functions
        btyp (str):
            Type of sorting for panel
        wvtype (str):
            Wavet type ('P', 'SV', or 'SH')
        pws (bool):
            Enables Phase-Weighted Stacking
        tmin (float):
            Lower bound of time axis (s)
        tmax (float):
            Upper bound of time axis (s)
        scale (float):
            Scaling factor
        save (bool):
            Whether or not to save the figure
        ftitle (str):
            Title of figure to be saved
        wvtype (str):
            Wave type ('P', 'SV', or 'SH')

    Returns:
        None
    """

    if not (btyp == 'baz' or btyp == 'slow' or btyp == 'dist'):
        raise ValueError('type has to be "baz" or "slow" or "dist"')

    if not fmt in ['png', 'PNG', 'jpg', 'JPG', 'eps', 'EPS', 'pdf', 'PDF']:
        raise ValueError("'fmt' has to be one of 'png', 'jpg', 'eps', 'pdf'")

    # Re-order streams in list
    str1 = Stream(traces=[st[0] for st in rflist])
    str2 = Stream(traces=[st[1] for st in rflist])

    # Get stacked traces
    tr1 = stack_all(str1, pws=pws)
    tr2 = stack_all(str2, pws=pws)

    # Time axis
    time = str1[0].stats.taxis

    # Initialize figure
    fig = plt.figure()
    plt.clf()

    # Get more control on subplots
    ax1 = fig.add_axes([0.1, 0.825, 0.3, 0.05])
    ax2 = fig.add_axes([0.1, 0.1, 0.3, 0.7])
    ax3 = fig.add_axes([0.45, 0.825, 0.3, 0.05])
    ax4 = fig.add_axes([0.45, 0.1, 0.3, 0.7])

    # Plot stack of all traces from str1 on top left
    ax1.fill_between(time, 0., tr1.data, where=tr1.data+1e-6 <= 0.,
                     facecolor='blue', linewidth=0)
    ax1.fill_between(time, 0., tr1.data, where=tr1.data+1e-6 >= 0.,
                     facecolor='red', linewidth=0)
    ax1.set_ylim(-np.max(np.abs(tr1.data)), np.max(np.abs(tr1.data)))
    ax1.set_yticks(())
    ax1.set_xticks(())
    ax1.set_title('Radial')
    ax1.set_xlim(tmin, tmax)

    # Plot stack of all SH traces on top right
    ax3.fill_between(time, 0., tr2.data, where=tr2.data+1e-6 <= 0.,
                     facecolor='blue', linewidth=0)
    ax3.fill_between(time, 0., tr2.data, where=tr2.data+1e-6 >= 0.,
                     facecolor='red', linewidth=0)
    ax3.set_xlim(tmin, tmax)
    ax3.set_ylim(-np.max(np.abs(tr1.data)), np.max(np.abs(tr1.data)))
    ax3.set_yticks(())
    ax3.set_xticks(())
    ax3.set_title('Transverse')

    axes = [ax2, ax4]
    streams = [str1, str2]

    for ax, st in zip(axes, streams):

        # Plot sorted traces from str1 on bottom left panel
        for tr in st:

            if scale:
                maxval = scale
                # Define y axis
                if btyp == 'baz':
                    y = tr.stats.baz
                elif btyp == 'slow':
                    y = tr.stats.slow
                elif btyp == 'dist':
                    y = tr.stats.slow
            else:
                # Define y axis
                if btyp == 'baz':
                    y = tr.stats.baz
                    maxval = 180
                elif btyp == 'slow':
                    y = tr.stats.slow
                    maxval = 0.02
                elif btyp == 'dist':
                    y = tr.stats.slow
                    maxval = 20

            # Fill positive in red, negative in blue
            ax.fill_between(
                time, y, y+tr.data*maxval, where=tr.data+1e-6 <= 0.,
                facecolor='blue', linewidth=0)
            ax.fill_between(
                time, y, y+tr.data*maxval, where=tr.data+1e-6 >= 0.,
                facecolor='red', linewidth=0)

        ax.set_xlim(tmin, tmax)

        if btyp == 'baz':
            ax.set_ylim(-5, 370)
            ax.set_ylabel('Back-azimuth (deg)')

        elif btyp == 'slow':
            if wvtype == 'P':
                ax.set_ylim(0.038, 0.082)
            elif wvtype == 'S':
                ax.set_ylim(0.07, 0.125)
            elif wvtype == 'SKS':
                ax.set_ylim(0.03, 0.06)
            ax.set_ylabel('Slowness (s/km)')
        elif btyp == 'dist':
            if wvtype == 'P':
                ax.set_ylim(28., 92.)
            elif wvtype == 'S':
                ax.set_ylim(53., 107.)
            elif wvtype == 'SKS':
                ax.set_ylim(83., 117.)
            ax.set_ylabel('Distance (deg)')

        ax.set_xlabel('Time (sec)')
        ax.grid(ls=':')

    # Remove labels on axis 4
    ax4.set_ylabel('')
    ax4.set_yticklabels(())

    if save:
        plt.savefig(ftitle+'.'+fmt, dpi=300, bbox_inches='tight', format=fmt)
    else:
        plt.show()

    return


def stream_wiggles(streamlist, btyp='baz', wvtype='P', tmin=-5., tmax=20.,
                   scale=None, save=False, ftitle='Figure_pw_wiggles',
                   fmt='png'):
    """
    Plots displacement seismograms sorted by back-azimuth or slowness.

    Args:
        streamlist (list):
            list of :class:`~obspy.core.Stream` objects containing displacement
            seismograms
        btyp (str):
            Type of sorting for panel
        wvtype (str):
            Wave type ('P', 'SV', or 'SH')
        tmin (float):
            Lower bound of time axis (s)
        tmax (float):
            Upper bound of time axis (s)
        scale (float):
            Scaling factor
        save (bool):
            Whether or not to save the figure
        ftitle (str):
            Title of figure to be saved

    Returns:
        None
    """

    if not (btyp == 'baz' or btyp == 'slow' or btyp == 'dist'):
        print('type has to be "baz" or "slow" or "dist"')
        return

    # Re-order streams in list
    str1 = Stream(traces=[st[0] for st in streamlist])
    str2 = Stream(traces=[st[1] for st in streamlist])
    str3 = Stream(traces=[st[2] for st in streamlist])

    # Time axis
    time = str1[0].stats.taxis

    # Initialize figure
    fig = plt.figure()
    plt.clf()

    # Get more control on subplots
    ax1 = fig.add_axes([0.1, 0.1, 0.25, 0.83])
    ax2 = fig.add_axes([0.4, 0.1, 0.25, 0.83])
    ax3 = fig.add_axes([0.7, 0.1, 0.25, 0.83])

    axes = [ax1, ax2, ax3]
    streams = [str1, str2, str3]

    for ax, st in zip(axes, streams):

        for tr in st:
            if scale:
                maxval = scale
                # Define y axis
                if btyp == 'baz':
                    y = tr.stats.baz
                elif btyp == 'slow':
                    y = tr.stats.slow
                elif btyp == 'dist':
                    y = tr.stats.slow
            else:
                # Define y axis
                if btyp == 'baz':
                    y = tr.stats.baz
                    maxval = 100
                elif btyp == 'slow':
                    y = tr.stats.slow
                    maxval = 0.02
                elif btyp == 'dist':
                    y = tr.stats.slow
                    maxval = 20

            # Fill positive in red, negative in blue
            ax.fill_between(
                time, y, y+tr.data*maxval, where=tr.data+1e-6 <= 0.,
                facecolor='blue', linewidth=0)
            ax.fill_between(
                time, y, y+tr.data*maxval, where=tr.data+1e-6 >= 0.,
                facecolor='red', linewidth=0)

        ax.set_xlim(tmin, tmax)

        if btyp == 'baz':
            ax.set_ylim(-5, 370)
            ax.set_ylabel('Back-azimuth (deg)')

        elif btyp == 'slow':
            if wvtype == 'P':
                ax.set_ylim(0.038, 0.082)
            elif wvtype == 'S':
                ax.set_ylim(0.07, 0.125)
            elif wvtype == 'SKS':
                ax.set_ylim(0.03, 0.06)
            ax.set_ylabel('Slowness (s/km)')
        elif btyp == 'dist':
            if wvtype == 'P':
                ax.set_ylim(28., 92.)
            elif wvtype == 'S':
                ax.set_ylim(53., 107.)
            elif wvtype == 'SKS':
                ax.set_ylim(83., 117.)
            ax.set_ylabel('Distance (deg)')

        ax.set_xlabel('Time (sec)')
        ax.set_title(st[0].stats.channel)
        ax.grid(ls=':')

    # Remove labels on axes 2 and 3
    ax2.set_ylabel('')
    ax3.set_ylabel('')
    ax2.set_yticklabels(())
    ax3.set_yticklabels(())

    if save:
        plt.savefig(ftitle+'.'+fmt, dpi=300, bbox_inches='tight', format=fmt)
    else:
        plt.show()

    return


def seis_wiggles(stream, tmin=-5., tmax=20., save=False,
                 ftitle='Figure_pw_wiggles_3c', fmt='png'):
    """
    Plots 3-component wiggles.

    Args:
        stream (:class:`~obspy.core.Stream`):
            Stream containing 3 traces
        tmin (float):
            Lower bound of time axis (s)
        tmax (float, optional):
            Upper bound of time axis (s)
        save (bool):
            Whether or not to save the figure
        ftitle (str):
            Title of figure to be saved
        fmt (str):
            Format of figure to be saved

    Returns:
        None
    """

    nt = stream[0].stats.npts
    dt = stream[0].stats.delta

    # Clear figure
    fig = plt.figure()

    # Time axis
    time = stream[0].stats.taxis

    max1 = np.max(np.abs(stream[0].data))
    max2 = np.max(np.abs(stream[1].data))
    max3 = np.max(np.abs(stream[2].data))

    ax = fig.add_subplot(313)
    maxv = np.max(np.array([max1, max2, max3]))
    ax.plot(time, stream[2].data/maxv, 'k',
            label='Vertical component', lw=0.75)
    ax.set_xlim(tmin, tmax)
    ax.set_ylim(-1.1, 1.1)
    ax.set_xlabel('Time following $P$-wave arrival (sec)')

    plt.legend()

    ax = fig.add_subplot(312)
    ax.plot(time, stream[0].data/maxv, 'k',
            label='North component', lw=0.75)
    ax.set_xlim(tmin, tmax)
    ax.set_ylim(-1.1, 1.1)
    ax.set_xticklabels(())

    plt.legend()

    ax = fig.add_subplot(311)
    ax.plot(time, stream[1].data/maxv, 'k', label='East component', lw=0.75)
    ax.set_xlim(tmin, tmax)
    ax.set_ylim(-1.1, 1.1)
    ax.set_xticklabels(())

    plt.legend()

    plt.legend()
    plt.tight_layout()

    if save:
        plt.savefig(ftitle+'.'+fmt, dpi=300, bbox_inches='tight', format=fmt)
    else:
        plt.show()

    return
