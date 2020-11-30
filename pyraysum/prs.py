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

Functions to interact with ``Raysum`` software

'''
import subprocess
import numpy as np
import pandas as pd
from obspy import Trace, Stream, UTCDateTime
from obspy.core import AttribDict


class Model(object):
    """
    ``model parameters``:
        - thickn (np.ndarray): Thickness of layers (m) (shape ``(nlay)``)
        - rho (np.ndarray): Density (kg/m^3) (shape ``(nlay)``)
        - vp (np.ndarray): P-wave velocity (m/s) (shape ``(nlay)``)
        - vs (np.ndarray): S-wave velocity (m/s) (shape ``(nlay)``)
        - isoflg (list of str, optional, defaut: ``'iso'``):
            Flags for type of layer material (dimension ``nlay``)
        - ani (np.ndarray, optional): Anisotropy (percent) (shape ``(nlay)``)
        - trend (np.ndarray, optional):
            Trend of symmetry axis (degree) (shape ``(nlay)``)
        - plunge (np.ndarray, optional):
            Plunge of symmetry axis (degree) (shape ``(nlay)``)
        - strike (np.ndarray, optional):
            azimuth of interface in RHR (degree) (shape ``(nlay)``)
        - dip (np.ndarray, optional):
            dip of interface in RHR (degree) (shape ``(nlay)``)

        - nlay (int): Number of layers
        - a (np.ndarray): Elastic thickness (shape ``(3, 3, 3, 3, nlay)``)
    """

    def __init__(self, thickn, rho, vp, vs, isoflg=1,
                 ani=None, trend=None, plunge=None,
                 strike=None, dip=None):

        def _get_val(v):
            return (np.array([v] * self.nlay if isinstance(v, (int, float))
                             else v) if v is not None else [0.]*self.nlay)
        self.nlay = len(thickn)
        self.thickn = np.array(thickn)
        self.rho = np.array(rho) if rho is not None else [None] * self.nlay
        self.vp = np.array(vp)
        self.vs = np.array(vs)
        self.isoflg = (list(isoflg) if not isinstance(isoflg, int)
                       else [isoflg] * self.nlay)
        self.ani = _get_val(ani)
        self.trend = _get_val(trend)
        self.plunge = _get_val(plunge)
        self.strike = _get_val(strike)
        self.dip = _get_val(dip)
        self.write_model()
        # self.update_tensor()

    # def update_tensor(self):
    #     """
    #     Update the elastic thickness tensor ``a``.

    #     Needs to be called when model parameters change.
    #     """
    #     self.nlay = len(self.thickn)
    #     self.a = np.zeros((3, 3, 3, 3, self.nlay))

    #     for j in range(self.nlay):
    #         if self.isoflg[j] == 1:
    #             cc = tensor.set_iso_tensor(self.vp[j], self.vs[j])
    #             self.a[:, :, :, :, j] = cc
    #         elif self.isoflg[j] == 0:
    #             cc = tensor.set_tri_tensor(self.vp[j], self.vs[j],
    #                                        self.trend[j], self.plunge[j],
    #                                        self.ani[j])
    #             self.a[:, :, :, :, j] = cc
    #         elif self.isoflg[j] in MINERALS or self.isoflg[j] in ROCKS:
    #             cc, rho = tensor.set_aniso_tensor(self.trend[j],
    #                                               self.plunge[j],
    #                                               typ=self.isoflg[j])
    #             self.a[:, :, :, :, j] = cc
    #             self.rho[j] = rho
    #         else:
    #             msg = ('\nFlag not defined: use either "iso", "tri" or one '
    #                    'among\n%s\n%s\n')
    #             raise ValueError(msg % (MINERALS, ROCKS))

    def write_model(self):
        """
        Writes model parameters to file to be processed by raysum

        """
        file = open('sample.mod', 'w')
        for i in range(self.nlay):
            file.writelines([
                str(self.thickn[i])+" "+str(self.rho[i])+" " +
                str(self.vp[i])+" "+str(self.vs[i])+" " +
                str(self.isoflg[i])+" "+str(self.ani[i])+" " +
                str(self.trend[i])+" "+str(self.plunge[i])+" " +
                str(self.strike[i])+" "+str(self.dip[i])+"\n"])
        file.close()


def read_model(modfile, encoding=None):
    """
    Reads model parameters from file and returns a Model object.

    Returns:
        Model object
    """
    values = np.genfromtxt(modfile, dtype=None, encoding=encoding)
    return Model(*zip(*values))


def write_params(verbose, wvtype, mults, npts, dt, gwidth, align, shift, rot):
    """
    Write parameters to `raysum-params` used by Raysum

    """

    file = open("raysum-params", "w")
    file.writelines("# Verbosity\n " + str(int(verbose)) + "\n")
    file.writelines("# Phase name\n " + wvtype + "\n")
    file.writelines("# Multiples: 0 for none, 1 for Moho, " +
                    "2 all first-order\n " + str(mults) + "\n")
    file.writelines("# Number of samples per trace\n " + str(npts) + "\n")
    file.writelines("# Sample rate (seconds)\n " + str(dt) + "\n")
    file.writelines("# Gaussian pulse width (seconds)\n " + str(gwidth) + "\n")
    file.writelines("# Alignment: 0 is none, 1 aligns on P\n " +
                    str(align) + "\n")
    file.writelines("# Shift or traces (seconds)\n " + str(shift) + "\n")
    file.writelines("# Rotation to output: 0 is NEZ, 1 is RTZ, 2 is PVH\n " +
                    str(rot) + "\n")
    file.close()

    return


def write_geom(baz, slow):

    # Check whether arguments are array-like; if not, store them in list
    if not hasattr(baz, "__len__"):
        baz = [baz]
    if not hasattr(slow, "__len__"):
        slow = [slow]

    # Write array_like objects to file to be used as input to Raysum
    file = open("sample.geom", "w")
    dat = [(bb, ss) for ss in slow for bb in baz]
    for dd in dat:
        file.writelines([str(dd[0]) + " " + str(dd[1]*1.e-3) + " 0. 0.\n"])
    file.close()

    # Return geometry to be used by `run_prs` to avoid too much i/o
    return dat


def read_traces(tracefile, dt, geom, rot):
    """
    Reads the traces produced by Raysum and stores them into a list
    of Stream objects

    Args:
        travefile (str): Name of file containing traces
        dt (float): Sample distance in seconds
        geom (np.ndarray): Array of [baz, slow] values
        rot (int): ID for rotation: 0 is NEZ, 1 is RTZ, 2 is PVH

    Returns:
        (list): streamlist: List of Stream objects

    """

    def _make_stats(net=None, sta=None, stime=None, dt=None,
                    slow=None, baz=None, wvtype=None, channel=None):
        """
        Updates the ``stats`` doctionary from an obspy ``Trace`` object.

        Args:
            tr (obspy.trace): Trace object to update
            dt (float): Sampling rate
            slow (float): Slowness value (s/km)
            baz (float): Back-azimuth value (degree)

        Returns:
            (obspy.trace): tr: Trace with updated stats
        """

        stats = AttribDict()
        stats.baz = baz
        stats.slow = slow
        stats.station = sta
        stats.network = net
        stats.starttime = stime
        stats.delta = dt
        stats.channel = channel
        stats.wvtype = wvtype

        return stats

    # Read traces from file
    try:
        df = pd.read_csv(tracefile)
    except:
        raise(Exception("Can't read "+str(tracefile)))

    # Component names
    if rot == 0:
        component = ['N', 'E', 'Z']
    elif rot == 1:
        component = ['R', 'T', 'Z']
    elif rot == 2:
        component = ['P', 'V', 'H']
    else:
        raise(Exception('invalid "rot" value: not in 0, 1, 2'))

    ntr = np.max(df.itr) + 1

    streamlist = []

    for itr in range(ntr):

        # Split by trace ID
        ddf = df[df.itr == itr]

        # Store into trace by channel with stats information
        # Channel 1
        stats = _make_stats(net='', sta='synt', stime=UTCDateTime(),
            dt=dt, slow=geom[itr][1], baz=geom[itr][0],
            channel='BH'+component[0])
        tr1 = Trace(data=ddf.trace1.values, header=stats)

        # Channel 2
        stats = _make_stats(net='', sta='synt', stime=UTCDateTime(),
            dt=dt, slow=geom[itr][1], baz=geom[itr][0],
            channel='BH'+component[1])
        tr2 = Trace(data=ddf.trace2.values, header=stats)

        # Channel 3
        stats = _make_stats(net='', sta='synt', stime=UTCDateTime(),
            dt=dt, slow=geom[itr][1], baz=geom[itr][0],
            channel='BH'+component[2])
        tr3 = Trace(data=ddf.trace3.values, header=stats)

        # Store into Stream object and append to list
        stream = Stream(traces=[tr1, tr2, tr3])
        streamlist.append(stream)

    return streamlist


def run_prs(model, verbose=False, wvtype='P', mults=2,
            npts=300, dt=0.1, gwidth=0.5, align=1, shift=0., rot=0,
            baz=[], slow=[]):
    """
    Reads the traces produced by Raysum and stores them into a list
    of Stream objects

    Args:
        model (:class:`~pyraysum.prs.Model`): Seismic velocity model
        verbose (bool): Whether or not to increase verbosity of Raysum
        wvtype (str): Wave type of incoming wavefield ('P', 'SV', or 'SH')
        mults (int): ID for calculating free surface multiples
                ('0': no multiples, '1': Moho only, '2': all first-order)
        npts (int): Number of samples in time series
        dt (float): Sampling distance in seconds
        gwidth (float): Width of Gaussian pulse in seconds
        align (int): ID for alignment of seismograms ('1': align at 'P',
                '2': align at 'SV' or 'SH')
        shift (float): Time shift in seconds (positive shift moves seismograms
                to greater lags)
        rot (int): ID for rotation: 0 is NEZ, 1 is RTZ, 2 is PVH
        baz (array_like): Array of input back-azimuth values in degrees
        slow (array_like): Array of input slowness values to model in s/km

    Returns:
        (list): streamlist: List of Stream objects

    """

    # Write parameter file to be used by Raysum
    write_params(verbose, wvtype, mults, npts, dt, gwidth, align, shift, rot)

    # Write geometry (baz, slow) to be used by Raysum
    geom = write_geom(baz, slow)

    # Call Raysum to produce the output 'sample.tr' containing synthetic traces
    subprocess.call(["seis-spread", "sample.mod",
                     "sample.geom", "sample.ph",
                     "sample.arr", "sample.tr"])

    # Read all traces and store them into a list of :class:`~obspy.core.Stream`
    # objects
    streamlist = read_traces('sample.tr', dt, geom, rot)

    return streamlist


# def tf_from_xyz(trxyz, pvh=False, vp=None, vs=None):
#     """
#     Function to generate transfer functions from displacement traces.

#     Args:
#         trxyz (obspy.stream):
#             Obspy ``Stream`` object in cartesian coordinate system
#         pvh (bool, optional):
#             Whether to rotate from Z-R-T coordinate system to
#             P-SV-SH wave mode
#         vp (float, optional):
#             Vp velocity at surface for rotation to P-SV-SH system
#         vs (float, optional):
#             Vs velocity at surface for rotation to P-SV-SH system

#     Returns:
#         (obspy.stream):
#             tfs: Stream containing Radial and Transverse transfer functions

#     """

#     # Extract East, North and Vertical
#     ntr = trxyz[0]
#     etr = trxyz[1]
#     ztr = trxyz[2]
#     baz = ntr.stats.baz
#     slow = ntr.stats.slow
#     wvtype = ntr.stats.wvtype

#     # Copy to radial and transverse
#     rtr = ntr.copy()
#     ttr = etr.copy()

#     # Rotate to radial and transverse
#     rtr.data, ttr.data = rotate_ne_rt(ntr.data, etr.data, baz)

#     if pvh:
#         trP, trV, trH = rotate_zrt_pvh(ztr, rtr, ttr, slow, vp=vp, vs=vs)

#         tfr = trV.copy()
#         tfr.data = np.zeros(len(tfr.data))
#         tft = trH.copy()
#         tft.data = np.zeros(len(tft.data))
#         ftfv = fft(trV.data)
#         ftfh = fft(trH.data)
#         ftfp = fft(trP.data)

#         if wvtype == 'P':
#             # Transfer function
#             tfr.data = fftshift(np.real(ifft(np.divide(ftfv, ftfp))))
#             tft.data = fftshift(np.real(ifft(np.divide(ftfh, ftfp))))
#         elif wvtype == 'Si':
#             tfr.data = fftshift(np.real(ifft(np.divide(-ftfp, ftfv))))
#             tft.data = fftshift(np.real(ifft(np.divide(-ftfp, ftfh))))
#         elif wvtype == 'SV':
#             tfr.data = fftshift(np.real(ifft(np.divide(-ftfp, ftfv))))
#         elif wvtype == 'SH':
#             tft.data = fftshift(np.real(ifft(np.divide(-ftfp, ftfh))))
#     else:
#         tfr = rtr.copy()
#         tfr.data = np.zeros(len(tfr.data))
#         tft = ttr.copy()
#         tft.data = np.zeros(len(tft.data))
#         ftfr = fft(rtr.data)
#         ftft = fft(ttr.data)
#         ftfz = fft(ztr.data)

#         if wvtype == 'P':
#             # Transfer function
#             tfr.data = fftshift(np.real(ifft(np.divide(ftfr, ftfz))))
#             tft.data = fftshift(np.real(ifft(np.divide(ftft, ftfz))))
#         elif wvtype == 'Si':
#             tfr.data = fftshift(np.real(ifft(np.divide(-ftfz, ftfr))))
#             tft.data = fftshift(np.real(ifft(np.divide(-ftfz, ftft))))
#         elif wvtype == 'SV':
#             tfr.data = fftshift(np.real(ifft(np.divide(-ftfz, ftfr))))
#         elif wvtype == 'SH':
#             tft.data = fftshift(np.real(ifft(np.divide(-ftfz, ftft))))

#     # Store in stream
#     tfs = Stream(traces=[tfr, tft])

#     # Return stream
#     return tfs
