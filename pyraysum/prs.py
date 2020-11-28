# Copyright 2020 Pascal Audet

# This file is part of Telewavesim.

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

Utility functions to interact with ``telewavesim`` modules.

'''
import numpy as np
from numpy import sin, cos
from scipy.signal import hilbert
from obspy.core import Trace, Stream
from pyraysum import elast
from pyraysum import tensor
import subprocess

MINERALS = ['atg', 'bt', 'cpx', 'dol', 'ep', 'grt', 'gln', 'hbl', 'jade',
            'lws', 'lz', 'ms', 'ol', 'opx', 'plag', 'qtz', 'zo']

ROCKS = ['BS_f', 'BS_m', 'EC_f', 'EC_m', 'HB', 'LHZ', 'SP_37', 'SP_80']


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

    def __init__(self, thickn, rho, vp, vs, isoflg='iso',
                 ani=None, trend=None, plunge=None, strike=None, dip=None):
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

    def update_tensor(self):
        """
        Update the elastic thickness tensor ``a``.

        Needs to be called when model parameters change.
        """
        self.nlay = len(self.thickn)
        self.a = np.zeros((3, 3, 3, 3, self.nlay))

        for j in range(self.nlay):
            if self.isoflg[j] == 1:
                cc = tensor.set_iso_tensor(self.vp[j], self.vs[j])
                self.a[:, :, :, :, j] = cc
            elif self.isoflg[j] == 0:
                cc = tensor.set_tri_tensor(self.vp[j], self.vs[j],
                                           self.trend[j], self.plunge[j],
                                           self.ani[j])
                self.a[:, :, :, :, j] = cc
            elif self.isoflg[j] in MINERALS or self.isoflg[j] in ROCKS:
                cc, rho = tensor.set_aniso_tensor(self.trend[j],
                                                  self.plunge[j],
                                                  typ=self.isoflg[j])
                self.a[:, :, :, :, j] = cc
                self.rho[j] = rho
            else:
                msg = ('\nFlag not defined: use either "iso", "tri" or one '
                       'among\n%s\n%s\n')
                raise ValueError(msg % (MINERALS, ROCKS))

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

    file = open("raysum-params", "w")
    file.writelines("# Verbosity\n" + str(int(verbose)) + "\n")
    file.writelines("# Phase name\n" + wvtype + "\n")
    file.writelines("# Multiples: 0 for none, 1 for Moho, " +
                    "2 all first-order\n" + str(mults) + "\n")
    file.writelines("# Number of samples per trace\n" + str(npts) + "\n")
    file.writelines("# Sample rate (seconds)\n" + str(dt) + "\n")
    file.writelines("# Gaussian pulse width (seconds)\n" + str(gwidth) + "\n")
    file.writelines("# Alignment: 0 is none, 1 aligns on P\n" +
                    str(align) + "\n")
    file.writelines("# Shift or traces (seconds)\n" + str(shift) + "\n")
    file.writelines("# Rotation to output: 0 is NEZ, 1 is RTZ, 2 is PVH\n" +
                    str(rot) + "\n")
    file.close()
    return


def write_geom(baz, slow):

    if not hasattr(baz, "__len__"):
        baz = [baz]
    if not hasattr(slow, "__len__"):
        slow = [slow]

    file = open("sample.geom", "w")
    dat = [(bb, ss) for ss in slow for bb in baz ]
    for dd in dat:
        file.writelines([str(dd[0]) + " " + str(dd[1]*1.e-3) + " 0. 0.\n"])
    file.close()
    return dat


def run_prs(model, verbose=False, wvtype='P', mults=2,
            npts=3000, dt=0.01, gwidth=1., align=1, shift=5., rot=0,
            baz=[], slow=[]):
    """
    """

    write_params(verbose, wvtype, mults, npts, dt, gwidth, align, shift, rot)
    dat = write_geom(baz, slow)

    subprocess.call(["seis-spread", "sample.mod",
                     "sample.geom", "sample.ph",
                     "sample.arr", "sample.tr"])



#     yx, yy, yz = raysum.get_arrivals(
#         npts, model.nlay, np.array(wvtype, dtype='c'))

#     # Transfer displacement seismograms to an ``obspy`` ``Stream`` object.
#     trxyz = get_trxyz(yx, yy, yz, npts, dt, slow, baz, wvtype)

#     return trxyz


# def get_trxyz(yx, yy, yz, npts, dt, slow, baz, wvtype):
#     """
#     Function to store displacement seismograms into ``obspy.Trace`` objects
#     and then an ``obspy`` ``Stream`` object.

#     Args:
#         ux (np.ndarray): x-component displacement seismogram
#         uy (np.ndarray): y-component displacement seismogram
#         uz (np.ndarray): z-component displacement seismogram

#     Returns:
#         (obspy.stream): trxyz: Stream containing 3-component displacement
#           seismograms

#     """

#     # Get displacements in time domain
#     ux = np.real(fft(yx))
#     uy = np.real(fft(yy))
#     uz = -np.real(fft(yz))

#     # Store in traces
#     tux = Trace(data=ux)
#     tuy = Trace(data=uy)
#     tuz = Trace(data=uz)

#     # Update trace header
#     tux = update_stats(tux, dt, slow, baz, wvtype, 'N')
#     tuy = update_stats(tuy, dt, slow, baz, wvtype, 'E')
#     tuz = update_stats(tuz, dt, slow, baz, wvtype, 'Z')

#     # Append to stream
#     trxyz = Stream(traces=[tux, tuy, tuz])

#     return trxyz


# def tf_from_xyz(trxyz, pvh=False, vp=None, vs=None):
#     """
#     Function to generate transfer functions from displacement traces.

#     Args:
#         trxyz (obspy.stream):
#             Obspy ``Stream`` object in cartesian coordinate system
#         pvh (bool, optional):
#             Whether to rotate from Z-R-T coordinate system to P-SV-SH wave mode
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


def update_stats(tr, dt, slow, baz, wvtype, cha):
    """
    Updates the ``stats`` doctionary from an obspy ``Trace`` object.

    Args:
        tr (obspy.trace): Trace object to update
        nt (int): Number of samples
        dt (float): Sampling rate
        slow (float): Slowness value (s/km)
        baz (float): Back-azimuth value (degree)

    Returns:
        (obspy.trace): tr: Trace with updated stats
    """

    tr.stats.delta = dt
    tr.stats.slow = slow
    tr.stats.baz = baz
    tr.stats.wvtype = wvtype
    tr.stats.channel = cha

    return tr


def stack_all(st1, st2, pws=False):
    """
    Stacks all traces in two ``Stream`` objects.

    Args:
        st1 (obspy.stream): Stream 1
        st2 (obspy.stream,): Stream 2
        pws (bool, optional): Enables Phase-Weighted Stacking

    Returns:
        (tuple): tuple containing:

            * stack1 (obspy.trace): Stacked trace for Stream 1
            * stack2 (obspy.trace): Stacked trace for Stream 2

    """

    print()
    print('Stacking ALL traces in streams')

    # Copy stats from stream
    str_stats = st1[0].stats

    # Initialize arrays
    tmp1 = np.zeros(len(st1[0].data))
    tmp2 = np.zeros(len(st2[0].data))
    weight1 = np.zeros(len(st1[0].data), dtype=complex)
    weight2 = np.zeros(len(st2[0].data), dtype=complex)

    # Stack all traces
    for tr in st1:
        tmp1 += tr.data
        hilb1 = hilbert(tr.data)
        phase1 = np.arctan2(hilb1.imag, hilb1.real)
        weight1 += np.exp(1j*phase1)

    for tr in st2:
        tmp2 += tr.data
        hilb2 = hilbert(tr.data)
        phase2 = np.arctan2(hilb2.imag, hilb2.real)
        weight2 += np.exp(1j*phase2)

    # Normalize
    tmp1 = tmp1/np.float(len(st1))
    tmp2 = tmp2/np.float(len(st2))

    # Phase-weighting
    if pws:
        weight1 = weight1/np.float(len(st1))
        weight2 = weight2/np.float(len(st2))
        weight1 = np.real(abs(weight1))
        weight2 = np.real(abs(weight2))
    else:
        weight1 = np.ones(len(st1[0].data))
        weight2 = np.ones(len(st1[0].data))

    # Put back into traces
    stack1 = Trace(data=weight1*tmp1, header=str_stats)
    stack2 = Trace(data=weight2*tmp2, header=str_stats)

    return stack1, stack2
