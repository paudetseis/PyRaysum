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

class Model(object):
    """
    ``model parameters``:
        - thickn (np.ndarray): Thickness of layers (km) (shape ``(nlay)``)
        - rho (np.ndarray): Density (kg/m^3) (shape ``(nlay)``)
        - vp (np.ndarray): P-wave velocity (km/s) (shape ``(nlay)``)
        - vs (np.ndarray): S-wave velocity (km/s) (shape ``(nlay)``)
        - isoflg (list of str, optional, defaut: ``'iso'``):
            Flags for type of layer material (dimension ``nlay``)
        - ani (np.ndarray, optional): Anisotropy (percent) (shape ``(nlay)``)
        - tr (np.ndarray, optional):
            Trend of symmetry axis (degree) (shape ``(nlay)``)
        - pl (np.ndarray, optional):
            Plunge of symmetry axis (degree) (shape ``(nlay)``)

        - nlay (int): Number of layers
        - a (np.ndarray): Elastic thickness (shape ``(3, 3, 3, 3, nlay)``)
    """

    def __init__(self, thickn, rho, vp, vs, isoflg='iso',
                 ani=None, tr=None, pl=None, strike=None, dip=None):
        def _get_val(v):
            return (np.array([v] * self.nlay if isinstance(v, (int, float))
                             else v) if v is not None else None)
        self.nlay = len(thickn)
        self.thickn = np.array(thickn)
        self.rho = np.array(rho) if rho is not None else [None] * self.nlay
        self.vp = np.array(vp)
        self.vs = np.array(vs)
        self.isoflg = (list(isoflg) if not isinstance(isoflg, str)
                       else [isoflg] * self.nlay)
        self.ani = _get_val(ani)
        self.tr = _get_val(tr)
        self.pl = _get_val(pl)
        self.strike = _get_val(strike)
        self.dip = _get_val(dip)
        self.update_tensor()

    def update_tensor(self):
        """
        Update the elastic thickness tensor ``a``.

        Needs to be called when model parameters change.
        """
        self.nlay = len(self.thickn)
        self.a = np.zeros((3, 3, 3, 3, self.nlay))

        for j in range(self.nlay):
            if self.isoflg[j] == 'iso':
                cc = tensor.set_iso_tensor(self.vp[j], self.vs[j])
                self.a[:, :, :, :, j] = cc
            elif self.isoflg[j] == 'tri':
                cc = tensor.set_tri_tensor(self.vp[j], self.vs[j],
                                    self.tr[j], self.pl[j], self.ani[j])
                self.a[:, :, :, :, j] = cc
            elif self.isoflg[j] in MINERALS or self.isoflg[j] in ROCKS:
                cc, rho = tensor.set_aniso_tensor(self.tr[j], self.pl[j],
                                           typ=self.isoflg[j])
                self.a[:, :, :, :, j] = cc
                self.rho[j] = rho
            else:
                msg = ('\nFlag not defined: use either "iso", "tri" or one '
                       'among\n%s\n%s\n')
                raise ValueError(msg % (MINERALS, ROCKS))


def read_model(modfile, encoding=None):
    """
    Reads model parameters from file and returns a Model object.

    Returns:
        Model object
    """
    values = np.genfromtxt(modfile, dtype=None, encoding=encoding)
    return Model(*zip(*values))




def rotate_zrt_pvh(trZ, trR, trT, slow, vp=None, vs=None):
    """
    Rotates traces from `Z-R-T` orientation to `P-SV-SH` wave mode.

    Args:
        trZ (obspy.trace): Vertical component
        trR (obspy.trace): Radial component
        trT (obspy.trace): Transverse component
        slow (float): slowness of wave
        vp (float, optional): P-wave velocity used for rotation
        vs (float, optional): S-wave velocity used for rotation

    Returns:
        (tuple): tuple containing:

            * trP (obspy.trace): Compressional (P) wave mode
            * trV (obspy.trace): Vertically polarized shear (SV) wave mode
            * trH (obspy.trace): Horizontally polarized shear (SH) wave mode

    """
    if vp is None:
        vp = 6.0
    if vs is None:
        vs = 3.5
    # Copy traces
    trP = trZ.copy()
    trV = trR.copy()
    trH = trT.copy()

    # Vertical slownesses
    qp = np.sqrt(1/vp**2 - slow**2)
    qs = np.sqrt(1/vs**2 - slow**2)

    # Elements of rotation matrix
    m11 = slow*vs*vs/vp
    m12 = -(1 - 2*vs*vs*slow*slow)/(2*vp*qp)
    m21 = (1 - 2*vs*vs*slow*slow)/(2*vs*qs)
    m22 = slow*vs

    # Rotation matrix
    rot = np.array([[-m11, m12], [-m21, m22]])

    # Vector of Radial and Vertical
    r_z = np.array([trR.data, trZ.data])

    # Rotation
    vec = np.dot(rot, r_z)

    # Extract P and SV components
    trP.data = vec[0, :]
    trV.data = vec[1, :]
    trH.data = -trT.data/2.

    return trP, trV, trH


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


def calc_ttime(model, slow, wvtype='P'):
    """
    Calculates total propagation time through model given the corresponding
    ``'P'`` or ``'S'`` wave type. The bottom layer is irrelevant in this
    calculation. All ``'S'`` wave types will return the same predicted time.
    This function is useful mainly for plotting purposes. For example, to show
    the first phase arrival at a time of 0, the traces can be shifted by the
    total propagation time through the model.

    Args:
        model (Model): Model object
        slow (float): Slowness value (s/km)
        wvtype (str): Incident wavetype (``'P'``, ``'SV'``, ``'SH'``, ``'Si'``)

    Returns:
        (float): t1: Time in seconds

    Example
    -------
    >>> from telewavesim import utils
    >>> # Define two-layer model with identical material
    >>> model = utils.Model([10, 0], None, 0, 0, 'atg', 0, 0, 0)
    >>> # Only topmost layer is useful for travel time calculation
    >>> wvtype = 'P'
    >>> slow = 0.06     # s/km
    >>> utils.calc_ttime(model, slow, wvtype)
    1.3519981570791182

    """

    t1 = 0.

    for i in range(model.nlay-1):
        if model.isoflg[i] == 'iso':
            a0 = model.a[2, 2, 2, 2, i]
            b0 = model.a[1, 2, 1, 2, i]
        else:
            cc = tensor.cc2voigt(model.a[:, :, :, :, i])
            rho = model.rho[i]
            K1, G1, K2, G2, K, G = tensor.VRH_average(cc*rho)
            a0, b0 = tensor.mod2vel(K, G, rho)
            a0 = a0**2
            b0 = b0**2
        if wvtype == 'P':
            t1 += 1000*model.thickn[i]*np.sqrt(1./a0 - (slow*1.e-3)**2)
        elif wvtype == 'Si' or wvtype == 'SV' or wvtype == 'SH':
            t1 += 1000*model.thickn[i]*np.sqrt(1./b0 - (slow*1.e-3)**2)
        else:
            raise ValueError('Invalid wave type')

    return t1



def run_prs(model, baz, slow, npts, dt, wvtype='P'):
    """
    """


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


# def update_stats(tr, dt, slow, baz, wvtype, cha):
#     """
#     Updates the ``stats`` doctionary from an obspy ``Trace`` object.

#     Args:
#         tr (obspy.trace): Trace object to update
#         nt (int): Number of samples
#         dt (float): Sampling rate
#         slow (float): Slowness value (s/km)
#         baz (float): Back-azimuth value (degree)

#     Returns:
#         (obspy.trace): tr: Trace with updated stats
#     """

#     tr.stats.delta = dt
#     tr.stats.slow = slow
#     tr.stats.baz = baz
#     tr.stats.wvtype = wvtype
#     tr.stats.channel = cha

#     return tr
