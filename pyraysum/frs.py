# Copyright 2022 Wasja Bloch, Pascal Audet

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

import numpy as np
from scipy import signal
from obspy import Trace, Stream, UTCDateTime
from numpy.fft import fft, ifft, fftshift

# Here to avoid circular import
_phnames = {1: "P", 2: "S", 3: "T", 4: "p", 5: "s", 6: "t"}


def read_traces(traces, rc, geometry, arrivals=None):
    """
    Create a :class:`Result` from the array produced by :meth:`fraysum.run_bare()`
    and :meth:`fraysum.run_full`.

    Parameters:
        rc (:class:`Control`):
            Run-control parameters
        geometry (:class:`Geometry`):
            Ray parameters
        arrivals (list):
            Output of :meth:`read_arrivals`. List of arrival times, amplitudes, and
            names

    Returns:
        list:
            List of :class:`obspy.Stream` objects
    """

    npts = rc.npts
    dt = rc.dt
    rot = rc.rot
    shift = rc.shift
    ntr = geometry.ntr

    # Crop unused overhang of oversized fortran arrays
    trs = [
        traces[0, :npts, :ntr].reshape(npts * ntr, order="F"),
        traces[1, :npts, :ntr].reshape(npts * ntr, order="F"),
        traces[2, :npts, :ntr].reshape(npts * ntr, order="F"),
    ]

    itr = np.array([npts * [tr] for tr in range(ntr)]).reshape(npts * ntr)

    # Component names
    if rot == 0:
        # Rotate to seismometer convention
        component = ["N", "E", "Z"]
        order = [2, 0, 1]
    elif rot == 1:
        component = ["R", "T", "Z"]
        order = [0, 1, 2]
    elif rot == 2:
        component = ["P", "V", "H"]
        order = [0, 1, 2]
    else:
        raise (ValueError('Invalid value for "rot": Must be 0, 1, 2'))

    taxis = np.arange(npts) * dt - shift

    streams = []

    for iitr in range(ntr):

        # Split by trace ID
        istr = itr == iitr

        # Store into trace by component with stats information
        stream = Stream()
        for ic in order:
            stats = {
                "baz": geometry[iitr][0],
                "slow": geometry[iitr][1],
                "station": "",
                "network": "",
                "starttime": UTCDateTime(0),
                "delta": dt,
                "channel": "SY" + component[ic],
                "taxis": taxis,
            }

            if arrivals:
                stats.update(
                    {
                        "phase_times": arrivals[iitr][ic][0],
                        "phase_amplitudes": arrivals[iitr][ic][1],
                        "phase_descriptors": arrivals[iitr][ic][2],
                        "phase_names": arrivals[iitr][ic][3],
                        "conversion_names": arrivals[iitr][ic][4],
                    }
                )

            if rot == 0 and component[ic] != "Z":
                # Raysum has z down, change here to z up
                tr = Trace(data=-trs[ic][istr], header=stats)
                tr.stats.phase_amplitudes *= -1
            else:
                tr = Trace(data=trs[ic][istr], header=stats)

            stream.append(tr)

        # Store into Stream object and append to list
        streams.append(stream)

    return streams


def read_arrivals(ttimes, amplitudes, phaselist, geometry):
    """
    Convert the :const:`phaselist`, :const:`amplitude` and :const:`traveltime` output of
    :meth:`fraysum.run_full()` to lists of phase arrival times, amplitudes, long phase
    descriptors, short phase names, and conversion names.

    Parameters:
        ttimes (np.array):
            Travel time array ...
        amplitudes (np.array):
            Amplitude array ...
        phaselist (np.array):
            Phase identifier array returned by :func:`fraysum.run_full`
        geometry (:class:`prs.Geometry`):
            Ray geometry

    Returns:
        list:
            3-component phase arrival lists, where indices translate to:

            * :const:`0`: phase arrival times
            * :const:`1`: phase amplitudes
            * :const:`2`: (long) phase descriptors, e.g. "1P0S"
            * :const:`3`: (short) phase names, e.g. "PS"
            * :const:`4`: (intermediate) conversion names, e.g. "P1S"
    """

    dscrs = []
    phnms = []
    convs = []
    for iph in range(len(phaselist[0, 0, :])):
        dscr = ""
        phnm = ""
        conv = ""

        for iseg in range(len(phaselist[:, 0, 0])):
            if phaselist[iseg, 0, iph] == 0:
                # No more reflections / conversions
                dscrs.append(dscr)
                phnms.append(phnm)
                convs.append(conv)
                break

            phid = phaselist[iseg, 1, iph]
            layn = phaselist[iseg, 0, iph] - 1  # Python indexing

            phn = _phnames[phid]
            dscr += str(layn) + phn

            if phn.isupper():
                # Upward-conversion at top of layer below
                con = str(layn + 1) + phn
            else:
                con = str(layn) + phn

            # Omit not-converted segment from phase name and conversions
            try:
                if phnm[-1] == phn:
                    phn = ""
                    con = ""
            except IndexError:
                con = ""

            # Always prefix incoming wavetype to conversions
            if iseg == 0:
                con = phn

            phnm += phn
            conv += con

        if phaselist[0, 0, iph + 1] == 0:
            break

    nphs = len(dscrs)

    dscrs = np.array(dscrs)
    phnms = np.array(phnms)
    convs = np.array(convs)

    tanss = []
    for itr in range(geometry.ntr):
        tans = []
        for comp in range(len(amplitudes[:, 0, 0])):
            amps = amplitudes[comp, :nphs, itr]
            ia = abs(amps) > 1e-6  # Dismiss 0 amplitude arrivals
            tans.append(
                np.array(
                    [ttimes[:nphs, itr][ia], amps[ia], dscrs[ia], phnms[ia], convs[ia]],
                    dtype=object,
                )
            )
        tanss.append(tans)

    return tanss


def make_array(geometry, rc):
    """
    Initialize array for `NumPy`-based post processing

    Parameters:
        geometry (:class:`prs.Geometry`):
            Recording geometry
        rc (:class:`prs.Control`):
            Run-control parameters

    Returns:
        :const:`numpy.zeros((geometry.ntr, 2, rc.npts))`:
            Array in shape to be used by :func:`filtered_rf_array` and
            :func:`filtered_array`.
    """
    return np.zeros((geometry.ntr, 2, rc.npts))


cached_coefficients = {}


def _get_cached_bandpass_coefs(order, corners):
    # from pyrocko.Trace.filter
    ck = (order, tuple(corners))
    if ck not in cached_coefficients:
        cached_coefficients[ck] = signal.butter(order, corners, btype="band")

    return cached_coefficients[ck]


def filtered_rf_array(traces, rfarray, ntr, npts, dt, fmin, fmax):
    """
    Fast, `NumPy`-based, receiver function computation and filtering of
    :meth:`fraysum.run_bare()` output

    Roughly equivalent to subsequent calls to :func:`read_traces()`,
    :meth:`Result.calculate_rfs()` and :meth:`Result.filter()`, stripped down
    for computational efficiency, for use in inversion/probabilistic approaches.

    Parameters:
        traces (np.ndarray):
            Output of :meth:`fraysum.run_bare()`
        rfarray (np.ndarray):
            Initialized array of shape (ntr, 2, npts) to store output (See:
            :func:`make_array()`)
        ntr (int):
            Number of traces (:attr:`Geometry.ntr`)
        npts (int):
            Number of points per trace (:attr:`Control.npts`)
        dt (float):
            Sampling interval (:attr:`Control.dt`)
        fmin (float):
            Lower bandpass frequency corner (Hz)
        fmax (float):
            Upper bandpass frequency corner (Hz)

    Returns:
        None:
            None. Output is written to :const:`rfarray`

    Warning:
        Assumes PVH alignment (ray-polarization), i.e. :attr:`Control.rot="PVH"`.
    """

    order = 2

    def _bandpass(arr):
        # from pyrocko.Trace.filter
        (b, a) = _get_cached_bandpass_coefs(order, (2 * dt * fmin, 2 * dt * fmax))
        arr -= np.mean(arr)
        firstpass = signal.lfilter(b, a, arr)
        return signal.lfilter(b, a, firstpass[::-1])[::-1]

    # Crop unused overhang of oversized fortran arrays and transpose to
    # [traces[components[samples]]] order
    data = np.array(
        [traces[0, :npts, :ntr], traces[1, :npts, :ntr], traces[2, :npts, :ntr]]
    ).transpose(2, 0, 1)

    for n, trace in enumerate(data):
        ft_ztr = fft(trace[0])  # P or R or N
        ft_rfr = fft(trace[1])  # V or T or E
        ft_rft = fft(trace[2])  # H or Z or Z

        # assuming PVH:
        rfarray[n, 0, :] = _bandpass(fftshift(np.real(ifft(np.divide(ft_rfr, ft_ztr)))))
        rfarray[n, 1, :] = _bandpass(fftshift(np.real(ifft(np.divide(ft_rft, ft_ztr)))))


def filtered_array(traces, rfarray, ntr, npts, dt, fmin, fmax):
    """
    Fast, `NumPy`-based filtering of :meth:`fraysum.run_bare()` output

    Roughly equivalent to subsequent calls to :func:`read_traces()`, and
    :meth:`Result.filter()`, stripped down for computational efficiency,
    for use in inversion/probabilistic approaches.

    Parameters:
        traces (np.ndarray):
            Output of :meth:`fraysum.run_bare()`
        rfarray (np.ndarray):
            Initialized array of shape (ntr, 2, npts) to store output (See:
            :func:`make_array()`)
        ntr (int):
            Number of traces (:attr:`Geometry.ntr`)
        npts (int):
            Number of points per trace (:attr:`Control.npts`)
        dt (float):
            Sampling intervall (:attr:`Control.dt`)
        fmin (float):
            Lower bandpass frequency corner (Hz)
        fmax (float):
            Upper bandpass frequency corner (Hz)

    Returns:
        None:
            None. Output is written to :const:`rfarray`

    Warning:
        Assumes PVH alignment (ray-polarization), i.e. :attr:`Control.rot="PVH"`.
    """

    npts2 = npts // 2
    rem = npts % 2

    order = 2

    def _bandpass(arr):
        # from pyrocko.Trace.filter
        (b, a) = _get_cached_bandpass_coefs(order, (2 * dt * fmin, 2 * dt * fmax))
        arr -= np.mean(arr)
        firstpass = signal.lfilter(b, a, arr)
        return signal.lfilter(b, a, firstpass[::-1])[::-1]

    # Crop unused overhang of oversized fortran arrays and transpose to
    # [traces[components[samples]]] order
    data = np.array(
        [traces[0, :npts, :ntr], traces[1, :npts, :ntr], traces[2, :npts, :ntr]]
    ).transpose(2, 0, 1)

    for n, trace in enumerate(data):
        # assuming PVH:
        rfarray[n, 0, npts2:] = _bandpass(trace[1][: npts2 + rem])  # SV
        rfarray[n, 1, npts2:] = _bandpass(trace[2][: npts2 + rem])  # SH
