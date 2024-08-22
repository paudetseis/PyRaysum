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

import re
import numpy as np
import re
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from obspy import Stream
from numpy.fft import fft, ifft, fftshift
from copy import deepcopy
import fraysum

from pyraysum import plot
from pyraysum.frs import read_arrivals, read_traces, _phnames

_alignn = {0: "none", 1: "P", 2: "SV", 3: "SH"}  # alignment name
_aligni = {_alignn[k]: k for k in _alignn}  # alignment ID
_rotn = {0: "ZNE", 1: "RTZ", 2: "PVH"}  # rotation name
_roti = {_rotn[k]: k for k in _rotn}  # rotation ID
_iphase = {"P": 1, "SV": 2, "SH": 3}
_phids = {_phnames[k]: k for k in _phnames}  # inverse dictionary, imported in core
_modhint = (
    "################################################\n"
    "#\n"
    "#   Model file to use with `PyRaysum` for \n"
    "#   modeling teleseismic body wave propagation \n"
    "#   through dippin anisotropic media.\n"
    "#\n"
    "#   Lines starting with '#' are ignored. Each \n"
    "#   line corresponds to a unique layer. The \n"
    "#   bottom layer is assumed to be a half-space\n"
    "#   (Thickness is irrelevant).\n"
    "#\n"
    "#   Format:\n"
    "#       Column  Contents\n"
    "#          0    Thickness (km)\n"
    "#          1    Density (kg/m^3)\n"
    "#          2    Layer P-wave velocity (km/s)\n"
    "#          3    Layer S-wave velocity (km/s)\n"
    "#          4    Layer flag \n"
    "#                   1: isotropic\n"
    "#                   0: transverse isotropy\n"
    "#          5    % Transverse anisotropy (if Layer flag is set to 0)\n"
    "#                   0: isotropic\n"
    "#                   +: fast symmetry axis\n"
    "#                   -: slow symmetry axis\n"
    "#          6    Trend of symmetry axis (degrees)\n"
    "#          7    Plunge of symmetry axis (degrees)\n"
    "#		    8	 Interface strike (degrees)\n"
    "#		    9	 Interface dip (degrees)\n"
    "#\n"
    "################################################\n"
)


class Model(object):
    """Model of the subsurface seismic velocity structure.

    .. note::
       This object holds the infromation of the .mod file in classic Raysum

    Parameters:
        thickn (array_like):
          Thickness of layers (m)
        rho (float or array_like):
          Density (kg/m^3)
        vp (float or array_like):
          P-wave velocity (m/s)
        vs (float or array_like, optional):
          S-wave velocity (m/s)
          If None, computed from :const:`vpvs`
        flag (array_like of str, optional):
          :const:`1` for isotropic, :const:`0` for anisotropic
        ani (float or array_like, optional):
          Anisotropy (percent)
        trend (float or array_like, optional):
          Trend of symmetry axis (degree)
        plunge (float or array_like, optional):
          Plunge of symmetry axis (degree)
        strike (float or array_like, optional):
          azimuth of interface in RHR (degree clockwise from north)
        dip (float or array_like, optional):
          dip of interface in RHR (degree down from horizontal)
        vpvs (float or array_like, optional):
          P-to-S velocity ratio.
          Defaults to 1.73. Ignored if :const:`vs` is set.
        maxlay (int):
          Maximum number of layers defined in params.h

    Warning:
        When setting `vpvs`, `vs` is adjusted to satisfy vs = vp / vpvs.

    The following attributes are set upon initialization, when setting a layer property
    via :meth:`Model.__setitem__()` (i.e. :attr:`model[0]` or :attr:`model[0,
    "thickn"]`), and when execuding :meth:`Model.update()` after setting a single
    element of the above model attributes. `f` prefixes indicate attributes used for
    interaction with `fraysum.run_bare()` and `fraysum.run_full()`. Set these directly
    for best performance.

        nlay
          Number of layers
        parameters
          Convenience attribute that collects the `f`-attributes in the order expected
          by `fraysum.run_bare()` and `fraysum.run_full()`
        fthickn
          Thickness of layers (m)
        frho
          Density (kg/m^3)
        fvp
          P-wave velocity (m/s)
        fvs
          S-wave velocity (m/s)
        fflag
          Flag indicating isotropy (1) or anisotropy (0) of layer material
        fani
          Anisotropy (percent)
        ftrend
          Trend of symmetry axis (radians)
        fplunge
          Plunge of symmetry axis (radians)
        fstrike
          azimuth of interface in RHR (radians)
        fdip
          dip of interface in RHR (radians)

    Example
    -------
        >>> from pyraysum import Model
        >>> model = Model([10000, 0], [3000, 4500], [6000, 8000], [3500, 4600])
        >>> print(model)
        #  thickn     rho      vp      vs  flag aniso   trend plunge strike   dip
          10000.0  3000.0  6000.0  3500.0    1    0.0     0.0    0.0    0.0   0.0
              0.0  4500.0  8000.0  4600.0    1    0.0     0.0    0.0    0.0   0.0
        >>> model[0]["thickn"]
        10000.0
        >>> model.thickn[0]
        10000.0
        >>> model.thickn[0] = 15000  # model.fthickn has not yet changed
        >>> model.update()  # Now everything works as expected
        >>> model.fthickn[0]
        15000.0
        >>> model[1, "vp"]
        8000.0
        >>> model[1, "vp"] = 9000.  # Does not require model.update()
        >>> model[1, "vp"]
        9000.0
        >>> model["vp", 1] = 7400.  # Indices can be swapped
        >>> model["vp", 1]
        7400.0
        >>> model[1] = {"vs": 4200., "rho": 4000.}
        >>> model[1]["vs"]
        4200.0
        >>> model[1, "rho"]
        4000.0
        >>> model[1] = {"thickn": 5000, "vp": 8000., "vpvs": 2.}
        >>> model[1]["vs"]
        4000.0
        >>> model += [5000, 3600, 8000, 4000]  # thickn, rho, vp, vs
        >>> print(model)
        #  thickn     rho      vp      vs  flag aniso   trend plunge strike   dip
          15000.0  3000.0  6000.0  3500.0    1    0.0     0.0    0.0    0.0   0.0
           5000.0  4000.0  8000.0  4000.0    1    0.0     0.0    0.0    0.0   0.0
           5000.0  3600.0  8000.0  4000.0    1    0.0     0.0    0.0    0.0   0.0
        >>> model += {"thickn": 0, "rho": 3800, "vp": 8500., "dip": 20, "strike": 90}
        >>> print(model)
        #  thickn     rho      vp      vs  flag aniso   trend plunge strike   dip
          15000.0  3000.0  6000.0  3500.0    1    0.0     0.0    0.0    0.0   0.0
           5000.0  4000.0  8000.0  4000.0    1    0.0     0.0    0.0    0.0   0.0
           5000.0  3600.0  8000.0  4000.0    1    0.0     0.0    0.0    0.0   0.0
              0.0  3800.0  8500.0  4913.3    1    0.0     0.0    0.0   90.0  20.0
    """

    def __init__(
        self,
        thickn,
        rho,
        vp,
        vs=None,
        flag=1,
        ani=None,
        trend=None,
        plunge=None,
        strike=None,
        dip=None,
        vpvs=1.73,
        maxlay=15,
    ):
        def _array(v):
            if v is not None:
                return np.array(
                    [v] * self.nlay if isinstance(v, (int, float)) else v, dtype=float
                )
            else:
                return np.array([0.0] * self.nlay)

        try:
            self.nlay = len(thickn)
        except TypeError:
            self.nlay = 1

        self._thickn = _array(thickn)
        self._rho = _array(rho)
        self._vp = _array(vp)

        if vs is None:
            self._vpvs = _array(vpvs)
            self._vs = self._vp / self._vpvs
        else:
            self._vs = _array(vs)
            self._vpvs = self._vp / self._vs

        self._flag = np.array(
            [flag] * self.nlay if isinstance(flag, int) else list(flag)
        )
        self._ani = _array(ani)
        self._trend = _array(trend)
        self._plunge = _array(plunge)
        self._strike = _array(strike)
        self._dip = _array(dip)

        self.maxlay = maxlay

        self.properties = [
            "thickn",
            "rho",
            "vp",
            "vs",
            "vpvs",
            "flag",
            "ani",
            "trend",
            "plunge",
            "strike",
            "dip",
        ]
        self._properties = ["_" + prop for prop in self.properties]

        self._set_fattributes()
        self._set_layers()

    @property
    def thickn(self):
        return self._thickn

    @thickn.setter
    def thickn(self, value):
        self._thickn = value
        self._set_fattributes()
        self._set_layers()

    @property
    def rho(self):
        return self._rho

    @rho.setter
    def rho(self, value):
        self._rho = value
        self._set_fattributes()
        self._set_layers()

    @property
    def vp(self):
        return self._vp

    @vp.setter
    def vp(self, value):
        self._vp = value
        self._set_fattributes()
        self._set_layers()

    @property
    def vpvs(self):
        return self._vpvs

    @vpvs.setter
    def vpvs(self, value):
        self._vpvs = value
        self._set_fattributes()
        self._set_layers()
        self.update("vs")

    @property
    def vs(self):
        return self._vs

    @vs.setter
    def vs(self, value):
        self._vs = value
        self._set_fattributes()
        self._set_layers()

    @property
    def flag(self):
        return self._flag

    @flag.setter
    def flag(self, value):
        self._flag = value
        self._set_fattributes()
        self._set_layers()

    @property
    def ani(self):
        return self._ani

    @ani.setter
    def ani(self, value):
        self._ani = value
        self._set_fattributes()
        self._set_layers()

    @property
    def trend(self):
        return self._trend

    @trend.setter
    def trend(self, value):
        self._trend = value
        self._set_fattributes()
        self._set_layers()

    @property
    def plunge(self):
        return self._plunge

    @plunge.setter
    def plunge(self, value):
        self._plunge = value
        self._set_fattributes()
        self._set_layers()

    @property
    def strike(self):
        return self._strike

    @strike.setter
    def strike(self, value):
        self._strike = value
        self._set_fattributes()
        self._set_layers()

    @property
    def dip(self):
        return self._dip

    @dip.setter
    def dip(self, value):
        self._dip = value
        self._set_fattributes()
        self._set_layers()

    def __getitem__(self, ilay):
        """Get a layer or property of the model"""
        try:
            if isinstance(ilay[0], int) and isinstance(ilay[1], str):
                # model[0, "thickn"]
                return self.layers[ilay[0]][ilay[1]]
            if isinstance(ilay[1], int) and isinstance(ilay[0], str):
                # model["thickn", 0]
                return self.layers[ilay[1]][ilay[0]]
        except TypeError:
            try:
                # model[0]
                return self.layers[ilay]
            except TypeError:
                # model["thickn"]
                return np.array([lay[ilay] for lay in self.layers])

    def __setitem__(self, layatt, value):
        """Set a layer or property of the model"""

        lays = []
        atts = []
        vals = []
        try:
            # model[0] = {"plunge": 10} syntax
            # layatt is layer and value is {"att": value}
            atts = value.keys()
            for att in atts:
                vals.append(value[att])
            lays = [layatt] * len(vals)
        except AttributeError:
            # model[0, "plunge"] = 10 syntax
            # layatt[0] is layer, layatt[1] is attribute, value is value
            if isinstance(layatt[1], str) and isinstance(layatt[0], int):
                lays = [layatt[0]]
                atts = [layatt[1]]
            elif isinstance(layatt[0], str) and isinstance(layatt[1], int):
                lays = [layatt[1]]
                atts = [layatt[0]]
            else:
                msg = "Cannot set model layer attribute: {:}".format(layatt)
                raise ValueError(msg)
            vals = [value]

        for lay, att, val in zip(lays, atts, vals):
            if att not in self.properties:
                msg = f"Unknown attribute: '{att}'. Must be one of: "
                msg += ", ".join(self.properties)
                raise ValueError(msg)

            self.__dict__["_" + att][lay] = val

            if att == "ani" and val != 0:
                self.__dict__["_flag"][lay] = 0
            if att == "ani" and val == 0:
                self.__dict__["_flag"][lay] = 1

            if att == "vpvs":
                self.update(change="vs")
            else:
                self.update()

    def __len__(self):
        return self.nlay

    def __str__(self):
        buf = "#  thickn     rho      vp      vs  flag aniso   trend "
        buf += "plunge strike   dip\n"

        f = "{: 9.1f} {: 7.1f} {: 7.1f} {: 7.1f} {: 4.0f} {: 6.1f} {: 7.1f} "
        f += "{: 6.1f} {: 6.1f} {: 5.1f}\n"

        for th, vp, vs, r, fl, a, tr, p, s, d in zip(
            self._thickn,
            self._vp,
            self._vs,
            self._rho,
            self._flag,
            self._ani,
            self._trend,
            self._plunge,
            self._strike,
            self._dip,
        ):
            buf += f.format(th, r, vp, vs, fl, a, tr, p, s, d)

        return buf.strip("\n")

    def __add__(self, other):
        if not isinstance(other, Model):
            try:
                other = Model(**other)
            except TypeError:
                other = Model(*other)
            except Exception:
                msg = "Can only add Model, or valid dict or list to Model."
                raise TypeError(msg)

        third = deepcopy(self)

        for att in third._properties:
            third.__dict__[att] = np.append(self.__dict__[att], other.__dict__[att])

        third.nlay += other.nlay
        third._set_fattributes()
        third._set_layers()
        return third

    def __eq__(self, other):
        issame = [
            slay[att] == olay[att]
            for slay, olay in zip(self.layers, other.layers)
            for att in self.properties
        ]
        return all(issame)

    def _set_layers(self):
        self.layers = [
            {
                prop: self.__dict__[_prop][lay]
                for prop, _prop in zip(self.properties, self._properties)
            }
            for lay in range(self.nlay)
        ]

    def _set_fattributes(self):
        if self.nlay > self.maxlay:
            msg = f"The object is larger (nlay={self.nlay}) than the memory allocated "
            msg += f"at compile time (maxlay={self.maxlay}). "
            msg += (
                f"Increase maxlay in params.h and when constucting this Model object."
            )
            raise IndexError(msg)

        tail = np.zeros(self.maxlay - self.nlay)
        self.fthickn = np.asfortranarray(np.append(self._thickn, tail))
        self.frho = np.asfortranarray(np.append(self._rho, tail))
        self.fvp = np.asfortranarray(np.append(self._vp, tail))
        self.fvs = np.asfortranarray(np.append(self._vs, tail))
        self.fflag = np.asfortranarray(np.append(self._flag, tail))
        self.fani = np.asfortranarray(np.append(self._ani, tail))
        self.ftrend = np.asfortranarray(np.append(self._trend, tail) * np.pi / 180)
        self.fplunge = np.asfortranarray(np.append(self._plunge, tail) * np.pi / 180)
        self.fstrike = np.asfortranarray(np.append(self._strike, tail) * np.pi / 180)
        self.fdip = np.asfortranarray(np.append(self._dip, tail) * np.pi / 180)
        self.parameters = [
            self.fthickn,
            self.frho,
            self.fvp,
            self.fvs,
            self.fflag,
            self.fani,
            self.ftrend,
            self.fplunge,
            self.fstrike,
            self.fdip,
            self.nlay,
        ]

    def _v12str(self):
        """Legacy Raysum .mod file convention"""
        buf = "#  thickn     rho      vp      vs  flag p-aniso  s-aniso   trend "
        buf += "plunge strike   dip\n"

        f = "{: 9.1f} {: 7.1f} {: 7.1f} {: 7.1f} {: 4.0f} {: 8.1f} {: 8.1f} {: 7.1f} "
        f += "{: 6.1f} {: 6.1f} {: 5.1f}\n"

        for th, vp, vs, r, fl, a, tr, p, s, d in zip(
            self._thickn,
            self._vp,
            self._vs,
            self._rho,
            self._flag,
            self._ani,
            self._trend,
            self._plunge,
            self._strike,
            self._dip,
        ):
            buf += f.format(th, r, vp, vs, fl, a, a, tr, p, s, d)

        return buf.strip("\n")

    def update(self, change="vpvs"):
        """
        Update all attributes after one of them was changed.

        Parameters:
            change (str):
                Decide which of :py:attr:`vp`, :py:attr:`vs` or :py:attr:`vpvs` should
                depend on the other two:

                    * `'vpvs'`: Calculate :py:attr:`vpvs` from `vpvs = vp / vs`
                    * `'vp'`: Calculate :py:attr:`vp` from `vp = vs * vpvs`
                    * `'vs'`: Calculate :py:attr:`vs` from `vs = vp / vpvs`

        Returns:
            None
        """

        if change == "vp":
            self._vp = self._vs * self._vpvs
        elif change == "vs":
            self._vs = self._vp / self._vpvs
        elif change == "vpvs":
            self._vpvs = self._vp / self._vs
        else:
            msg = "Unknown value for keyword: " + change
            raise ValueError(msg)

        self._set_fattributes()
        self._set_layers()

    def copy(self):
        """Return a copy of the model"""
        return deepcopy(self)

    def change(self, command, verbose=True):
        """
        Change layer properties using a command string.

        Args:
            command (str):
                An arbitrary number of command substrings separated by ';'

            verbose (bool):
                Print changed parameters to screen

        Returns:
            list of 3*tuple:
               List of changes applied of the form:
               (attribute, layer, new value)

        Note:
            In the :data:`command` argument, each substring has the form:

            KEY LAYER SIGN VAL;

            where

            KEY [t|vp|vs|psp|pss|s|d|a|tr|pl] is the attribute to change

                * t: thickness (km)
                * vp: P wave velocity (km/s)
                * vs: S wave velocity (km/s)
                * psp: P to S wave velocity ratio with fixed vs (changing vp)
                * pss: P to S wave velocity ratio with fixed vp (changing vs)
                * s: strike (degree)
                * d: dip (degree)
                * a: anisotropy (%)
                * tr: trend of the anisotropy axis (degree)
                * pl: plunge ot the anisotropy axis (degree)

            LAYER (int) is the index of the layer

            SIGN [=|+|-] is to set / increase / decrease the attribute

            VAL (float) is the value to set / increase / decrease

        Hint:
            For example, ``Model.change('t0+10;psp0-0.2;d1+5;s1=45')`` does:

                1. Increase the thickness of the first layer by 10 km
                2. Decrease Vp/Vs of the of the first layer by 0.2, holding Vs
                   fixed
                3. Increase the dip of the second layer by 5 degree
                4. Set the strike of the second layer to 45 degree
        """

        ATT = {
            "t": "thickn",
            "thickn": "thickn",
            "vp": "vp",
            "vs": "vs",
            "psp": "vpvs",
            "pss": "vpvs",
            "s": "strike",
            "strike": "strike",
            "d": "dip",
            "dip": "dip",
            "a": "ani",
            "ani": "ani",
            "tr": "trend",
            "trend": "trend",
            "pl": "plunge",
            "plunge": "plunge",
        }
        _ATT = {key: "_" + val for key, val in zip(ATT.keys(), ATT.values())}

        changed = []
        for com in command.split(";"):
            com = com.strip()
            if not com:
                continue

            # split by sign
            for sign in "=+-":
                ans = com.split(sign, maxsplit=1)
                if len(ans) == 2:
                    break

            (attlay, val) = ans

            # Split attribute and layer
            for n, char in enumerate(attlay):
                if char in "0123456789":
                    break

            att = attlay[:n]
            lay = int(attlay[n:])
            val = float(val)

            # convert thicknes and velocities from kilometers
            if att in ["t", "vp", "vs"]:
                val *= 1000

            # Which velocity to fix
            change = "vpvs"
            if att == "pss":
                change = "vs"
            if att == "psp":
                change = "vp"

            attribute = ATT[att]
            _attribute = _ATT[att]

            # Apply
            if sign == "=":
                self.__dict__[_attribute][lay] = val
                sign = ""  # to print nicely below
            elif sign == "+":
                self.__dict__[_attribute][lay] += val
            elif sign == "-":
                self.__dict__[_attribute][lay] -= val

            # Set isotropy flag iff layer is isotropic
            self._flag[lay] = 1
            if self._ani[lay] != 0:
                self._flag[lay] = 0

            self.update(change=change)
            changed.append((attribute, lay, self.__dict__[_attribute][lay]))

            if verbose:
                msg = "Changed: {:}[{:d}] {:}= {:}".format(attribute, lay, sign, val)
                print(msg)

        return changed

    def split_layer(self, n):
        """
        Split layer `n` into two with half the thickness each, but otherwise
        identical parameters.

        Args:
            n : (int)
                Index of the layer to split
        """

        for att in self._properties:
            self.__dict__[att] = np.insert(self.__dict__[att], n, self.__dict__[att][n])

        self._thickn[n] /= 2
        self._thickn[n + 1] /= 2
        self.nlay += 1

        self.update()

    def remove_layer(self, n):
        """
        Remove layer `n`

        Args:
            n (int):
                Index of the layer to remove
        """

        for att in self._properties:
            self.__dict__[att] = np.delete(self.__dict__[att], n)

        self.nlay -= 1
        self.update()

    def average_layers(self, top, bottom):
        """
        Combine layers between top and bottom indices into one with summed
        thicknesses and averaged vp, vs, and rho.

        Args:
            top (int):
                Index before top-most layer to include in combination
            bottom (int):
                Index after bottom-most layer to include in combination

        Raises:
            IndexError: if bottom is less or equal to top
            ValueError: if any layer is anisotropic
            ValueError: if any layer has a difference in strike or dip
        """
        if bottom <= top:
            raise IndexError("bottom must be larger than top.")

        if not all(self._flag[top:bottom]):
            raise ValueError("Can only combine isotropic layers")

        if not all(self._dip[top:bottom][0] == self._dip[top:bottom]):
            raise ValueError("All layers must have the same dip")

        if not all(self._strike[top:bottom][0] == self._strike[top:bottom]):
            raise ValueError("All layers must have the same strike")

        thickn = sum(self._thickn[top:bottom])
        weights = self._thickn[top:bottom] / thickn

        layer = {
            "_thickn": thickn,
            "_vp": sum(self._vp[top:bottom] * weights),
            "_vs": sum(self._vs[top:bottom] * weights),
            "_rho": sum(self._rho[top:bottom] * weights),
        }

        for att in self._properties:
            try:
                self.__dict__[att][top] = layer[att]
            except KeyError:
                pass
            self.__dict__[att] = np.delete(self.__dict__[att], range(top + 1, bottom))

        self.nlay -= bottom - top - 1
        self.update()

    def save(self, fname="sample.mod", comment="", hint=False, version="prs"):
        """
        Alias for :class:`write()`
        """
        self.write(fname=fname, comment=comment, hint=hint, version=version)

    def write(self, fname="sample.mod", comment="", hint=False, version="prs"):
        """
        Write seismic velocity model to disk as Raysum ASCII model file

        Args:
            fname (str): Name of the output file (including extension)
            comment (str): String to write into file header
            hint (bool): Include usage comment to model file
            version ("prs" or "raysum"): Use PyRaysum or Raysum file format
        """

        if not comment.startswith("#"):
            comment = "# " + comment
        if not comment.endswith("\n"):
            comment += "\n"

        if not isinstance(fname, str):
            print("Warning: filename reverts to default 'sample.mod'")
            fname = "sample.mod"

        buf = "# Raysum velocity model created with PyRaysum\n"
        buf += "# on: {:}\n".format(datetime.now().isoformat(" ", "seconds"))

        if hint:
            buf += _modhint

        buf += comment

        if version == "prs":
            buf += self.__str__()
        elif version == "raysum":
            buf += self._v12str()
        else:
            msg = f"Unknown version: {version}"
            raise ValueError(msg)

        with open(fname, "w") as fil:
            fil.write(buf)

    def plot(self, zmax=75.0, show=True):
        """
        Plot model as a staircase, layers and labelled interfaces, and show it

        Args:
            zmax (float): Maximum depth of model to plot (km)

        """

        # Initialize new figure
        fig = plt.figure(figsize=(10, 5))

        # Add subplot for profile
        ax1 = fig.add_subplot(1, 4, 1)
        self.plot_profile(zmax=zmax, ax=ax1)

        # Add subplot for layers
        ax2 = fig.add_subplot(1, 4, 2)
        self.plot_layers(zmax=zmax, ax=ax2)

        ax3 = fig.add_subplot(1, 4, (3, 4))
        self.plot_interfaces(zmax=zmax, ax=ax3)

        # Tighten the plot and show it
        # TODO: tight layout does not work with colorbar.
        # Do this manually.
        # plt.tight_layout()
        if show:
            plt.show()
        else:
            plt.close(fig)
        return fig
        

    def plot_profile(self, zmax=75.0, ax=None):
        """
        Plot model as stair case and show it

        Args:
            zmax (float):
                Maximum depth of model to plot (km)
            ax (plt.axis):
                Axis handle for plotting. If ``None``, show the plot

        Returns:
            plt.axis:
                Axis handle for plotting
        """

        # Defaults to not show the plot
        show = False

        # Find depths of all interfaces in km
        thickn = self._thickn.copy()
        if thickn[-1] == 0.0:
            thickn[-1] = 50000.0
        depths = np.concatenate(([0.0], np.cumsum(thickn))) / 1000.0

        # Get corner coordinates of staircase representation of model
        depth = np.array(list(zip(depths[:-1], depths[1:]))).flatten()
        vs = np.array(list(zip(self._vs, self._vs))).flatten()
        vp = np.array(list(zip(self._vp, self._vp))).flatten()
        rho = np.array(list(zip(self._rho, self._rho))).flatten()
        ani = np.array(list(zip(self._ani, self._ani))).flatten()

        # Generate new plot if an Axis is not passed
        if ax is None:
            fig = plt.figure(figsize=(5, 5))
            ax = fig.add_subplot(111)
            show = True

        # Plot background model
        ax.plot(vs, depth, color="C0", label=r"Vs (m s$^{-1}$)")
        ax.plot(vp, depth, color="C1", label=r"Vp (m s$^{-1}$)")
        ax.plot(rho, depth, color="C2", label=r"Density (kg m$^{-3}$)")

        # If there is anisotropy, show variability
        if np.any([flag == 0 for flag in self._flag]):
            ax.plot(vs * (1.0 - ani / 100.0), depth, "--", color="C0")
            ax.plot(vs * (1.0 + ani / 100.0), depth, "--", color="C0")
            ax.plot(vp * (1.0 - ani / 100.0), depth, "--", color="C1")
            ax.plot(vp * (1.0 + ani / 100.0), depth, "--", color="C1")

        # Fix axes and add labels
        ax.legend(fontsize=8)
        ax.set_xlabel("Velocity or Density")
        ax.set_ylabel("Depth (km)")
        ax.set_ylim(0.0, zmax)
        ax.invert_yaxis()
        ax.grid(ls=":")

        if show:
            plt.show()

        return ax

    def plot_layers(self, zmax=75.0, ax=None):
        """
        Plot model as horizontal layers and show it

        Args:
            zmax (float):
                Maximum depth of model to plot (km)
            ax (plt.axis):
                Axis handle for plotting. If ``None``, show the plot

        Returns:
            plt.axis:
                Axis handle for plotting
        """

        # Defaults to not show the plot
        show = False

        # Find depths of all interfaces
        thickn = self._thickn.copy()
        if thickn[-1] == 0.0:
            thickn[-1] = 50000.0
        depths = np.concatenate(([0.0], np.cumsum(thickn))) / 1000.0

        # Generate new plot if an Axis is not passed
        if ax is None:
            fig = plt.figure(figsize=(2, 5))
            ax = fig.add_subplot(111)
            show = True
        else:
            fig = ax.get_figure()

        # Define color palette
        norm = plt.Normalize()
        colors = plt.cm.GnBu(norm(self._vs))

        # Cycle through layers
        for i in range(len(depths) - 1):

            # If anisotropic, add texture - still broken hatch
            if not self._flag[i] == 1:
                cax = ax.axhspan(depths[i], depths[i + 1], color=colors[i])
                cax.set_hatch("o")
            # Else isotropic
            else:
                cax = ax.axhspan(depths[i], depths[i + 1], color=colors[i])

        # Fix axes and labels
        ax.set_ylim(0.0, zmax)
        ax.set_xticks(())
        ax.invert_yaxis()
        pos = ax.get_position()

        ax2 = fig.add_axes([pos.x0, pos.y0 - 0.02, pos.width, 0.015])

        cmap = plt.cm.ScalarMappable(norm=norm, cmap="GnBu")
        plt.colorbar(cmap, cax=ax2, orientation="horizontal", label="$V_S$ (m/s)")

        if show:
            ax.set_ylabel("Depth (km)")
            plt.tight_layout()
            plt.show()

        return ax

    def plot_interfaces(self, zmax=75, info=["vp", "vs", "vpvs", "rho"], ax=None):
        """
        Plot model as labeled interfaces with possibly dipping layers

        Args:
            zmax (float):
                Maximum depth of model to plot (km)
            info (list of str):
                Which :class:`Model` attributes to list for each layer
            ax (plt.axis):
                Axis handle for plotting. If ``None``, show the plot

        Returns:
            plt.axis:
                Axis handle for plotting

        Raises:
            ValueError: If element in :const:`info` is not recognized

        """

        # Defaults to not show the plot
        show = False

        # Find depths of all interfaces
        depths = np.concatenate(([0.0], np.cumsum(self._thickn))) / 1000.0
        maxdep = depths[-1] + 50
        xs = np.array([-maxdep / 2, maxdep / 2])

        # Generate new plot if an Axis is not passed
        if ax is None:
            fig = plt.figure(figsize=(4, 4))
            ax = fig.add_subplot(111)
            show = True

        ax.scatter(0, -0.6, 60, marker="v", c="black")
        # Cycle through layers
        for i, depth in enumerate(depths[:-1]):
            dzdx = np.sin(self._dip[i] * np.pi / 180)
            zs = depth + xs * dzdx
            ax.plot(xs, zs, color="black")
            dipdir = (self._strike[i] + 90) % 360

            if i == 0 or self._strike[i] != self._strike[i - 1]:
                ax.text(
                    xs[-1], zs[-1], ">{:.0f}°".format(dipdir), ha="left", va="center"
                )

            if info:
                msg = "{: 2d}: ".format(i)

            for n, inf in enumerate(info):
                if inf == "vp":
                    msg += "$V_P={:.1f}$km/s".format(self._vp[i] / 1000)
                elif inf == "vs":
                    msg += "$V_S={:.1f}$km/s".format(self._vs[i] / 1000)
                elif inf == "vpvs":
                    msg += "$V_P/V_S={:.2f}$".format(self._vpvs[i])
                elif inf == "rho":
                    msg += "$\\rho={:.1f}$kg/m$^3$".format(self._rho[i] / 1000)
                else:
                    err = "Unknown argument to 'info': " + inf
                    raise ValueError(err)
                if len(info[n:]) > 1:
                    msg += ", "

            ax.text(
                0,
                depth + 1,
                msg,
                rotation=-self._dip[i],
                rotation_mode="anchor",
                ha="center",
                va="top",
            )

            if self._flag[i] == 0:
                aninfo = "-{:.0f}%-{:.0f}°".format(self._ani[i], self._trend[i])
                ax.text(
                    xs[-1],
                    zs[-1] + self._thickn[i] / 2000,
                    aninfo,
                    rotation=-self._plunge[i],
                    ha="left",
                    va="center",
                )

        # Fix axes and labels
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.xaxis.set_ticks_position("bottom")
        ax.yaxis.set_ticks_position("left")
        ax.set_ylim(0.0, zmax)
        ax.axis("equal")
        ax.invert_yaxis()

        if show:
            ax.set_ylabel("Depth (km)")
            plt.tight_layout()
            plt.show()

        return ax


class Geometry(object):
    """Ray geometry of seismic events and station configuration.

    One set of synthetic traces will be computed for each array element.

    .. note::
       This object holds the infromation of the .geom file in classic Raysum

    Parameters:
        baz (float or array_like):
          Ray back-azimuths (deg)
        slow (float or array_like):
          Ray slownesses (s/km)
        dn (float):
          North-offset of the seismic station (m)
        de (float):
          East-offset of the seismic station (m)
        maxtr (int):
          Maximum number of traces defined in params.h

    .. hint::
       If :attr:`baz` and :attr:`slow` do not have a different length, one ray will
       be computed for each *combination* of :attr:`baz` and :attr:`slow` (See examples)

    The following attributes are set upon initialization. `f` prefixes indicate
    attributes used for interaction with `fraysum.run_bare()` and
    `fraysum.run_full()`.

        ntr (int):
          Number of traces
        parameters (list):
          Convenience attribute that collects the `f`-attributes in the order expected
          by `fraysum.run_bare()` and `fraysum.run_full()`
        fbaz (np.ndarray):
          Ray back-azimuth (radians)
        fslow (np.ndarray):
          Ray slowness (m/s)
        fdn (np.ndarray):
          North-offset of the seismic station (m)
        fde (np.ndarray):
          East-offset of the seismic station (m)

    Example
    -------
        >>> from pyraysum import Geometry
        >>> geom = Geometry(60, 0.06)
        >>> print(geom)
        #Back-azimuth, Slowness, N-offset, E-offset
                60.00    0.0600      0.00      0.00
        >>> geom += [[90, 120], 0.04]
        >>> print(geom)
        #Back-azimuth, Slowness, N-offset, E-offset
                60.00    0.0600      0.00      0.00
                90.00    0.0400      0.00      0.00
               120.00    0.0400      0.00      0.00
        >>> geom = Geometry(range(0, 360, 90), [0.04, 0.08], 10, 35)
        >>> print(geom)
        #Back-azimuth, Slowness, N-offset, E-offset
                 0.00    0.0400     10.00     35.00
                90.00    0.0400     10.00     35.00
               180.00    0.0400     10.00     35.00
               270.00    0.0400     10.00     35.00
                 0.00    0.0800     10.00     35.00
                90.00    0.0800     10.00     35.00
               180.00    0.0800     10.00     35.00
               270.00    0.0800     10.00     35.00
    """

    def __init__(self, baz, slow, dn=0, de=0, maxtr=500):

        self.maxtr = maxtr
        self.properties = {"baz": 0, "slow": 1, "dn": 2, "de": 3}

        if type(baz) == int or type(baz) == float:
            baz = [baz]

        if type(slow) == int or type(slow) == float:
            slow = [slow]

        if len(baz) != len(slow):
            self._geom = [(bb, ss) for ss in slow for bb in baz]
        else:
            self._geom = [(bb, ss) for bb, ss in zip(baz, slow)]

        baz, slow = zip(*self._geom)
        self._baz = np.array(baz)
        self._slow = np.array(slow)
        self.ntr = len(self._baz)

        self._dn = np.full(self.ntr, dn)
        self._de = np.full(self.ntr, de)

        self.rays = [*zip(self._baz, self._slow, self._dn, self._de)]

        self._set_fattributes()

    @property
    def baz(self):
        return self._baz

    @property
    def slow(self):
        return self._slow

    @property
    def dn(self):
        return self._dn

    @property
    def de(self):
        return self._de

    def __getitem__(self, iray):
        if isinstance(iray, tuple):
            # geometry[0, "baz"]
            iprop = self.properties[iray[1]]
            return self.rays[iray[0]][iprop]
        else:
            try:
                # geometry["baz"]
                iprop = self.properties[iray]
                return [tup[iprop] for tup in self.rays]
            except KeyError:
                # geometry[0]
                return self.rays[iray]

    def __setitem__(self, iray, value):

        if not isinstance(value, tuple) or len(value) != 4:
            msg = "Can only set tuple(baz, slow, dn, de)"
            raise TypeError(msg)
        self.rays[iray] = value

        self._baz[iray] = value[0]
        self._slow[iray] = value[1]
        self._geom[iray] = (value[0], value[1])
        self._dn[iray] = value[2]
        self._de[iray] = value[3]

        self._set_fattributes()

    def __len__(self):
        return self.ntr

    def __add__(self, other):
        if not isinstance(other, Geometry):
            try:
                other = Geometry(**other)
            except TypeError:
                other = Geometry(*other)
            except Exception:
                msg = "Can only add Geometry, or valid dict or list to Geometry."
                raise TypeError(msg)

        third = deepcopy(self)

        third._baz = np.append(self._baz, other._baz)
        third._slow = np.append(self._slow, other._slow)
        third._dn = np.append(self._dn, other._dn)
        third._de = np.append(self._de, other._de)
        third._geom = np.append(self._geom, other._geom)
        third.ntr += other.ntr
        third.rays += other.rays
        third._set_fattributes()

        return third

    def __eq__(self, other):
        return (
            all(np.equal(self._baz, other._baz))
            and all(np.equal(self._slow, other._slow))
            and all(np.equal(self._dn, other._dn))
            and all(np.equal(self._de, other._de))
        )

    def __str__(self):
        out = "#Back-azimuth, Slowness, N-offset, E-offset\n"
        form = "{: 13.2f} {: 9.4f} {:9.2f} {:9.2f}\n"
        for bb, ss, xx, yy in zip(self._baz, self._slow, self._dn, self._de):
            out += form.format(bb, ss, xx, yy)
        return out.strip("\n")

    def _fstr(self):
        out = "#Back-azimuth, Slowness, N-offset, E-offset\n"
        form = "{: 13.2f} {: 9.1e} {:9.2f} {:9.2f}\n"
        for bb, ss, xx, yy in zip(self._baz, self._slow, self._dn, self._de):
            out += form.format(bb, ss * 1e-3, xx, yy)
        return out.strip("\n")

    def _set_fattributes(self):
        if self.ntr > self.maxtr:
            msg = f"The object is larger (ntr={self.ntr}) than the memory allocated "
            msg += f"at compile time (maxtr={self.maxtr}). "
            msg += (
                f"Increase maxtr in params.h and when constucting this Geometry object."
            )
            raise IndexError(msg)
        tail = np.zeros(self.maxtr - self.ntr)
        self.fbaz = np.asfortranarray(np.append(self._baz, tail) * np.pi / 180)
        self.fslow = np.asfortranarray(np.append(self._slow, tail) * 1e-3)
        self.fdn = np.asfortranarray(np.append(self._dn, tail))
        self.fde = np.asfortranarray(np.append(self._de, tail))
        self.parameters = [self.fbaz, self.fslow, self.fdn, self.fde, self.ntr]

    def copy(self):
        """Return a copy of the geometry"""
        return deepcopy(self)

    def save(self, fname="sample.geom"):
        """
        Alias for :meth:`write()`
        """
        self.write(fname=fname)

    def write(self, fname="sample.geom"):
        """
        Write ray geometry to disk as ascii file

        Parameters:
            fname (str):
                Name of file
        """

        with open(fname, "w") as f:
            f.write(self._fstr())

        print("Geometry written to: " + fname)

    def plot(self, show=True):
        """
        Plot ray geometry in polar coordinates

        Returns:
            plt.axis:
                Axis handle for plotting

        """

        fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
        ax.set_theta_zero_location("N")
        ax.set_theta_direction("clockwise")
        ax.scatter(self._baz * np.pi / 180, self._slow, s=200, color="black", zorder=10)
        for n, (b, s) in enumerate(zip(self._baz, self._slow)):
            t = "{:d}".format(n)
            b *= np.pi / 180
            ax.text(b, s, t, color="white", ha="center", va="center", zorder=12)
        ax.set_title("Ray Backazimuth and Slowness indices")

        if show:
            plt.show()

        return ax


class Control(object):
    """
    Run Control parameters for :meth:`prs.run()`.

    Parameters:
        wvtype (str):
            Wave type of incoming wavefield (:const:`'P'`, :const:`'SV'`, or :const:`'SH'`)
        mults (int):
            ID for calculating free surface multiples

                * 0: no multiple
                * 1: First interface multiples only
                * 2: all first-order multiples (once reflected from the surface)
                * 3: supply phases to be computed via :meth:`Control.set_phaselist`

        npts (int):
            Number of samples in time series
        dt (float):
            Sampling interval in seconds
        align (int):
            ID for time alignment of seismograms

                * 0 or "none": do not align
                * 1 or "P": align at `P`
                * 2 or "SV": align at `SV`
                * 3 or "SH": align at `SH`

        shift (float or None):
            Time shift in seconds (positive shift moves seismograms
            to greater lags). If ``None``, set to :attr:`dt`, to include direct wave when
            :const:`align` > 0.
        rot (int or str):
            ID for rotation:

                * 0 or "ZNE": up, north, east [left-handed]
                * 1 or "RTZ": radial, transverse, down [positive towards source]
                * 2 or "PVH": P-, SH-, SV-polarization [positive in ray direction]

        verbose (int or bool):
            Verbosity:

                * :const:`0` or :const:`False`: silent
                * :const:`1` or :const:`True`: verbose

        maxseg (int):
            Maximum number of segments per ray, as defined in param.h
        maxtr (int):
            Maximum number of traces, as defined in param.h
        maxph (int):
            Maximum number of phases per trace, as defined in param.h
        maxsamp (int):
            Maximum number of samples per trace, as defined in param.h

    Attributes:

        parameters (list):
          Convenience attribute that collects the attributes in the order expected
          by :meth:`fraysum.run_bare()` and :meth:`fraysum.run_full()`

    Warning:
        The option :const:`rot=0` returns the left-handed ZNE coordinates (Z positive
        up), conventionally used for seismometers. The internal Fortran rountines
        (everything imported through :mod:`fraysum`) return right-handed NEZ
        coordinates (Z positive down).
    """

    def __init__(
        self,
        verbose=0,
        wvtype="P",
        mults=0,
        npts=300,
        dt=0.025,
        align=1,
        shift=None,
        rot=1,
        maxseg=45,
        maxtr=500,
        maxph=40000,
        maxsamp=100000,
    ):

        self.parameters = [0] * 11
        self.verbose = verbose
        self.wvtype = wvtype
        self.mults = mults
        self.npts = npts
        self.dt = dt
        self.align = align
        self.rot = rot
        if shift is None:
            self.shift = self.dt
        else:
            self.shift = shift

        self._iphase = _iphase[self.wvtype]
        self._numph = np.int32(0)
        self._maxseg = maxseg
        self._maxph = maxph
        self._maxtr = maxtr
        self._maxsamp = maxsamp
        self._phaselist = np.asfortranarray(
            np.zeros((maxseg, 2, maxph), dtype=np.int32)
        )
        self._nseg = np.asfortranarray(np.zeros(maxph), dtype=np.int32)
        self._update()

    @property
    def verbose(self):
        return self._verbose

    @verbose.setter
    def verbose(self, value):
        if value not in [0, 1, "0", "1", True, False]:
            msg = "verbose must be 0 or 1, not: " + str(value)
            raise ValueError(msg)
        self._verbose = int(value)
        self.parameters[7] = self._verbose

    @property
    def wvtype(self):
        return self._wvtype

    @wvtype.setter
    def wvtype(self, value):
        if value not in ["P", "SV", "SH"]:
            msg = "wvtype must be 'P', 'SV', or 'SH', not: " + str(value)
            raise ValueError(msg)
        self._wvtype = value
        self._iphase = _iphase[self.wvtype]
        self.parameters[0] = self._iphase

    @property
    def mults(self):
        return self._mults

    @mults.setter
    def mults(self, value):
        if value not in [0, 1, 2, 3, "0", "1", "2", "3"]:
            msg = "mults must be 0, 1, 2, or 3, not: " + str(value)
            raise ValueError(msg)

        self._mults = int(value)
        self.parameters[1] = self.mults

    @property
    def npts(self):
        return self._npts

    @npts.setter
    def npts(self, value):
        self._npts = int(value)
        self.parameters[2] = self._npts

    @property
    def dt(self):
        return self._dt

    @dt.setter
    def dt(self, value):
        self._dt = float(value)
        self.parameters[3] = self._dt

    @property
    def align(self):
        return self._align

    @align.setter
    def align(self, value):
        alignids = [*_aligni] + [*_alignn]
        if value not in alignids:
            msg = (
                "align must be "
                + ", ".join([str(alignid) for alignid in alignids])
                + ". Not: "
                + str(value)
            )
            raise ValueError(msg)
        try:
            self._align = _aligni[value]
        except KeyError:
            self._align = int(value)
        self.parameters[4] = self._align

    @property
    def rot(self):
        return self._rot

    @rot.setter
    def rot(self, value):
        rotids = [*_roti] + [*_rotn]
        if value not in rotids:
            msg = (
                "rot must be: "
                + ", ".join([str(rotid) for rotid in rotids])
                + ". Not: "
                + str(value)
            )
            raise ValueError(msg)
        try:
            self._rot = _roti[value]
        except KeyError:
            self._rot = int(value)
        self.parameters[6] = self._rot

    @property
    def shift(self):
        return self._shift

    @shift.setter
    def shift(self, value):
        self._shift = float(value)
        self.parameters[5] = self._shift

    def __str__(self):
        out = "Run control parameters:\n\n"
        out += "Verbosity: "
        out += "{:}\n".format(self.verbose)
        out += "Incoming wave type: "
        out += "{:}\n".format(self.wvtype)
        out += "Multiples: "
        out += "{:}\n".format(self.mults)
        out += "Number of samples per trace: "
        out += "{:}\n".format(self.npts)
        out += "Sample rate (seconds): "
        out += "{:}\n".format(self.dt)
        out += "Alignment: "
        out += "{:}\n".format(_alignn[self.align])
        out += "Shift of traces (seconds): "
        out += "{:}\n".format(self.shift)
        out += "Rotation to output: "
        out += "{:}".format(_rotn[self.rot])
        return out

    def _prs_str(self):
        out = "# Verbosity\n"
        out += "{:}\n".format(int(self.verbose))
        out += "# Phase name\n"
        out += "{:}\n".format(self.wvtype)
        out += "# Multiples: 0 for none, 1 for Moho, 2 all first-order\n"
        out += "{:}\n".format(self.mults)
        out += "# Number of samples per trace\n"
        out += "{:}\n".format(self.npts)
        out += "# Sample rate (seconds)\n"
        out += "{:}\n".format(self.dt)
        out += "# Alignment: 0 is none, 1 aligns on P, 2 on SV, 3 on SH\n"
        out += "{:}\n".format(self.align)
        out += "# Shift of traces (seconds)\n"
        out += "{:}\n".format(self.shift)
        out += "# Rotation to output: 0 is ENZ, 1 is RTZ, 2 is PVH\n"
        out += "{:}\n".format(self.rot)
        return out

    def _v12_str(self):
        out = "# Multiples: 0 for none, 1 for Moho, 2 all first-order\n"
        out += "{:}\n".format(self.mults)
        out += "# Number of samples per trace\n"
        out += "{:}\n".format(self.npts)
        out += "# Sample rate (seconds)\n"
        out += "{:}\n".format(self.dt)
        out += "# Gaussian pulse width (seconds)\n"
        out += "1.\n"
        out += "# Alignment: 0 is none, 1 aligns on P, 2 on SV, 3 on SH\n"
        out += "{:}\n".format(self.align)
        out += "# Shift of traces (seconds)\n"
        out += "{:}\n".format(self.shift)
        out += "# Rotation to output: 0 is ENZ, 1 is RTZ, 2 is PVH\n"
        out += "{:}\n".format(self.rot)
        return out

    def set_phaselist(self, descriptors, equivalent=False, kinematic=False):
        """
        Explicitly set phases to be calculated.

        Parameters:

            descriptors (list of str):
                List of phase descriptors, where each descriptor is a string of
                consecutive PHASE SEGMENT pairs, where
                PHASE is:

                    * P upgoing P-wave
                    * S upgoing fast S-wave
                    * T upgoing slow S-wave
                    * p downgoing P-wave
                    * s downgoing fast S-wave
                    * t downgoing slow S-wave

                SEGMENT is the index of the model layer.

            equivalent (bool):
                Augment phaselist by equivalent phases (e.g., `'1P0P0s0P'`;
                `'1P0S0p0P'`)

            kinematic (bool):
                Limit equivalent phases to those with same polarity as input phases

        Raises:
            IndexError: If resulting phaselist is longer than :const:`maxph`.

        Caution:
            At every interface, S-wave energy split up into a slow and a fast S-wave,
            `S` and `T`, regardless of layer anisotropy. Make sure to include all
            contributions when setting a explicit phaselist. It is best to inspect the
            output of :meth:`Result.descriptors()` of a previous run with
            :attr:`mults=2` for the ray combinations relevant for your problem.

        Caution:
            Using :const:`equivalent=False` may produce incorrect amplitudes, because
            contributions of equivalent phases may be missing. At the same time,
            amplitudes are oftentimes correct within 10%. Using :const:`equivalent=True`
            may cause performance issues.

        Hint:
            In a two layer model (index `0` is the topmost subsurface layer, index `1`
            is the underlying half-space) the :const:`descriptors` compute the phases:

            * ``['1P0P']``: direct P-wave
            * ``['1P0S']``: P-to-S1 converted wave
            * ``['1P0T']``: P-to-S2 converted wave
            * ``['1P0P', '1P0S', '1P0T]``: direct P-wave and all P-to-S converted waves
            * ``['1P0P0s0S']``: P reflected s at the surface, reflected S at the top of the half-space

        See also:
            :meth:`Result.descriptors()` returns all phases computed in a previous run.

        See also:
            :meth:`~pyraysum.prs.equivalent_phases()` outputs equivalent phases explicitly
        """

        self.mults = 3
        phre = "[" + "".join(_phids.keys()) + "]"
        descriptors = list(descriptors)

        if equivalent:
            descriptors += equivalent_phases(descriptors, kinematic=kinematic)

        self._numph = np.int32(len(descriptors))
        for iph, dscr in enumerate(descriptors):
            self._nseg[iph] = sum([dscr.count(c) for c in "".join(_phids.keys())])
            for iseg, (mlay, mphs) in enumerate(
                zip(re.finditer(r"\d+", dscr), re.finditer(phre, dscr))
            ):
                layn = int(mlay.group(0))
                phid = _phids[mphs.group(0)]

                self._phaselist[iseg, 0, iph] = layn + 1  # Fortran indexing
                self._phaselist[iseg, 1, iph] = phid

        self._update()

    def null_phaselist(self, mults=2):
        """
        Do not use phaselist, but compute using :const:`mults` keyword.

        Args:
            mults (int):
                Set :attr:`Control.mults` to this value.
        """

        if mults > 2:
            msg = "Invalid value : " + str(mults)
            msg += ". To compute only specific phases, use set_phaselist()"
            raise ValueError(msg)

        self._phaselist = np.asfortranarray(
            np.zeros((self.maxseg, 2, self.maxph), dtype=np.int32)
        )
        self.mults = mults

    def _update(self):
        """
        Explicitly update :attr:`parameters`
        """

        self.parameters = [
            self._iphase,
            self.mults,
            self.npts,
            self.dt,
            self.align,
            self.shift,
            self.rot,
            self.verbose,
            self._nseg,
            self._numph,
            self._phaselist,
        ]

    def save(self, fname="raysum-param", version="prs"):
        """
        Alias for :class:`~pyraysum.prs.Control.write()`
        """

        self.write(fname=fname, version=version)

    def write(self, fname="raysum-param", version="prs"):
        """
        Write parameter file to disk

        Args:
            fname: (str)
                Name of file
            version: ("prs" or "raysum")
                Pyraysum or Raysum raysum compatible
        """

        if version == "prs":
            buf = self._prs_str()
        elif version == "raysum":
            buf = self._v12_str()
        else:
            msg = f"Unknown version: {version}"
            raise ValueError(msg)

        with open(fname, "w") as f:
            f.write(buf)


class Result(object):
    """
    Result of a PyRaysum wavefield simulation. 3-component synthetic seismograms
    are strored as a list of streams in the stream attribute. Includes methods
    to calculate receiver functions, filter and plot the results. See
    :func:`run()` for examples.

    Parameters:
        model (:class:`~pyraysum.prs.Model`):
            Subsurface velocity model
        geometry (:class:`~pyraysum.prs.Geometry`):
            Recording geometry
        rc (:class:`~pyraysum.prs.Control`):
            Run-control parameters
        streams (list):
            List of :class:`~obspy.core.Stream` objects.

    If created with :const:`mode='full'` in :meth:`run()`, the :attr:`stats` attribute
    of each :class:`~obspy.core.Trace` in each :class:`~obspy.core.Stream` in
    :attr:`Result.streams` holds the additional attributes:

        phase_times
            Arrival times of seismic phases
        phase_amplitudes
            Amplitudes of seismic phases
        phase_descriptors
            Descriptors of seismic phases. Index indicates layer through which phase
            propagates. (e.g. "1P0S", see also :meth:`descriptors()`)
        phase_names
            Short names of seismic phases (e.g. "PS")
        conversion_names
            Conversion names of seismic phases. Index indicates top of layer at which
            conversion occurs. (e.g. "P1S")

    """

    def __init__(self, model=None, geometry=None, rc=[], streams=[]):

        self.model = model
        self.geometry = geometry
        self.streams = streams
        self.rc = rc
        self.rfs = []

    def __str__(self):
        msg = "Result contains:\n"
        msg += "{:d} synthetic {:}-seismogram(s)\n".format(
            len(self.streams), _rotn[self.rc.rot]
        )
        msg += "{:d} synthetic receiver function(s)\n".format(len(self.rfs))
        if self.model:
            msg += "\nFrom subsurface model:\n"
            msg += self.model.__str__()
            msg += "\n"
        else:
            msg += "\nNo subsurface model stored\n"
        if self.geometry:
            baz = self.geometry._baz
            slow = self.geometry._slow
            msg += "\nBack-azimuth and slowness range is:\n"
            msg += "{:f} - {:f}° and {:f} - {:f}s/km".format(
                min(baz), max(baz), min(slow), max(slow)
            )
        else:
            msg += "\nNo geometry stored."

        return msg

    def __getitem__(self, iray):
        istreams = ["stream", "streams", "seis", "seismogram", "seismograms"]
        irfs = ["rf", "rfs"]

        if iray in istreams:
            return self.streams

        elif iray in irfs:
            return self.rfs

        elif isinstance(iray, int):
            stream = self.streams[iray]
            try:
                rf = self.rfs[iray]
            except IndexError:
                rf = Stream()
            return stream, rf

        else:
            msg = "Can only return integer ray index or key in: "
            msg += ", ".join(istreams, irfs)
            raise ValueError(msg)

    def __len__(self):
        return len(self.streams)

    def calculate_rfs(self):
        """
        Generate receiver functions from spectral division of displacement traces.
        Will be stored in :attr:`rflist`.

        Raises:
            ValueError: In case :class:`Control` parameters a unsuitable.
        """

        if self.rc.rot == 0:
            msg = "Receiver functions cannot be calculated in geographical "
            msg += "coordinates, i.e. rc.rot must not be 0"
            raise (ValueError(msg))
        elif self.rc.rot == 1:
            cmpts = ["R", "T", "Z"]
        elif self.rc.rot == 2:
            cmpts = ["V", "H", "P"]

        rflist = []

        # Cycle through list of displacement streams
        for stream in self.streams:

            # Calculate time axis
            npts = stream[0].stats.npts
            taxis = np.arange(-npts / 2.0, npts / 2.0) * stream[0].stats.delta

            # Extract 3-component traces from stream
            rtr = stream.select(component=cmpts[0])[0]
            ttr = stream.select(component=cmpts[1])[0]
            ztr = stream.select(component=cmpts[2])[0]

            # Deep copy and re-initialize data to 0.
            rfr = rtr.copy()
            rfr.data = np.zeros(len(rfr.data))
            rft = ttr.copy()
            rft.data = np.zeros(len(rft.data))

            # Fourier transform
            ft_rfr = fft(rtr.data)
            ft_rft = fft(ttr.data)
            ft_ztr = fft(ztr.data)

            # Spectral division to calculate receiver functions
            if self.rc.wvtype == "P":
                rfr.data = fftshift(np.real(ifft(np.divide(ft_rfr, ft_ztr))))
                rft.data = fftshift(np.real(ifft(np.divide(ft_rft, ft_ztr))))
            elif self.rc.wvtype == "SV":
                rfr.data = fftshift(np.real(ifft(np.divide(-ft_ztr, ft_rfr))))
            elif self.rc.wvtype == "SH":
                rft.data = fftshift(np.real(ifft(np.divide(-ft_ztr, ft_rft))))

            # Update stats
            rfr.stats.channel = "RF" + cmpts[0]
            rft.stats.channel = "RF" + cmpts[1]
            rfr.stats.taxis = taxis
            rft.stats.taxis = taxis

            # Store in Stream
            rfstream = Stream(traces=[rfr, rft])

            # Append to list
            rflist.append(rfstream)

        self.rfs = rflist

        return

    def descriptors(self):
        """
        Returns:
            list
                Unique list of all phase descriptors present

        Raises:
            AttributeError: If no descriptors are present

        Example
        -------
        >>> from pyraysum import Model, Control, Geometry, run
        >>> model = Model([30000., 0], [2800., 3300.], [6000., 8000.], [3600., 4500.])
        >>> geom = Geometry(0., 0.06)
        >>> rc = Control(mults=0)
        >>> seismogram = run(model, geom, rc)
        >>> seismogram.descriptors()
        ['1P0P', '1P0S']
        """

        phl = []
        for st in self.streams:
            for tr in st:
                try:
                    phl.extend(tr.stats.phase_descriptors)
                except AttributeError:
                    msg = "No phase descriptors present. Did you run in 'full' mode?"
                    raise AttributeError(msg)

        return sorted(list(set(phl)))

    def write(self, fname="sample.tr"):
        """
        Write streams to file, using legacy Raysum trace format

        Parameters:
            fname (str):
                Name of the output file (including extension)
        """

        buf = "#{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}\n".format(
            "#traces", "#samples", "dt (s)", "align", "shift"
        )
        buf += " {:11d}{:11d}{:11.3f}{:11d}{:11.3f}\n".format(
            len(self.streams), self.rc.npts, self.rc.dt, self.rc.align, self.rc.shift
        )
        for itr, stream in enumerate(self.streams):
            arr = np.vstack((stream[0].data, stream[1].data, stream[2].data)).T
            buf += "#-------------------\n"
            buf += "# Trace number {:5d}\n".format(itr + 1)  # Fortran indexing
            buf += "#-------------------\n"
            buf += (
                np.array2string(
                    arr,
                    threshold=3 * arr.shape[0] + 1,
                    formatter={"float": lambda i: "{:15.7e}".format(i)},
                )
                .replace("[", " ")
                .replace("]", " ")
            )
            buf += "\n"

        with open(fname, "w") as fid:
            fid.write(buf)

    def plot(self, typ, **kwargs):
        """
        Plot the displacement seismograms and/or receiver functions stored in
        :class:`~pyraysum.prs.Result` streams.

        Parameters:
            typ (str):
                Type of plot to show. Options are :const:`'streams'`, :const:`'rfs'`, or
                :const:`'all'` for the displacement seismograms, receiver functions, or
                both
            **kwargs:
                Are passed to the underlying :meth:`~pyraysum.plot.stream_wiggles` or
                :meth:`~pyraysum.plot.rf_wiggles`
        Hint:
            If :class:`Result` contains only one event (i.e., single :attr:`baz` and
            :attr:`slow` values), it is more appropriate to plot the waveforms using the
            default :meth:`Stream.plot()` method on :attr:`Result.streams[0]`.


        """

        if typ == "streams":
            self.plot_streams(**kwargs)
        elif typ == "rfs":
            self.plot_rfs(**kwargs)
        elif typ == "all":
            self.plot_streams(**kwargs)
            self.plot_rfs(**kwargs)
        else:
            msg = "'typ' has to be either 'streams', 'rfs' or 'all'"
            raise (ValueError(msg))

    def filter(self, typ, ftype, **kwargs):
        """

        Filters the displacement seismograms and/or receiver functions stored
        in :class:`~pyraysum.prs.Result` streams.

        Parameters:
            typ (str):
                Type of plot to show. Options are :const:`'streams'`,
                :const:`'rfs'`, or :const:`'all'` for the displacement seismograms,
                receiver functions, or both
            ftype (str):
                Type of filter to use in
                `obspy.Trace.filter <https://tinyurl.com/45dkyvwy>`_
            **kwargs:
                Keyword arguments passed to
                `obspy.Trace.filter <https://tinyurl.com/45dkyvwy>`_

        """
        if typ == "streams":
            self.filter_streams(ftype, **kwargs)
        elif typ == "rfs":
            self.filter_rfs(ftype, **kwargs)
        elif typ == "all":
            self.filter_streams(ftype, **kwargs)
            try:
                self.filter_rfs(ftype, **kwargs)
            except:
                print("Cannot filter 'rfs'. Continuing.")
        else:
            msg = "'typ' has to be either 'streams', 'rfs' or 'all'"
            raise (TypeError(msg))

    def plot_streams(self, scale=1.0e3, tmin=-5.0, tmax=20.0):
        plot.stream_wiggles(self, scale=scale, tmin=tmin, tmax=tmax)

    def plot_rfs(self, scale=1.0e3, tmin=-5.0, tmax=20.0):
        plot.rf_wiggles(self, scale=scale, tmin=tmin, tmax=tmax)

    def filter_streams(self, ftype, **kwargs):
        [stream.filter(ftype, **kwargs) for stream in self.streams]

    def filter_rfs(self, ftype, **kwargs):
        [rf.filter(ftype, **kwargs) for rf in self.rfs]


def run(model, geometry, rc, mode="full", rf=False):
    """
    Run a wave-field simulation. This function calls the compiled :mod:`fraysum`
    binaries.

    Parameters:
        model (:class:`~pyraysum.prs.Model`):
            Subsurface velocity model
        geometry (:class:`~pyraysum.prs.Geometry`):
            Recording geometry
        rc (:class:`~pyraysum.prs.Control`):
            Run-control parameters
        mode (str):
            * :const:`'full'`: Compute seismograms, phase arrivals and descriptors (slower)
            * :const:`'bare'`: Only compute seismograms (faster)
        rf (bool):
            Whether or not to calculate receiver functions

    Returns:
        :class:`~pyraysum.prs.Result`:
            Synthetic seimograms

    Hint:
        Best performance is achived by calling :func:`fraysum.run_bare`

    Example
    -------
    >>> from pyraysum import prs, Model, Geometry, Control
    >>> # Define two-layer model with isotropic crust over isotropic half-space
    >>> model = Model([30000., 0], [2800., 3300.], [6000., 8000.], [3600., 4500.])
    >>> geom = Geometry(0., 0.06) # baz = 0 deg; slow = 0.06 s/km
    >>> rc = Control(npts=1500, dt=0.025)
    >>> result = prs.run(model, geom, rc, rf=True)
    >>> type(result.streams[0])
    <class 'obspy.core.stream.Stream'>
    >>> type(result.rfs[0])
    <class 'obspy.core.stream.Stream'>
    >>> seis, rf = result[0]
    >>> print(seis)
    3 Trace(s) in Stream:
    ...SYR | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:37.475000Z | 40.0 Hz, 1500 samples
    ...SYT | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:37.475000Z | 40.0 Hz, 1500 samples
    ...SYZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:37.475000Z | 40.0 Hz, 1500 samples
    >>> print(rf)
    2 Trace(s) in Stream:
    ...RFR | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:37.475000Z | 40.0 Hz, 1500 samples
    ...RFT | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:37.475000Z | 40.0 Hz, 1500 samples
    """

    if rf and rc.rot == 0:
        msg = "Receiver functions cannot be calculated in ZNE coordinates, "
        msg += "i.e. in rc, 'rot' must not be '0'"
        raise ValueError(msg)

    if mode == "full":
        traces, traveltimes, amplitudes, phaselist = fraysum.run_full(
            *model.parameters, *geometry.parameters, *rc.parameters
        )

        if rc.mults == 3:
            #  phaselist has not been populated by fraysum
            phaselist = rc._phaselist

        arrivals = read_arrivals(traveltimes, amplitudes, phaselist, geometry)

    elif mode == "bare":
        traces = fraysum.run_bare(
            *model.parameters, *geometry.parameters, *rc.parameters
        )

        arrivals = None

    else:
        raise ValueError("Unknown mode: " + str(mode))

    # Read all traces and store them into a list of :class:`~obspy.core.Stream`
    streams = read_traces(traces, rc, geometry, arrivals=arrivals)

    # Store everything into Result object
    seismogram = Result(model=model, geometry=geometry, rc=rc, streams=streams)

    if rf:
        seismogram.calculate_rfs()

    return seismogram


def read_model(modfile, encoding=None, version="prs"):
    """
    Reads model parameters from file

    Parameters:
        version ("prs" or "raysum"):
            File has PyRaysum or Raysum format

    Returns:
        :class:`~pyraysum.prs.Model`:
            Seismic velocity model
    """

    vals = np.genfromtxt(modfile, encoding=encoding).T
    if version == "raysum":
        # ignore s-anisotropy
        vals = np.vstack((vals[:6], vals[7:]))

    return Model(*vals)


def read_geometry(geomfile, encoding=None):
    """
    Reads geometry parameters from file

    Returns:
        :class:`~pyraysum.prs.Geometry`:
            Ray geometry
    """

    vals = np.genfromtxt(geomfile, dtype=None, encoding=encoding)
    vals[:, 1] *= 1e3

    try:
        geom = Geometry(*zip(*vals))
    except TypeError:
        # *values contains single line, which is not iterable
        geom = Geometry([vals[0]], [vals[1]], [vals[2]], [vals[3]])
    return geom


def read_control(paramfile, version="prs"):
    """
    Read Raysum run control parameters from file.

    Parameters:
        version: ("prs" or "raysum")
            File has PyRaysum or Raysum format

    Returns:
        :class:`~pyraysum.prs.Control`:
            Run control parameters
    """

    with open(paramfile, "r") as f:
        lines = f.readlines()

    buf = []
    for line in lines:
        line = line.strip()
        if line.startswith("#"):
            continue
        buf.append(line)

    if version == "prs":
        return Control(*buf)
    elif version == "raysum":
        rc = Control(
            mults=buf[0], npts=buf[1], dt=buf[2], align=buf[4], shift=buf[5], rot=buf[6]
        )
        return rc


def equivalent_phases(descriptors, kinematic=False):
    """
    Return descriptors of equivalent phases, i.e. those arriving at the same time as
    the input phases

    Parameters:
        descriptors (list of strings):
            Phase descriptors as in :meth:`Control.set_phaselist`

        kinematic (boolean):
            If True, restrict to kinematically equivalent phases, i.e. those that have
            the same polarization as input phases

    Returns:
        list of strings:
            List of unique descriptors of equivalent phases

    Example
    -------
    >>> from pyraysum import prs
    >>> prs.equivalent_phases(['1P0P0p0S'])
    ['1P0P0s0P', '1P0S0p0P']
    """

    ndscrs = []

    for dscr in descriptors:

        lays = np.array([match.group(0) for match in re.finditer(r"\d+", dscr)])

        starts = np.array(
            [match.start() for match in re.finditer(r"\d+", dscr)] + [len(dscr)]
        )

        nsegs = {lay: len(lays[lays == lay]) for lay in set(lays)}

        for lay in lays:

            if nsegs[lay] <= 1:
                continue

            for iseg, i0 in enumerate(starts[:-1]):

                # phase is last char before next segment
                iph = starts[iseg + 1] - 1

                # the actural phase name
                ipha = dscr[iph]

                # layer is the chars before that
                ilay = dscr[starts[iseg] : iph]

                for jseg, j0 in enumerate(starts[:-1]):

                    # See above comments
                    jph = starts[jseg + 1] - 1
                    jpha = dscr[jph]
                    jlay = dscr[starts[jseg] : jph]

                    # Only proceed for phases in same layer
                    # but not for the same ray segment
                    if ilay != lay or jlay != lay or iseg == jseg:
                        continue

                    ndscr = np.array(list(dscr))

                    # Exchange wavetypes
                    # preserve propagation direction
                    if ipha.isupper():
                        kpha = jpha.upper()
                    else:
                        kpha = jpha.lower()

                    if jpha.isupper():
                        lpha = ipha.upper()
                    else:
                        lpha = ipha.lower()

                    ndscr[iph] = kpha
                    ndscr[jph] = lpha

                    ndscr = "".join(ndscr)

                    if ndscr in descriptors:
                        # Bail out solution already in input
                        continue

                    if kinematic:
                        # Bail out dissimilar last phase
                        lpin = dscr[-1]
                        lpout = ndscr[-1]
                        if lpin == "P" and (lpout == "S" or lpout == "T"):
                            continue
                        if (lpin == "S" or lpin == "T") and lpout == "P":
                            continue

                    ndscrs.append(ndscr)

    return sorted(list(set(ndscrs)))
