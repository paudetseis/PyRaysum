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

import re
from datetime import datetime
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from obspy import Trace, Stream, UTCDateTime
from numpy.fft import fft, ifft, fftshift
from pyraysum import plot
import fraysum

_align = {0: "None", 1: "P", 2: "SV", 3: "SH"}
_rot = {0: "ZNE", 1: "RTZ", 2: "PVH"}
_iphase = {"P": 1, "SV": 2, "SH": 3}
_phnames = {1: "P", 2: "S", 3: "T", 4: "p", 5: "s", 6: "t"}
_phids = {_phnames[k]: k for k in _phnames}  # inverse dictionary


class Model(object):
    """
    Model of the subsurface seismic velocity structure

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

    The following attributes are set upon initialization and when executing
    :meth:`update()`. `f` prefixes indicate attributes used for interaction with
    `fraysum.run_bare()` and `fraysum.run_full()`.

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
        def _get_val(v):
            if v is not None:
                return np.array(
                    [v] * self.nlay if isinstance(v, (int, float)) else v, dtype=float
                )
            else:
                return np.array([0.0] * self.nlay)

        self.nlay = len(thickn)
        self.thickn = np.array(thickn)
        self.rho = _get_val(rho)
        self.vp = _get_val(vp)

        if vs is None:
            self.vpvs = _get_val(vpvs)
            self.vs = self.vp / self.vpvs
        else:
            self.vs = np.array(vs)
            self.vpvs = self.vp / self.vs

        self.flag = np.array(
            [flag] * self.nlay if isinstance(flag, int) else list(flag)
        )
        self.ani = _get_val(ani)
        self.trend = _get_val(trend)
        self.plunge = _get_val(plunge)
        self.strike = _get_val(strike)
        self.dip = _get_val(dip)

        self.maxlay = maxlay

        self._set_fattributes()

        self._useratts = [
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

    def __len__(self):
        return self.nlay

    def __str__(self):
        buf = "# thickn     rho      vp      vs  flag aniso   trend "
        buf += "plunge strike   dip\n"

        f = "{: 8.1f} {: 7.1f} {: 7.1f} {: 7.1f} {: 4.0f} {: 6.1f} {: 7.1f} "
        f += "{: 6.1f} {: 6.1f} {: 5.1f}\n"

        for th, vp, vs, r, fl, a, tr, p, s, d in zip(
            self.thickn,
            self.vp,
            self.vs,
            self.rho,
            self.flag,
            self.ani,
            self.trend,
            self.plunge,
            self.strike,
            self.dip,
        ):
            buf += f.format(th, r, vp, vs, fl, a, tr, p, s, d)

        return buf

    def _set_fattributes(self):
        tail = np.zeros(self.maxlay - self.nlay)
        self.fthickn = np.asfortranarray(np.append(self.thickn, tail))
        self.frho = np.asfortranarray(np.append(self.rho, tail))
        self.fvp = np.asfortranarray(np.append(self.vp, tail))
        self.fvs = np.asfortranarray(np.append(self.vs, tail))
        self.fflag = np.asfortranarray(np.append(self.flag, tail))
        self.fani = np.asfortranarray(np.append(self.ani, tail))
        self.ftrend = np.asfortranarray(np.append(self.trend, tail) * np.pi / 180)
        self.fplunge = np.asfortranarray(np.append(self.plunge, tail) * np.pi / 180)
        self.fstrike = np.asfortranarray(np.append(self.strike, tail) * np.pi / 180)
        self.fdip = np.asfortranarray(np.append(self.dip, tail) * np.pi / 180)
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

    def update(self, change="vpvs"):
        """
        Update all attributes after one of them was changed by the user.

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
            self.vp = self.vs * self.vpvs
        elif change == "vs":
            self.vs = self.vp / self.vpvs
        elif change == "vpvs":
            self.vpvs = self.vp / self.vs
        else:
            msg = "Unknown value for keyword: " + change
            raise ValueError(msg)

        self._set_fattributes()

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
            In the ``command`` argument, each substring has the form:

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
            For example, :py:meth:`Model.change('t0+10;psp0-0.2;d1+5;s1=45')` does:

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

            # Apply
            if sign == "=":
                self.__dict__[attribute][lay] = val
                sign = ""  # to print nicely below
            elif sign == "+":
                self.__dict__[attribute][lay] += val
            elif sign == "-":
                self.__dict__[attribute][lay] -= val

            # Set isotropy flag iff layer is isotropic
            self.flag[lay] = 1
            if self.ani[lay] != 0:
                self.flag[lay] = 0

            self.update(change=change)
            changed.append((attribute, lay, self.__dict__[attribute][lay]))

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

        for att in self._useratts:
            self.__dict__[att] = np.insert(self.__dict__[att], n, self.__dict__[att][n])

        self.thickn[n] /= 2
        self.thickn[n + 1] /= 2
        self.nlay += 1

        self.update()

    def remove_layer(self, n):
        """
        Remove layer `n`

        Args:
            n (int):
                Index of the layer to remove
        """

        for att in self._useratts:
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

        if not all(self.flag[top:bottom]):
            raise ValueError("Can only combine isotropic layers")

        if not all(self.dip[top:bottom][0] == self.dip[top:bottom]):
            raise ValueError("All layers must have the same dip")

        if not all(self.strike[top:bottom][0] == self.strike[top:bottom]):
            raise ValueError("All layers must have the same strike")

        thickn = sum(self.thickn[top:bottom])
        weights = self.thickn[top:bottom] / thickn

        layer = {
            "thickn": thickn,
            "vp": sum(self.vp[top:bottom] * weights),
            "vs": sum(self.vs[top:bottom] * weights),
            "rho": sum(self.rho[top:bottom] * weights),
        }

        for att in self._useratts:
            try:
                self.__dict__[att][top] = layer[att]
            except KeyError:
                pass
            self.__dict__[att] = np.delete(self.__dict__[att], range(top + 1, bottom))

        self.nlay -= bottom - top - 1
        self.update()

    def save(self, fname="sample.mod", comment=""):
        """
        Alias for :class:`write()`
        """
        self.write(fname=fname, comment=comment)

    def write(self, fname="sample.mod", comment=""):
        """
        Write seismic velocity model to disk as Raysum ASCII model file

        Args:
            fname (str): Name of the output file (including extension)
            comment (str): String to write into file header

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
        buf += comment
        buf += self.__str__()

        with open(fname, "w") as fil:
            fil.write(buf)

    def plot(self, zmax=75.0):
        """
        Plot model as stair case, layers and labeled interfaces, and show it

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
        plt.tight_layout()
        plt.show()

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
        thickn = self.thickn.copy()
        if thickn[-1] == 0.0:
            thickn[-1] = 50000.0
        depths = np.concatenate(([0.0], np.cumsum(thickn))) / 1000.0

        # Get corner coordinates of staircase representation of model
        depth = np.array(list(zip(depths[:-1], depths[1:]))).flatten()
        vs = np.array(list(zip(self.vs, self.vs))).flatten()
        vp = np.array(list(zip(self.vp, self.vp))).flatten()
        rho = np.array(list(zip(self.rho, self.rho))).flatten()
        ani = np.array(list(zip(self.ani, self.ani))).flatten()

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
        if np.any([flag == 0 for flag in self.flag]):
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
        thickn = self.thickn.copy()
        if thickn[-1] == 0.0:
            thickn[-1] = 50000.0
        depths = np.concatenate(([0.0], np.cumsum(thickn))) / 1000.0

        # Generate new plot if an Axis is not passed
        if ax is None:
            fig = plt.figure(figsize=(2, 5))
            ax = fig.add_subplot(111)
            show = True

        # Define color palette
        norm = plt.Normalize()
        colors = plt.cm.GnBu(norm(self.vs))

        # Cycle through layers
        for i in range(len(depths) - 1):

            # If anisotropic, add texture - still broken hatch
            if not self.flag[i] == 1:
                cax = ax.axhspan(depths[i], depths[i + 1], color=colors[i])
                cax.set_hatch("o")
            # Else isotropic
            else:
                cax = ax.axhspan(depths[i], depths[i + 1], color=colors[i])

        # Fix axes and labels
        ax.set_ylim(0.0, zmax)
        ax.set_xticks(())
        ax.invert_yaxis()

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
        depths = np.concatenate(([0.0], np.cumsum(self.thickn))) / 1000.0
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
            dzdx = np.sin(self.dip[i] * np.pi / 180)
            zs = depth + xs * dzdx
            ax.plot(xs, zs, color="black")
            dipdir = (self.strike[i] + 90) % 360

            if i == 0 or self.strike[i] != self.strike[i - 1]:
                ax.text(
                    xs[-1], zs[-1], ">{:.0f}°".format(dipdir), ha="left", va="center"
                )

            if info:
                msg = "{: 2d}: ".format(i)

            for n, inf in enumerate(info):
                if inf == "vp":
                    msg += "$V_P={:.1f}$km/s".format(self.vp[i] / 1000)
                elif inf == "vs":
                    msg += "$V_S={:.1f}$km/s".format(self.vs[i] / 1000)
                elif inf == "vpvs":
                    msg += "$V_P/V_S={:.2f}$".format(self.vpvs[i])
                elif inf == "rho":
                    msg += "$\\rho={:.1f}$kg/m$^3$".format(self.rho[i] / 1000)
                else:
                    err = "Unknown argument to 'info': " + inf
                    raise ValueError(err)
                if len(info[n:]) > 1:
                    msg += ", "

            ax.text(
                0,
                depth + 1,
                msg,
                rotation=-self.dip[i],
                rotation_mode="anchor",
                ha="center",
                va="top",
            )

            if self.flag[i] == 0:
                aninfo = "-{:.0f}%-{:.0f}°".format(self.ani[i], self.trend[i])
                ax.text(
                    xs[-1],
                    zs[-1] + self.thickn[i] / 2000,
                    aninfo,
                    rotation=-self.plunge[i],
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
    """
    Recording geometry of rays and events at the seismic station. 
    One set of synthetic traces will be computed for each array element.

    Parameters:
        baz (float or array_like):
          Ray back-azimuths (deg)
        slow (float or array_like):
          Ray slownesses (s/km)
        dn (array_like):
          North-offset of the seismic station (m)
        de (array_like):
          East-offset of the seismic station (m)
        maxtr (int):
          Maximum number of traces defined in params.h

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
    """

    def __init__(self, baz, slow, dn=[0], de=[0], maxtr=500):

        if type(baz) == int or type(baz) == float:
            baz = [baz]

        if type(slow) == int or type(slow) == float:
            slow = [slow]

        if len(baz) != len(slow):
            self.geom = [(bb, ss) for ss in slow for bb in baz]
        else:
            self.geom = [(bb, ss) for bb, ss in zip(baz, slow)]

        baz, slow = zip(*self.geom)
        self.baz = np.array(baz)
        self.slow = np.array(slow)

        self.ntr = len(self.baz)

        self.dn = np.array(dn)
        self.de = np.array(de)

        if len(self.dn) != self.ntr:
            self.dn = np.full(self.ntr, self.dn[0])

        if len(self.de) != self.ntr:
            self.de = np.full(self.ntr, self.de[0])

        tail = np.zeros(maxtr - self.ntr)
        self.fbaz = np.asfortranarray(np.append(self.baz, tail) * np.pi / 180)
        self.fslow = np.asfortranarray(np.append(self.slow, tail) * 1e-3)
        self.fdn = np.asfortranarray(np.append(self.dn, tail))
        self.fde = np.asfortranarray(np.append(self.de, tail))
        self.parameters = [self.fbaz, self.fslow, self.fdn, self.fde, self.ntr]

    def __len__(self):
        return self.ntr

    def __str__(self):
        out = "#Back-azimuth, Slowness, N-offset, E-offset"
        form = "{: 7.2f} {: 8.4f} {:7.2f} {:7.2f}\n"
        for bb, ss, xx, yy in zip(self.baz, self.slow, self.dn, self.de):
            out += form.format(bb, ss, xx, yy)
        return out

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
            f.write(self.__str__())

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
        ax.scatter(self.baz * np.pi / 180, self.slow, s=200, color="black", zorder=10)
        for n, (b, s) in enumerate(zip(self.baz, self.slow)):
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
                * 3: supply phases to be computed via meth:`Control.set_phaselist()`

        npts (int):
            Number of samples in time series
        dt (float):
            Sampling interval in seconds
        align (int):
            ID for time alignment of seismograms

                * 0: do not align
                * 1: align at `P`
                * 2: align at `SV`
                * 3: align at `SH`

        shift (float or None):
            Time shift in seconds (positive shift moves seismograms
            to greater lags). If ``None``, set to :attr:`dt`, to include direct wave when
            :const:`align` > 0.
        rot (int):
            ID for rotation:

                * 0: ZNE (up, north, east [left-handed]);
                * 1: RTZ (radial, transverse, down [positive towards source])
                * 2: PVH (P-, SH-, SV-polarization [positive in ray direction])

        verbose (int):
            Verbosity. :const:`0` - silent; :const:`1` - verbose
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

        if verbose not in [0, 1, "0", "1"]:
            msg = "verbose must be 0 or 1, not: " + str(verbose)
            raise ValueError(msg)

        if wvtype not in ["P", "SV", "SH"]:
            msg = "wvtype must be 'P', 'SV', or 'SH', not: " + str(wvtype)
            raise ValueError(msg)

        if mults not in [0, 1, 2, 3, "0", "1", "2", "3"]:
            msg = "mults must be 0, 1, 2, or 3, not: " + str(mults)
            raise ValueError(msg)

        if align not in [0, 1, 2, 3, "0", "1", "2", "3"]:
            msg = "align must be 0, 1, 2, or 3, not: " + str(align)
            raise ValueError(msg)

        if rot not in [0, 1, 2, "0", "1", "2"]:
            msg = "rot must be 0, 1, or 2, not: " + str(rot)
            raise ValueError(msg)

        self.verbose = int(verbose)
        self.wvtype = wvtype
        self.mults = int(mults)
        self.npts = int(npts)
        self.dt = float(dt)
        self.align = int(align)
        self.rot = int(rot)
        if shift is None:
            self.shift = self.dt
        else:
            self.shift = float(shift)

        self.iphase = _iphase[self.wvtype]

        self._numph = np.int32(0)
        self._maxseg = maxseg
        self._maxph = maxph
        self._maxtr = maxtr
        self._maxsamp = maxsamp
        self._phaselist = np.asfortranarray(
            np.zeros((maxseg, 2, maxph), dtype=np.int32)
        )
        self._nseg = np.asfortranarray(np.zeros(maxph), dtype=np.int32)
        self.update()

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
        out += "{:}\n".format(_align[self.align])
        out += "Shift of traces (seconds): "
        out += "{:}\n".format(self.shift)
        out += "Rotation to output: "
        out += "{:}".format(_rot[self.rot])
        return out

    def __repr__(self):
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
                Augment phaselist by equivalent phases (e.g., '1P0P0s0P'; '1P0S0p0P')

            kinematic (bool):
                Limit equivalent phases to those with same polarity as input phases

        Raises:
            IndexError: If resulting phaselist is longer than :const:`maxph`.

        Hint:
            In a two layer model (index `0` is the topmost subsurface layer, index `1` 
            is the underlying half-space) the :const:`descriptors` compute the phases:

            * ``['1P0P']``: direct P-wave
            * ``['1P0S']``: P-to-S converted wave
            * ``['1P0P', '1P0S']``: direct P-wave and P-to-S converted wave
            * ``['1P0P0s0S']``: P reflected s at the surface, reflected S at the top of the half-space

        Caution:
            Using :const:`equivalent=False` may produce incorrect amplitudes, because
            contributions of equivalent phases may be missing.

        See also:
            To output equivalent phases, see :meth:`~pyraysum.prs.equivalent_phases()`
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

        self.update()

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

    def update(self):
        """
        Explicitly update :attr:`parameters`
        """

        self.parameters = [
            self.iphase,
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

    def save(self, fname="raysum-param"):
        """
        Alias for :class:`~pyraysum.prs.Control.write()`
        """

        self.write(fname=fname)

    def write(self, fname="raysum-param"):
        """
        Write parameter file to disk

        Args:
            fname: (str)
                Name of file
        """

        with open(fname, "w") as f:
            f.write(self.__repr__())


class Seismogram(object):
    """
    List of streams of 3-component synthetic seismograms produced by Raysum.
    Includes methods to calculate receiver functions, filter and plot the
    streams.

    Parameters:
        model (:class:`~pyraysum.prs.Model`):
            Subsurface velocity model
        geom (:class:`~pyraysum.prs.Geometry`):
            Recording geometry
        rc (:class:`~pyraysum.prs.Control`):
            Run-control parameters
        streams (List):
            List of :class:`~obspy.core.Stream` objects.

    If created with :const:`mode='full'` in :meth:`run()`, the :attr:`stats` attribute
    of each :class:`~obspy.core.Trace` in each :class:`~obspy.core.Stream` holds the
    additional attributes:

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

    def __init__(self, model=None, geom=None, rc=None, streams=None):

        self.model = model
        self.geom = geom
        self.streams = streams
        self.rc = rc

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
        >>> rc = Control(mults=1)
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
            buf += "# Trace number {:5d}\n".format(itr+1)  # Fortran indexing
            buf += "#-------------------\n"
            buf += (
                np.array2string(
                    arr,
                    threshold=3*arr.shape[0] + 1,
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
        :class:`~pyraysum.prs.Seismogram` streams.

        Parameters:
            typ (str):
                Type of plot to show. Options are :const:`'streams'`, :const:`'rfs'`, or
                :const:`'all'` for the displacement seismograms, receiver functions, or
                both
            **kwargs:
                Are passed to the underlying :meth:`~pyraysum.plot.stream_wiggles` or
                :meth:`~pyraysum.plot.rf_wiggles`
        Hint:
            If :class:`Seismogram` contains only one event (i.e., single :attr:`baz` and
            :attr:`slow` values), it is more appropriate to plot the waveforms using the
            default :meth:`Stream.plot()` method on :attr:`Seismogram.streams[0]`.


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
        in :class:`~pyraysum.prs.Seismogram` streams.

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
        :class:`~pyraysum.prs.Seismogram`:
            Synthetic seimograms

    Example
    -------
    >>> from pyraysum import prs, Model, Geometry, Control
    >>> # Define two-layer model with isotropic crust over isotropic half-space
    >>> model = Model([30000., 0], [2800., 3300.], [6000., 8000.], [3600., 4500.])
    >>> geom = Geometry(0., 0.06) # baz = 0 deg; slow = 0.06 x/km
    >>> rc = Control(npts=1500, dt=0.025)
    >>> seismogram = prs.run(model, geom, rc)
    >>> type(seismogram.streams[0])
    <class 'obspy.core.stream.Stream'>
    >>> print(seismogram.streams[0])
    3 Trace(s) in Stream:
    .prs..BHN | 2020-11-30T21:04:43.890339Z - 2020-11-30T21:05:21.365339Z | 40.0 Hz, 1500 samples
    .prs..BHE | 2020-11-30T21:04:43.891418Z - 2020-11-30T21:05:21.366418Z | 40.0 Hz, 1500 samples
    .prs..BHZ | 2020-11-30T21:04:43.891692Z - 2020-11-30T21:05:21.366692Z | 40.0 Hz, 1500 samples
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

    # Store everything into Seismogram object
    seismogram = Seismogram(model=model, geom=geometry.geom, rc=rc, streams=streams)

    if rf:
        seismogram.calculate_rfs()

    return seismogram


def read_model(modfile, encoding=None):
    """
    Reads model parameters from file

    Returns:
        :class:`~pyraysum.prs.Model`:
            Seismic velocity model
    """

    vals = np.genfromtxt(modfile, encoding=encoding).T
    return Model(*vals)


def read_geometry(geomfile, encoding=None):
    """
    Reads geometry parameters from file

    Returns:
        :class:`~pyraysum.prs.Geometry`:
            Ray geometry
    """

    vals = np.genfromtxt(geomfile, dtype=None, encoding=encoding)

    try:
        geom = Geometry(*zip(*vals))
    except TypeError:
        # *values contains single line, which is not iterable
        geom = Geometry([vals[0]], [vals[1]], [vals[2]], [vals[3]])
    return geom


def read_rc(paramfile):
    """
    Read Raysum run control parameters from file.

    Returns:
        :class:`~pyraysum.prs.Control`:
            Run control parameters
    """

    with open(paramfile, "r") as f:
        lines = f.readlines()

    values = []
    for line in lines:
        line = line.strip()
        if line.startswith("#"):
            continue
        values.append(line)
    return Control(*values)


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

    return list(set(ndscrs))


def read_traces(traces, rc, geometry, arrivals=None):
    """
    Create a :class:`Seismogram` from the array produced by :meth:`fraysum.run_bare()`
    and :meth:`fraysum.run_full`.

    Parameters:
        rc (:class:`Control`):
            Run-control parameters
        geometry (:class:`Geometry`):
            Geometry parameters
        arrivals (list):
            Output of :meth:`read_arrivals`. List of arrival times, amplitudes, and
            names

    Returns:
        :class:`~pyraysum.prs.Seismogram`:
            List of Stream objects
    """

    npts = rc.npts
    dt = rc.dt
    rot = rc.rot
    shift = rc.shift
    ntr = geometry.ntr
    geom = geometry.geom

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
                "baz": geom[iitr][0],
                "slow": geom[iitr][1],
                "station": "prs",
                "network": "",
                "starttime": UTCDateTime(0),
                "delta": dt,
                "channel": component[ic],
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


def rfarray(geometry, rc):
    """
    Initialize array for `NumPy`-based post processing

    Parameters:
        geometry (:class:`prs.Geometry`):
            Recording geometry
        rc (:class:`prs.Control`):
            Run-control parameters

    Returns:
        :const:`numpy.zeros((geometry.ntr, 2, rc.npts))`:
            Array in shape to be used by :meth:`filterd_rf_array` and
            :meth:`filtered_array`.
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
    :meth:`Seismogram.calculate_rfs()` and :meth:`Seismogram.filter()`, stripped down
    for computational efficiency, for use in inversion/probabilistic approaches.

    Parameters:
        traces (np.ndarray):
            Output of :meth:`fraysum.run_bare()`
        rfarray (np.ndarray):
            Initialized array of shape (ntr, 2, npts) to store output (See:
            :func:`rfarray()`)
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
            Output is written to :const:`rfarray`

    Warning:
        Assumes PVH alignment (ray-polarization), i.e. :attr:`Control.rot=2`.
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
    :meth:`Seismogram.filter()`, stripped down for computational efficiency, 
    for use in inversion/probabilistic approaches.

    Parameters:
        traces (np.ndarray):
            Output of :meth:`fraysum.run_bare()`
        rfarray (np.ndarray):
            Initialized array of shape (ntr, 2, npts) to store output (See:
            :func:`rfarray()`)
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
            Output is written to :const:`rfarray`

    Warning:
        Assumes PVH alignment (ray-polarization), i.e. :attr:`Control.rot=2`.
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
