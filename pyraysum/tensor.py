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

Utility functions to interact with ``pyraysum`` module.

'''
import itertools
import numpy as np
from numpy import sin, cos
from pyraysum import elast


def set_iso_tensor(a, b):
    """
    Function to generate tensor for isotropic material.

    Args:
        a (float): P-wave velocity (km/s)
        b (float): S-wave velocity (km/s)

    Returns:
        (np.ndarray): cc: Elastic tensor (GPa /density) \
        (shape ``(3, 3, 3, 3)``)

    """

    a = a*1.e3
    b = b*1.e3
    C = elastiso_tensor(a, b)

    # Convert Voigt to full tensor
    cc = voigt2cc(C)

    return cc


def set_tri_tensor(a, b, tr, pl, ani):
    """
    Function to generate tensor for transverse isotropy. The
    tensor is rotated using the trend and plunge of the symmetry
    axis.

    Args:
        a (float): P-wave velocity (km/s)
        b (float): S-wave velocity (km/s)
        tr (float): Trend angle of symmetry axis (degree)
        pl (float): Plunge angle of symmetry axis (degree)
        ani (float): Percent anisotropy

    Returns:
        (np.ndarray): cc: Elastic tensor (GPa / kg/m^3) \
        (shape ``(3, 3, 3, 3)``)

    """

    # Trend and plunge of symmetry axis
    tr = -tr*np.pi/180.
    pl = (90. - pl)*np.pi/180.

    # Percent anisotropy
    da = (a*1.e3)*ani/100.
    db = (b*1.e3)*ani/100.

    # Set up matrix elements
    AA = (a*1.e3 - da/2.)**2
    CC = (a*1.e3 + da/2.)**2
    LL = (b*1.e3 + db/2.)**2
    NN = (b*1.e3 - db/2.)**2
    AC = (a*1.e3)**2
    FF = -LL + np.sqrt((2.*AC)**2 - 2.*AC*(AA + CC + 2.*LL) +
                       (AA + LL)*(CC + LL))
    # eta = FF/(AA - 2.*LL)

    # Get tensor with horizontal axis
    cc = np.zeros((3, 3, 3, 3))

    cc[0, 0, 0, 0] = AA
    cc[1, 1, 1, 1] = AA
    cc[2, 2, 2, 2] = CC

    cc[0, 0, 1, 1] = (AA - 2.*NN)
    cc[1, 1, 0, 0] = (AA - 2.*NN)

    cc[0, 0, 2, 2] = FF
    cc[2, 2, 0, 0] = FF

    cc[1, 1, 2, 2] = FF
    cc[2, 2, 1, 1] = FF

    cc[1, 2, 1, 2] = LL
    cc[2, 1, 2, 1] = LL
    cc[2, 1, 1, 2] = LL
    cc[1, 2, 2, 1] = LL

    cc[2, 0, 2, 0] = LL
    cc[0, 2, 0, 2] = LL
    cc[2, 0, 0, 2] = LL
    cc[0, 2, 2, 0] = LL

    cc[0, 1, 0, 1] = NN
    cc[1, 0, 1, 0] = NN
    cc[0, 1, 1, 0] = NN
    cc[1, 0, 0, 1] = NN

    # Rotate tensor using trend and plunge
    cc = rot_tensor(cc, pl, tr, 0.)

    # Return tensor
    return cc


def set_aniso_tensor(tr, pl, typ='atg'):
    """
    Function to generate tensor for anisotropic minerals. The \
    tensor is rotated using the trend and plunge of the symmetry \
    axis.

    Args:
        tr (float): Trend angle of symmetry axis (degree)
        pl (float): Plunge angle of symmetry axis (degree)
        type (str, optional): Type of elastic material

    Returns:
        (tuple): Tuple containing:
            * cc (np.ndarray): Elastic tensor (GPa / kg/m^3)\
            (shape ``(3, 3, 3, 3)``)
            * rho (float): Density (kg/m^3)

    """

    # Trend and plunge of symmetry axis
    tr = -tr*np.pi/180.
    pl = (90. - pl)*np.pi/180.

    # Get tensor with horizontal axis

    # Minerals
    if typ == 'atg':
        C, rho = elastantigorite()
    elif typ == 'bt':
        C, rho = elastbiotite()
    elif typ == 'cpx':
        C, rho = elastclinopyroxene_92()
    elif typ == 'dol':
        C, rho = elastdolomite()
    elif typ == 'ep':
        C, rho = elastepidote()
    elif typ == 'grt':
        C, rho = elastgarnet()
    elif typ == 'gln':
        C, rho = elastglaucophane()
    elif typ == 'hbl':
        C, rho = elasthornblende()
    elif typ == 'jade':
        C, rho = elastjadeite()
    elif typ == 'lws':
        C, rho = elastlawsonite()
    elif typ == 'lz':
        C, rho = elastlizardite()
    elif typ == 'ms':
        C, rho = elastmuscovite()
    elif typ == 'ol':
        C, rho = elastolivine()
    elif typ == 'opx':
        C, rho = elastorthopyroxene()
    elif typ == 'plag':
        C, rho = elastplagioclase_06()
    elif typ == 'qtz':
        C, rho = elastquartz()
    elif typ == 'zo':
        C, rho = elastzoisite()

    # Rocks
    elif typ == 'BS_f':
        C, rho = elastblueschist_felsic()
    elif typ == 'BS_m':
        C, rho = elastblueschist_mafic()
    elif typ == 'EC_f':
        C, rho = elasteclogite_foliated()
    elif typ == 'EC_m':
        C, rho = elasteclogite_massive()
    elif typ == 'HB':
        C, rho = elastharzburgite()
    elif typ == 'SP_37':
        C, rho = elastserpentinite_37()
    elif typ == 'SP_80':
        C, rho = elastserpentinite_80()
    elif typ == 'LHZ':
        C, rho = elastlherzolite()

    else:
        print('type of mineral/rock not implemented')
        return

    # Convert Voigt to full tensor
    cc = voigt2cc(C)*1.e9/rho

    # Rotate tensor using trend and plunge
    cc = rot_tensor(cc, pl, tr, 0.)

    # Return tensor
    return cc, rho


def full_3x3_to_Voigt_6_index(i, j):
    """
    Conversion of tensor to Voigt notation for indices

    """
    if i == j:
        return i
    return 6-i-j


def voigt2cc(C):
    """
    Convert the Voigt representation of the stiffness matrix to the full
    3x3x3x3 tensor representation.

    Args:
        C (np.ndarray): Stiffness matrix (shape ``(6, 6)``)

    Returns:
        (np.ndarray): cc: Elastic tensor (shape ``(3, 3, 3, 3)``)

    """

    C = np.asarray(C)
    cc = np.zeros((3, 3, 3, 3), dtype=float)
    for i, j, k, l in itertools.product(range(3), repeat=4):
        Voigt_i = full_3x3_to_Voigt_6_index(i, j)
        Voigt_j = full_3x3_to_Voigt_6_index(k, l)
        cc[i, j, k, l] = C[Voigt_i, Voigt_j]
    return cc


def cc2voigt(cc):
    """
    Convert from the full 3x3x3x3 tensor representation
    to the Voigt notation of the stiffness matrix.

    Args:
        cc (np.ndarray): Elastic tensor (shape ``(3, 3, 3, 3)``)

    Returns:
        (np.ndarray): C: Stiffness matrix (shape ``(6, 6)``)

    """

    Voigt_notation = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]

    cc = np.asarray(cc)
    C = np.zeros((6, 6))
    for i in range(6):
        for j in range(6):
            k, l = Voigt_notation[i]
            m, n = Voigt_notation[j]
            C[i, j] = cc[k, l, m, n]

    return C


def VRH_average(C):
    """
    Performs a Voigt-Reuss-Hill average of the anisotropic
    stifness matrix to the bulk modulus K and the shear modulus
    G.

    Args:
        C (np.ndarray): Stiffness matrix (shape ``(6, 6)``)

    Returns:
        (tuple): Tuple containing:
            * Kvoigt (float): Voigt average bulk modulus (GPa)
            * Gvoigt (float): Voigt average shear modulus (GPa)
            * Kreuss (float): Reuss average bulk modulus (GPa)
            * Greuss (float): Reuss average shear modulus (GPa)
            * Kvrh (float):   Voigt-Reuss-Hill average bulk modulus (GPa)
            * Gvrh (float):   Voigt-Reuss-Hill average shear modulus (GPa)

    Example
    -------
    >>> from telewavesim import utils
    >>> cc, rho = utils.set_aniso_tensor(0., 0., typ='atg')
    >>> C = utils.cc2voigt(cc)
    >>> utils.VRH_average(C*rho)
    (75655555555.555557, 48113333333.333336, 61245706544.967415,
    28835098086.844658, 68450631050.26149, 38474215710.088997)

    """

    # Compliance matrix
    S = np.linalg.inv(C)

    # Voigt averaging
    Kvoigt = (C[0, 0] + C[1, 1] + C[2, 2] + 2 * C[0, 1] + 2 * C[0, 2] +
              2 * C[1, 2]) / 9
    Gvoigt = (C[0, 0] + C[1, 1] + C[2, 2] - C[0, 1] - C[0, 2] - C[1, 2] +
              3 * C[3, 3] + 3 * C[4, 4] + 3 * C[5, 5]) / 15

    # Reuss averaging
    Kreuss = 1 / (S[0, 0] + S[1, 1] + S[2, 2] + 2 * S[0, 1] + 2 * S[0, 2] +
                  2 * S[1, 2])
    Greuss = 15 / (4 * S[0, 0] + 4 * S[1, 1] + 4 * S[2, 2] - 4 * S[0, 1] -
                   4 * S[0, 2] - 4 * S[1, 2] + 3 * S[3, 3] + 3 * S[4, 4] +
                   3 * S[5, 5])

    # Voigt-Reuss-Hill average
    Kvrh = (Kvoigt + Kreuss) / 2
    Gvrh = (Gvoigt + Greuss) / 2

    return Kvoigt, Gvoigt, Kreuss, Greuss, Kvrh, Gvrh


def rot_tensor(a, alpha, beta, gam):
    """

    Performs a rotation of the tensor cc (c_ijkl) about three angles (alpha,
    beta, gamma)

    Args:
        a (np.ndarray): Elastic tensor with shape ``(3, 3, 3, 3)``
        alpha (float): Angle in radians
        beta (float): Angle in radians
        gam (float): Angle in radians

    Returns:
        (np.ndarray): aa: Rotated tensor with shape ``(3, 3, 3, 3)``

    .. note::

        The three angles (``alpha``, ``beta``, ``gam``) correspond to rotation
        about the x_2, x_3, x_1 axelast Note that the sequence of the rotation
        is important: (AB ~= BA). In this case we rotate about x_2 first,
        x_3 second and x_1 third.

        For trend and plunge of symmetry axis (e.g., tri_tensor):

            ``alpha`` = plunge

            ``beta`` = trend

    """

    rot = np.zeros((3, 3))
    aa = np.zeros((3, 3, 3, 3))

    rot[0, 0] = cos(alpha) * cos(beta)
    rot[0, 1] = sin(beta)
    rot[0, 2] = sin(alpha) * cos(beta)

    rot[1, 0] = -cos(gam) * sin(beta) * cos(alpha) - sin(gam) * sin(alpha)
    rot[1, 1] = +cos(gam) * cos(beta)
    rot[1, 2] = -cos(gam) * sin(beta) * sin(alpha) + sin(gam) * cos(alpha)
    rot[2, 0] = +sin(gam) * sin(beta) * cos(alpha) - cos(gam) * sin(alpha)
    rot[2, 1] = -sin(gam) * cos(beta)
    rot[2, 2] = +sin(gam) * sin(beta) * sin(alpha) + cos(gam) * cos(alpha)

    #
    #  c_ijkl ---> c_mnrs
    #
    for m, n, r, s in itertools.product(range(3), repeat=4):
        asum = 0.0
        for i, j, k, l in itertools.product(range(3), repeat=4):
            rr = rot[m, i] * rot[n, j] * rot[r, k] * rot[s, l]
            asum = asum + rr * a[i, j, k, l]
        aa[m, n, r, s] = asum

    return aa


def mod2vel(K, G, rho):
    """

    Calculates the isotropic P and S wave velocities from given
    bulk (K) and shear (G) moduli and density (rho) in kg/m^3

    Args:
        K (float): Bulk modulus (GPa)
        G (float): Shear modulus (GPa)
        rho (float): Density (kg/m^3)
    Returns:
        (tuple): tuple containing:

            * Vp (float): P-wave velocity (m/s)
            * Vs (float): S-wave velocity (m/s)

    Example
    -------
    >>> from telewavesim import utils
    >>> cc, rho = utils.set_aniso_tensor(0., 0., typ='atg')
    >>> C = utils.cc2voigt(cc)
    >>> K, G = utils.VRH_average(C*rho)[4:6]
    >>> utils.mod2vel(K, G, rho)
    (6760.617471753726, 3832.0771334254896)

    """

    Vp = np.sqrt((K + 4.*G/3.)/rho)
    Vs = np.sqrt(G/rho)

    return Vp, Vs
