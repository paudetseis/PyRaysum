.. figure:: ../pyraysum/examples/picture/PyRaysum_logo.png
   :align: center

Synthetic seismograms are computed from three classes: :class:`Model`,
:class:`Geometry`, and :class:`RC` using :func:`prs.run`, and are stored in
an instance of class :class:`Seismogram`.

Classes
=======

:mod:`~pyraysum` defines the following classes:

- :class:`~pyraysum.prs.Model`
- :class:`~pyraysum.prs.Geometry`
- :class:`~pyraysum.prs.Seismogram`
- :class:`~pyraysum.prs.RC`

Model
-----

.. autoclass:: pyraysum.prs.Model
   :members:

Geometry
--------

.. autoclass:: pyraysum.prs.Geometry
   :members:

RC
--

.. autoclass:: pyraysum.prs.RC
   :members:

Seismogram
----------

.. autoclass:: pyraysum.prs.Seismogram
   :members:

Modules
=======

Module :mod:`prs`
-----------------

.. automodule:: pyraysum.prs
   :members:
   :exclude-members: Model, Geometry, RC, Seismogram

Module :mod:`plot`
------------------

.. automodule:: pyraysum.plot
	:members:

Module :mod:`fraysum`
---------------------
These are the the access points to the low-level fortran routines.

.. tip::
   The corresponding functions can be more succinctly called with the unpacked 
   lists of :attr:`parameters`, like this:
   ``fraysum.run_bare(*model.parameters, *geometry.parameters, *rc.parameters)``.

.. py:function:: tr_ph = fraysum.run_bare(thick, rho, alpha, beta, isoflag, pct, trend,\
    plunge, strike, dip, nlay, baz, slow, sta_dx, sta_dy, ntr, iphase, mults, nsamp,\
    dt, align, shift, out_rot, verb, nsegin, numphin, phaselistin)

    :param thick: :attr:`Model.fthickn`
    :param rho: :attr:`Model.frho`
    :param alpha: :attr:`Model.fvp`
    :param beta: :attr:`Model.fvs`
    :param isoflag: :attr:`Model.fflag`
    :param pct: :attr:`Model.fani`
    :param trend: :attr:`Model.ftrend`
    :param plunge: :attr:`Model.fplunge`
    :param strike: :attr:`Model.fstrike`
    :param dip: :attr:`Model.fdip`
    :param nlay: :attr:`Model.nlay`
    :param baz: :attr:`Geometry.fbaz`
    :param slow: :attr:`Geometry.fslow`
    :param sta_dx: :attr:`Geometry.fdn`
    :param sta_dy: :attr:`Geometry.fde`
    :param ntr: :attr:`Geometry.ntr`
    :param iphase: :attr:`RC.wvtype`
    :param mults: :attr:`RC.mults`
    :param nsamp: :attr:`RC.nsamp`
    :param dt: :attr:`RC.dt`
    :param align: :attr:`RC.align`
    :param shift: :attr:`RC.shift`
    :param out_rot: :attr:`RC.rot`
    :param verb: :attr:`RC.verbose`
    :param nsegin: :attr:`RC._nseg`
    :param numphin: :attr:`RC._numph`
    :param phaselistin: :attr:`RC.phaselist`

    :returns: tr_ph
    :rtype: rank-3 array('f') with bounds (3, `maxsamp`, `maxtr`)

.. py:function:: tr_ph, travel_time, amplitudeout, phaselist = fraysum.run_full(thick,\
    rho, alpha, beta, isoflag, pct, trend, plunge, strike, dip, nlay, baz, slow,\
    sta_dx, sta_dy, ntr, iphase, mults, nsamp, dt, align, shift, out_rot, verb, nsegin,\
    numphin, phaselistin)

    :param thick: :attr:`Model.fthickn`
    :param rho: :attr:`Model.frho`
    :param alpha: :attr:`Model.fvp`
    :param beta: :attr:`Model.fvs`
    :param isoflag: :attr:`Model.fflag`
    :param pct: :attr:`Model.fani`
    :param trend: :attr:`Model.ftrend`
    :param plunge: :attr:`Model.fplunge`
    :param strike: :attr:`Model.fstrike`
    :param dip: :attr:`Model.fdip`
    :param nlay: :attr:`Model.nlay`
    :param baz: :attr:`Geometry.fbaz`
    :param slow: :attr:`Geometry.fslow`
    :param sta_dx: :attr:`Geometry.fdn`
    :param sta_dy: :attr:`Geometry.fde`
    :param ntr: :attr:`Geometry.ntr`
    :param iphase: :attr:`RC.wvtype`
    :param mults: :attr:`RC.mults`
    :param nsamp: :attr:`RC.nsamp`
    :param dt: :attr:`RC.dt`
    :param align: :attr:`RC.align`
    :param shift: :attr:`RC.shift`
    :param out_rot: :attr:`RC.rot`
    :param verb: :attr:`RC.verbose`
    :param nsegin: :attr:`RC._nseg`
    :param numphin: :attr:`RC._numph`
    :param phaselistin: :attr:`RC.phaselist`

    :returns: tr_ph
    :rtype: rank-3 array('f') with bounds (3, `maxsamp`, `maxtr`):
    :returns: travel_time
    :rtype: rank-2 array('f') with bounds (`maxph`, `maxtr`)
    :returns: amplitudeout
    :rtype: rank-3 array('f') with bounds (3, `maxph`,`maxtr`)
    :returns: phaselist
    :rtype: rank-3 array('i') with bounds (`maxseg`, 2, `maxph`)
