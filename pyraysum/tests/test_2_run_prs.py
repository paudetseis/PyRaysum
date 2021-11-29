from pyraysum import prs, Model, Geometry
from fraysum import call_seis_spread
import numpy as np
import pytest
import matplotlib.pyplot as mp

try:
    #  pytest
    from . import test_1_model as mod
except ImportError:
    #  ipython
    import test_1_model as mod

def test_Porter2011():
    # Define range of slowness and back-azimuths
    baz = np.arange(0., 360., 10.)
    slow = 0.06

    # Read first model with dipping lower crustal layer
    model = mod.test_read_model_dip()
    geom = Geometry(baz, slow)

    print('model', model._dict__)
    print('geom', geom.__dict__)
    print('running run_frs')
    # Run Raysum with most default values and `rot=1` and `mults=0`
    # to reproduce the results of Porter et al., 2011
    streamlist = prs.run_frs(model, geom, rot=1, mults=0)

    # Calculate receiver functions
    streamlist.calculate_rfs()

    # Filter and plot receiver functions with argument `rfs` (other option
    # is 'streams')
    streamlist.filter('streams', 'lowpass', freq=1., zerophase=True, corners=2)
    streamlist.filter('rfs', 'lowpass', freq=1., zerophase=True, corners=2)
    streamlist.plot('streams', tmin=-0.5, tmax=8.)
    streamlist.plot('rfs', tmin=-0.5, tmax=8.)

    # Now load a different model and repeat (lower crustal anisotropic layer)
    model = mod.test_read_model_aniso()

    streamlist = prs.run_frs(model, geom, rot=1, mults=0)
    streamlist.calculate_rfs()
    streamlist.filter('all', 'lowpass', freq=1., zerophase=True, corners=2)
    streamlist.plot('all', tmin=-0.5, tmax=8.)

def test_frs():
    # Define range of slowness and back-azimuths
    baz = np.arange(0., 360., 10.)
    slow = 0.06

    # Read first model with dipping lower crustal layer
    model = mod.test_read_model_dip()
    geom = Geometry(baz, slow)

    # Run Raysum with most default values and `rot=1` and `mults=0`
    # to reproduce the results of Porter et al., 2011
    fstreamlist = prs.run_frs(model, geom, rot=1, mults=0)

    # Now load a different model and repeat (lower crustal anisotropic layer)
    model = mod.test_read_model_aniso()

    fstreamlist = prs.run_frs(model, geom, rot=2, mults=2)

    fstreamlist.calculate_rfs()


def test_filtered_rf_array():
    from timeit import timeit
    # Define range of slowness and back-azimuths
    baz = range(0, 360, 15)
    slow = 0.06 
    dt = 0.05
    rot = 2
    mults = 0
    verbose = 0
    wvtype = 'P'
    align = 2
    npts = 1000
    fmin = 0.05
    fmax = 0.5

    # Read first model with dipping lower crustal layer
    model = mod.test_read_model_dip()
    geom = Geometry(baz, slow)
    rfarray = np.zeros((geom.ntr, 2, npts))


    # To time this, do:
    # print(timeit('_run_frs()', number=50, globals=globals()))
    # >> 27.262143349274993
    def _run_frs():
        streams = prs.run_frs(model, geom, dt=dt, rot=rot, mults=mults,
                              align=align, wvtype=wvtype, verbose=verbose,
                              npts=npts)
        streams.calculate_rfs()
        streams.filter('rfs', 'bandpass', freqmin=fmin, freqmax=fmax,
                       zerophase=True, corners=2)
        return streams

    # To time this, do:
    # print(timeit('_run_sspread()', number=50, globals=globals()))
    # >> 16.743131840601563
    def _run_sspread():
        tr_ph, _ = call_seis_spread(
                model.fthickn, model.frho, model.fvp, model.fvs, model.fflag,
                model.fani, model.ftrend, model.fplunge, model.fstrike, model.fdip,
                model.nlay,
                geom.fbaz, geom.fslow, geom.fdx, geom.fdy, geom.ntr,
                wvtype, mults, npts, dt, align, dt, rot, verbose)

        prs.filtered_rf_array(tr_ph, rfarray, geom.ntr, npts, dt, fmin, fmax)

    streams = _run_frs()

    _run_sspread()

    fig, ax = mp.subplots(len(geom.geom), 2, tight_layout=True)
    for l, (stream, array)  in enumerate(zip(streams.rfs, rfarray)):
        for r, (scomp, acomp) in enumerate(zip(stream, array)):
            ax[l, r].plot(scomp.data)
            ax[l, r].plot(acomp)
    fig.show()
    #input('Press key to continue')
    


def test_single_event():

    model = mod.test_def_model()
    slow = 0.06     # s/km
    baz = 0.
    npts = 1500
    dt = 0.025      # s
    geom = Geometry(baz, slow)
    streamlist = prs.run_frs(model, geom, npts=npts, dt=dt, rot=2)
    with pytest.raises(Exception):
        assert prs.run_frs(model, geom, npts=npts, dt=dt, rot=3)
    streamlist = prs.run_frs(model, geom, npts=npts, dt=dt, rot=1, wvtype='SV')
    streamlist = prs.run_frs(model, geom, npts=npts, dt=dt, rot=1, wvtype='SH')


def test_rfs():

    model = mod.test_def_model()
    slow = 0.06     # s/km
    baz = 0.
    npts = 1500
    dt = 0.025      # s

    geom = Geometry(baz, slow)

    # test 1
    with pytest.raises(Exception):
        assert prs.run_frs(model, geom, npts=npts, dt=dt, rot=0, rf=True)
    streamlist1 = prs.run_frs(model, geom, npts=npts, dt=dt, rot=1, rf=True)
    streamlist1 = prs.run_frs(model, geom, npts=npts, dt=dt, rot=2, rf=True)
    streamlist1.filter('rfs', 'lowpass', freq=1., corners=2, zerophase=True)
    streamlist1.filter('streams', 'lowpass', freq=1., corners=2, zerophase=True)

    # test 2
    streamlist2 = prs.run_frs(model, geom, npts=npts, dt=dt)
    with pytest.raises(Exception):
        assert streamlist2.calculate_rfs()

    streamlist2 = prs.run_frs(model, geom, npts=npts, dt=dt, rot=1)
    rflist = streamlist2.calculate_rfs()
    [rf.filter('lowpass', freq=1., corners=2, zerophase=True) for rf in rflist]
