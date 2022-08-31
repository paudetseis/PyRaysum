from pyraysum import prs, Model, Geometry, RC
from fraysum import run_full, run_bare
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
    rc = RC(rot=1, mults=0)

    print('model', model.__dict__)
    print('geom', geom.__dict__)
    print('running run()')
    # Run Raysum with most default values and `rot=1` and `mults=0`
    # to reproduce the results of Porter et al., 2011
    streamlist = prs.run(model, geom, rc)

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

    streamlist = prs.run(model, geom, rc)
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
    rc = RC(rot=1, mults=0)

    # Run Raysum with most default values and `rot=1` and `mults=0`
    # to reproduce the results of Porter et al., 2011
    fstreamlist = prs.run(model, geom, rc)

    # Now load a different model and repeat (lower crustal anisotropic layer)
    model = mod.test_read_model_aniso()

    fstreamlist = prs.run(model, geom, rc)

    fstreamlist.calculate_rfs()


def test_filtered_rf_array():
    from timeit import timeit
    # Define range of slowness and back-azimuths
    baz = range(0, 360, 15)
    slow = 0.06 
    dt = 0.05
    rot = 2
    mults = 1
    verbose = False
    wvtype = 'P'
    align = 2
    npts = 1000
    fmin = 0.05
    fmax = 0.5

    # Read first model with dipping lower crustal layer
    model = mod.test_read_model_dip()
    geom = Geometry(baz, slow)
    rc = RC(dt=dt, rot=rot, mults=mults, align=align, wvtype=wvtype, verbose=verbose,
            npts=npts)

    rfarray = np.zeros((geom.ntr, 2, npts))

    # To time this, do:
    # print(timeit('_run()', number=50, globals=globals()))
    # >> 21.00839215517044
    def _run():
        streams = prs.run(model, geom, rc)
        streams.calculate_rfs()
        streams.filter('rfs', 'bandpass', freqmin=fmin, freqmax=fmax,
                       zerophase=True, corners=2)
        return streams

    # To time this, do:
    # print(timeit('_run_bare()', number=50, globals=globals()))
    # >> 9.778649725019932
    def _run_bare():
        ph_traces = run_bare(
            *model.parameters, *geom.parameters, *rc.parameters)

        prs.filtered_rf_array(ph_traces, rfarray, geom.ntr, rc.npts, rc.dt, fmin, fmax)

    streams = _run()

    _run_bare()

    fig, ax = mp.subplots(len(geom.geom), 2, tight_layout=True)
    for l, (stream, array)  in enumerate(zip(streams.rfs, rfarray)):
        for r, (scomp, acomp) in enumerate(zip(stream, array)):
            ax[l, r].plot(scomp.data)
            ax[l, r].plot(acomp)
    fig.show()
    #input('Press key to continue')
    

def test_single_event():

    model = mod.test_def_model()
    geom = Geometry(baz=0, slow=0.06)
    for wvtype in ['P', 'SV', 'SH']:
        rc = RC(npts=1500, dt=0.025, rot=2, wvtype=wvtype)
        streamlist = prs.run(model, geom, rc)

def test_rfs():

    model = mod.test_def_model()
    rc = RC(npts=1500, dt=0.025, rot=0)
    geom = Geometry(slow=0.06, baz=0.)

    # test 1
    with pytest.raises(ValueError):
        assert prs.run(model, geom, rc, rf=True)

    rc.rot = 1
    streamlist1 = prs.run(model, geom, rc, rf=True)

    rc.rot = 2
    streamlist1 = prs.run(model, geom, rc, rf=True)
    streamlist1.filter('rfs', 'lowpass', freq=1., corners=2, zerophase=True)
    streamlist1.filter('streams', 'lowpass', freq=1., corners=2, zerophase=True)

    # test 2
    rc.rot = 0
    streamlist2 = prs.run(model, geom, rc)
    with pytest.raises(ValueError):
        assert streamlist2.calculate_rfs()

    rc.rot = 1
    streamlist2 = prs.run(model, geom, rc)
    rflist = streamlist2.calculate_rfs()
    [rf.filter('lowpass', freq=1., corners=2, zerophase=True) for rf in streamlist2.rfs]
