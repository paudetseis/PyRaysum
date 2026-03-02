from pyraysum import Geometry, Control, run
from pyraysum import prs, frs
from fraysum import run_bare
from obspy import Stream
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
    rc = Control(rot=1, mults=0)

    print('model', model.__dict__)
    print('geom', geom.__dict__)
    print('running run()')
    # Run Raysum with most default values and `rot=1` and `mults=0`
    # to reproduce the results of Porter et al., 2011
    streamlist = run(model, geom, rc)

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

    streamlist = run(model, geom, rc)
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
    rc = Control(rot=1, mults=0)

    # Run Raysum with most default values and `rot=1` and `mults=0`
    # to reproduce the results of Porter et al., 2011
    fstreamlist = run(model, geom, rc)

    # Now load a different model and repeat (lower crustal anisotropic layer)
    model = mod.test_read_model_aniso()

    fstreamlist = run(model, geom, rc)
    assert len(fstreamlist[0][0]) == 3  # 3-components
    assert len(fstreamlist[0][1]) == 0  # No RFs
    assert len(fstreamlist) == len(geom)

    fstreamlist.calculate_rfs()
    assert len(fstreamlist[0][1]) == 2  # No RFs


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
    rc = Control(dt=dt, rot=rot, mults=mults, align=align, wvtype=wvtype, verbose=verbose,
            npts=npts)

    rfarray = frs.make_array(geom, rc)

    # To time this, do:
    # print(timeit('_run()', number=50, globals=globals()))
    # >> 21.00839215517044
    def _run():
        streams = run(model, geom, rc)
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

        frs.filtered_rf_array(ph_traces, rfarray, geom.ntr, rc.npts, rc.dt, fmin, fmax)

    streams = _run()

    _run_bare()

    fig, ax = mp.subplots(len(geom._geom), 2, tight_layout=True)
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
        rc = Control(npts=1500, dt=0.025, rot=2, wvtype=wvtype)
        result = run(model, geom, rc)

    assert isinstance(result[0][0], Stream)
    assert len(result[0][0]) == 3  # 3-components
    assert len(result[0][1]) == 0  # No RFs
    assert len(result) == 1
    assert len(result["seis"]) == 1
    assert len(result["streams"]) == 1
    assert len(result["seismograms"]) == 1

def test_rfs():

    model = mod.test_def_model()
    rc = Control(npts=1500, dt=0.025, rot=0)
    geom = Geometry(slow=0.06, baz=0.)

    # test 1
    with pytest.raises(ValueError):
        assert run(model, geom, rc, rf=True)

    rc.rot = 1
    streamlist1 = run(model, geom, rc, rf=True)

    rc.rot = 2
    streamlist1 = run(model, geom, rc, rf=True)
    streamlist1.filter('rfs', 'lowpass', freq=1., corners=2, zerophase=True)
    streamlist1.filter('streams', 'lowpass', freq=1., corners=2, zerophase=True)

    # test 2
    rc.rot = 0
    streamlist2 = run(model, geom, rc)
    with pytest.raises(ValueError):
        assert streamlist2.calculate_rfs()

    rc.rot = 1
    streamlist2 = run(model, geom, rc)
    rflist = streamlist2.calculate_rfs()
    [rf.filter('lowpass', freq=1., corners=2, zerophase=True) for rf in streamlist2.rfs]

    assert isinstance(streamlist2[0][1], Stream)
    assert len(streamlist2[0][1]) == 2
    assert len(streamlist2["rfs"]) == len(geom)
    assert len(streamlist2["rf"]) == len(geom)

def test_bailout():
    """
    The given model produces the 'WARNING in evec_check'. Make sure phases are bailed out
    """
    from pkg_resources import resource_filename

    modf = resource_filename('pyraysum',
                               'tests/evec_warn_model.txt')
    mod = prs.read_model(modf)
    # modf = "evec_warn_model.txt"
    # try:
    #     mod = prs.read_model(modf)
    # except OSError:
    #     # VSCode starts tests from root directory
    #     modf = "pyraysum/tests/" + modf
    #     mod = prs.read_model(modf)
        
    rc = Control(rot=2)

    geom1 = Geometry(7.5, 0.045)

    seis = run(mod, geom1, rc)

    amps = seis.streams[0][1].stats.phase_amplitudes

    assert all(abs(amps) < 2.5)

    geom2 = Geometry(135, 0.06)
    seis = run(mod, geom2, rc)

    assert all(abs(amps) < 2.5)
