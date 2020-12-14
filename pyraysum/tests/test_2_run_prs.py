from pyraysum import prs, Model
from . import test_1_model as mod
import numpy as np
import pytest

def test_Porter2011():
    # Define range of slowness and back-azimuths
    baz = np.arange(0., 360., 10.)
    slow = 0.06

    # Read first model with dipping lower crustal layer
    model = mod.test_read_model_dip()

    # Run Raysum with most default values and `rot=1` and `mults=0`
    # to reproduce the results of Porter et al., 2011
    streamlist = prs.run_prs(model, baz, slow, rot=1, mults=0)

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

    streamlist = prs.run_prs(model, baz, slow, rot=1, mults=0)
    streamlist.calculate_rfs()
    streamlist.filter('all', 'lowpass', freq=1., zerophase=True, corners=2)
    streamlist.plot('all', tmin=-0.5, tmax=8.)

def test_single_event():

    model = mod.test_def_model()
    slow = 0.06     # s/km
    baz = 0.
    npts = 1500
    dt = 0.025      # s
    streamlist = prs.run_prs(model, baz, slow, npts=npts, dt=dt, rot=2)
    with pytest.raises(Exception):
        assert prs.run_prs(model, baz, slow, npts=npts, dt=dt, rot=3)
    streamlist = prs.run_prs(model, baz, slow, npts=npts, dt=dt, rot=1, wvtype='SV')
    streamlist = prs.run_prs(model, baz, slow, npts=npts, dt=dt, rot=1, wvtype='SH')


def test_rfs():

    model = mod.test_def_model()
    slow = 0.06     # s/km
    baz = 0.
    npts = 1500
    dt = 0.025      # s

    # test 1
    with pytest.raises(Exception):
        assert prs.run_prs(model, baz, slow, npts=npts, dt=dt, rot=0, rf=True)
    streamlist1 = prs.run_prs(model, baz, slow, npts=npts, dt=dt, rot=1, rf=True)
    streamlist1.filter('rfs', 'lowpass', freq=1., corners=2, zerophase=True)
    streamlist1.filter('streams', 'lowpass', freq=1., corners=2, zerophase=True)

    # test 2
    streamlist2 = prs.run_prs(model, baz, slow, npts=npts, dt=dt)
    with pytest.raises(Exception):
        assert streamlist2.calculate_rfs()

    streamlist2 = prs.run_prs(model, baz, slow, npts=npts, dt=dt, rot=1)
    rflist = streamlist2.calculate_rfs()
    [rf.filter('lowpass', freq=1., corners=2, zerophase=True) for rf in rflist]
