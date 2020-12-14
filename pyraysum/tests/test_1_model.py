import numpy as np
from pkg_resources import resource_filename
from pyraysum import prs, Model

dipfile = resource_filename('pyraysum',
                           'examples/models/model_Porter2011_dip.txt')
anisofile = resource_filename('pyraysum',
                           'examples/models/model_Porter2011_aniso.txt')


def test_read_model_dip():
    model = prs.read_model(dipfile)
    assert isinstance(model, Model)
    return model

def test_read_model_aniso():
    model = prs.read_model(anisofile)
    assert isinstance(model, Model)
    return model

def test_def_model():
    thick = [20000., 10000., 0.]
    rho = [2800., 2950., 3300.]
    vp = [4600., 5000., 6000.]
    vs = [2600., 3000., 3600.]
    flag = [1, 0, 1]
    ani = [0, 5, 0]
    model = Model(thick, rho, vp, vs, flag=flag, ani=ani)
    assert isinstance(model, Model)
    return model

def test_plot_model():
    model = test_def_model()
    model.plot()

