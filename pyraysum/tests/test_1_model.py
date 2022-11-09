from pkg_resources import resource_filename
from pyraysum import prs, Model
import pytest
import numpy as np

dipfile = resource_filename('pyraysum',
                           'examples/data/model_Porter2011_dip.txt')
anisofile = resource_filename('pyraysum',
                           'examples/data/model_Porter2011_aniso.txt')


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

def test_getitem_setitem_add():
    model1 = prs.read_model(dipfile)
    assert model1[0]["thickn"] == 20000.
    assert model1[1]["rho"] == 2800.
    assert model1[2]["vp"] == 7800.
    assert model1[2]["vs"] == 4480.
    assert model1[2]["flag"] == 1
    assert model1[1]["dip"] == 20.
    assert model1[2]["strike"] == 90

    model2 = prs.read_model(anisofile)
    assert model2[1]["flag"] == 0
    assert model2[1]["ani"] == -20
    assert model2[1]["trend"] == 180
    assert model2[1]["plunge"] == 45

    model2[1, "plunge"] = 10
    assert model2[1]["plunge"] == 10
    assert model2[1, "plunge"] == 10
    assert model2.plunge[1] == 10
    assert model2.fplunge[1] == 10 * np.pi / 180

    model2[2] = {"plunge": 15, "thickn": 10000}
    assert model2[2, "plunge"] == 15
    assert model2[2, "thickn"] == 10000

    with pytest.raises(ValueError):
        # Only can set user attributes
        model2[1, "fplunge"] = 10

    model3 = model1 + model2
    assert model3[0, "thickn"] == model1[0, "thickn"]
    assert model3[3, "thickn"] == model2[0, "thickn"]
    assert model3[5, "thickn"] == model2[2, "thickn"]

    model4 = model1 + {"thickn": 5000, "rho": 3000, "vp": 6000}
    assert model4[3, "thickn"] == 5000

    model5 = model1 + [5000, 3000, 6000]
    assert model5[3, "thickn"] == 5000

    with pytest.raises(TypeError):
        # Only can add another model
        model1 + 5
        

    
    
