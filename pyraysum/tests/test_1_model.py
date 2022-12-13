from pkg_resources import resource_filename
from pyraysum import prs, Model
import pytest
import tempfile
from os import remove
import numpy as np

dipfile = resource_filename("pyraysum", "examples/data/model_Porter2011_dip.txt")
anisofile = resource_filename("pyraysum", "examples/data/model_Porter2011_aniso.txt")


def test_read_model_dip():
    model = prs.read_model(dipfile)
    assert isinstance(model, Model)
    return model


def test_read_model_aniso():
    model = prs.read_model(anisofile)
    assert isinstance(model, Model)
    return model


def test_write_read_model_v1_2():
    """Write and read a raysum v1.2 model"""
    model = prs.read_model(anisofile)
    tmpf = tempfile.NamedTemporaryFile(delete=False)
    model.write(tmpf.name, version="1.2")
    model = prs.read_model(tmpf.name, version="1.2")
    remove(tmpf.name)
    assert isinstance(model, Model)
    assert model[1, "ani"] == -20.0


def test_def_model():
    thick = [20000.0, 10000.0, 0.0]
    rho = [2800.0, 2950.0, 3300.0]
    vp = [4600.0, 5000.0, 6000.0]
    vs = [2600.0, 3000.0, 3600.0]
    flag = [1, 0, 1]
    ani = [0, 5, 0]
    model = Model(thick, rho, vp, vs, flag=flag, ani=ani)
    assert isinstance(model, Model)
    return model


def test_plot_model():
    model = test_def_model()
    model.plot()


def test_getitem_setitem_add():
    # __getitem__
    model1 = prs.read_model(dipfile)
    assert model1[0]["thickn"] == 20000.0
    assert model1[1]["rho"] == 2800.0
    assert model1[2]["vp"] == 7800.0
    assert model1[2]["vs"] == 4480.0
    assert model1[2]["flag"] == 1
    assert model1[1]["dip"] == 20.0
    assert model1[2]["strike"] == 90

    model2 = prs.read_model(anisofile)
    assert model2[1]["flag"] == 0
    assert model2[1]["ani"] == -20
    assert model2[1]["trend"] == 180
    assert model2[1]["plunge"] == 45

    # __eq__
    model22 = prs.read_model(anisofile)
    assert model2 == model22

    # __setitem__
    model2[1, "plunge"] = 10
    assert model2[1]["plunge"] == 10
    assert model2[1, "plunge"] == 10
    assert model2.plunge[1] == 10
    assert model2._plunge[1] == 10
    assert model2.fplunge[1] == 10 * np.pi / 180
    
    # @property.setter
    model2.plunge[1] = 20
    assert model2.plunge[1] == 20

    # TODO: The above did not trigger the setter yet
    model2.update()
    assert model2.fplunge[1] == 20 * np.pi / 180
    assert model2[1, "plunge"] == 20

    # Does anisotropy flag get (de-)activated automatically?
    model1[0, "ani"] = 10
    assert model1[0]["ani"] == 10
    assert model1[0]["flag"] == 0
    model1[0, "ani"] = 0
    assert model1[0]["ani"] == 0
    assert model1[0]["flag"] == 1

    model1.ani[0] = 10
    assert model1.ani[0] == 10

    model1.ani[0] = 0
    assert model1.ani[0] == 0

    model2[2] = {"plunge": 15, "thickn": 10000}
    assert model2[2, "plunge"] == 15
    assert model2[2, "thickn"] == 10000

    with pytest.raises(ValueError):
        # Only can set user attributes
        model2[1, "fplunge"] = 10

    # __add__
    model3 = model1 + model2
    assert model3[0, "thickn"] == model1[0, "thickn"]
    assert model3[3, "thickn"] == model2[0, "thickn"]
    assert model3[5, "thickn"] == model2[2, "thickn"]

    assert model3.thickn[0] == model1.thickn[0]
    assert model3.thickn[3]== model2.thickn[0]
    assert model3.thickn[5] == model2.thickn[2]

    model4 = model1 + {"thickn": 5000, "rho": 3000, "vp": 6000}
    assert model4[3, "thickn"] == 5000

    model5 = model1 + [5000, 3000, 6000]
    assert model5[3, "thickn"] == 5000

    with pytest.raises(TypeError):
        # Only can add another model
        model1 + 5

def test_average_model():
    model = Model([1000, 1000, 0], [1000, 3000, 4000], [1000, 3000, 4000], [1000, 3000, 4000])
    model.average_layers(0, 2)
    assert model[0, "vp"] == 2000
    assert model[0, "vs"] == 2000
    assert model[0, "rho"] == 2000
    assert model[0, "thickn"] == 2000
    assert model.nlay == 2