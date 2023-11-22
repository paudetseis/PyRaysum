from pkg_resources import resource_filename
from pyraysum import prs, Model
import pytest
import tempfile
from os import remove
import numpy as np

dipfile = resource_filename("pyraysum", "examples/data/model_Porter2011_dip.txt")
anisofile = resource_filename("pyraysum", "examples/data/model_Porter2011_aniso.txt")

MODEL = Model(
    [1000, 1000, 0], [1000, 3000, 4000], [1000, 3000, 4000], [1000, 3000, 4000]
)


def test_read_model_dip():
    model = prs.read_model(dipfile)
    assert isinstance(model, Model)
    return model


def test_read_model_aniso():
    model = prs.read_model(anisofile)
    assert isinstance(model, Model)
    return model


def test_write_read_model_v1_2():
    """Write and read a legavy raysum model"""
    model = prs.read_model(anisofile)
    tmpf = tempfile.NamedTemporaryFile(delete=False)
    model.write(tmpf.name, version="raysum")
    model = prs.read_model(tmpf.name, version="raysum")
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
    assert model1[2, "strike"] == 90
    assert model1["strike", 2] == 90

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

    # Set it the other way round
    model1["ani", 0] = 10
    assert model1[0]["ani"] == 10
    assert model1[0]["flag"] == 0
    model1["ani", 0] = 0
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

    with pytest.raises(ValueError):
        # Tuples must be valid
        model2[1, 1] = 10
        model2["thickn", "ani"] = 10

    # __add__
    model3 = model1 + model2
    assert model3[0, "thickn"] == model1[0, "thickn"]
    assert model3[3, "thickn"] == model2[0, "thickn"]
    assert model3[5, "thickn"] == model2[2, "thickn"]

    assert model3.thickn[0] == model1.thickn[0]
    assert model3.thickn[3] == model2.thickn[0]
    assert model3.thickn[5] == model2.thickn[2]

    model4 = model1 + {"thickn": 5000, "rho": 3000, "vp": 6000}
    assert model4[3, "thickn"] == 5000

    model5 = model1 + [5000, 3000, 6000]
    assert model5[3, "thickn"] == 5000

    with pytest.raises(TypeError):
        # Only can add another model
        model1 + 5


def test_average_model():
    mod = MODEL.copy()
    mod.average_layers(0, 2)
    assert mod[0, "vp"] == 2000
    assert mod[0, "vs"] == 2000
    assert mod[0, "rho"] == 2000
    assert mod[0, "thickn"] == 2000
    assert mod.nlay == 2


def assert_phaselist(phl1, phl2):
    """Assert that all and only all elements in phl1 are present in phl2"""
    assert len(phl1) == len(phl2)

    for ph1, ph2 in zip(sorted(phl1), sorted(phl2)):
        assert ph1 == ph2


def test_get_conversions():
    # Get phases from original fortran routine.
    # Convert to descriptors using frs.read_arrivals

    # Direct phase
    pph = ["2P1P0P"]
    phl = MODEL.get_direct_wave()
    assert_phaselist(phl, pph)

    # Conversions at interface 1
    phl1 = ["2P1P0S", "2P1P0T"]
    phl = MODEL.get_conversions(interfaces=1)
    assert_phaselist(phl, pph + phl1)

    # Conversions at interface 2
    phl2 = ["2P1S0S", "2P1T0T"]

    # Fast <-> Slow conversions
    phl3 = ["2P1S0T", "2P1T0S"]

    phl = MODEL.get_conversions(interfaces=2)
    assert_phaselist(phl, pph + phl2 + phl3)

    phl = MODEL.get_conversions(interfaces=2, direct=False)
    assert_phaselist(phl, phl2 + phl3)

    # Conversions at interface 2 get converted again at interface 1:
    phl = MODEL.get_conversions(interfaces=1, isotropic=False)
    assert_phaselist(phl, pph + phl1 + phl3)

    phl = MODEL.get_conversions()
    assert_phaselist(phl, pph + phl1 + phl2 + phl3)

    phl = MODEL.get_conversions(direct=False)
    assert_phaselist(phl, phl1 + phl2 + phl3)


def test_get_reflections():
    pph = ["2P1P0P"]
    sph = ["2S1S0S"]
    tph = ["2T1T0T"]
    Pp1P = ["2P1P0P0p0P"]
    Pp1S = ["2P1P0P0p0S", "2P1P0P0p0T"]
    Ps1P = ["2P1P0P0s0P", "2P1P0P0t0P"]
    Ps1S = ["2P1P0P0s0S", "2P1P0P0s0T", "2P1P0P0t0S", "2P1P0P0t0T"]
    Sp1S = [
        "2S1S0S0p0S",
        "2S1S0S0p0T",
        "2S1S0T0p0S",
        "2S1S0T0p0T",
        "2S1T0S0p0S",
        "2S1T0S0p0T",
        "2S1T0T0p0S",
        "2S1T0T0p0T",
    ]
    Ss1P = [
        "2S1S0S0s0P",
        "2S1S0S0t0P",
        "2S1S0T0s0P",
        "2S1S0T0t0P",
        "2S1T0S0s0P",
        "2S1T0S0t0P",
        "2S1T0T0s0P",
        "2S1T0T0t0P",
    ]
    Ts1P = [
        "2T1S0S0s0P",
        "2T1S0S0t0P",
        "2T1S0T0s0P",
        "2T1S0T0t0P",
        "2T1T0S0s0P",
        "2T1T0S0t0P",
        "2T1T0T0s0P",
        "2T1T0T0t0P",
    ]

    phl = MODEL.get_reflections(1, modes=["pP"])
    assert_phaselist(phl, pph + Pp1P)

    phl = MODEL.get_reflections(1, modes=["pS"])
    assert_phaselist(phl, pph + Pp1S)

    # No converted conversions
    phl = MODEL.get_reflections(1, modes=["sP"], direct=False)
    assert_phaselist(phl, [])

    phl = MODEL.get_reflections(1, modes=["sP"], equivalent=True)
    assert_phaselist(phl, pph + Ps1P)

    phl = MODEL.get_reflections(1, wvtype="SV", modes=["sP"], direct=False)
    assert_phaselist(phl, Ss1P)

    phl = MODEL.get_reflections(1, modes=["sS"], direct=False)
    assert_phaselist(phl, Ps1S)

    phl = MODEL.get_reflections(1, modes=["sP", "pS"])
    assert_phaselist(phl, pph + Pp1S)

    phl = MODEL.get_reflections(1, wvtype="SV", modes=["sP", "pS"])
    assert_phaselist(phl, sph + Ss1P)

    phl = MODEL.get_reflections(1, wvtype="SH", modes=["sP", "pS"])
    assert_phaselist(phl, tph + Ts1P)

    phl = MODEL.get_reflections(1, wvtype="S", modes=["sP", "pS"], direct=False, equivalent=True)
    assert_phaselist(phl, Ss1P + Sp1S)

    phl = MODEL.get_reflections(1)
    assert_phaselist(phl, pph + Pp1P + Pp1S + Ps1S)

    Pp2P = ["2P1P0P0p1p1P0P"]
    Pp2S = ["2P1P0P0p1p1S0S", "2P1P0P0p1p1S0T", "2P1P0P0p1p1T0S", "2P1P0P0p1p1T0T"]
    Ps2S = [
        "2P1P0P0s1s1S0S",
        "2P1P0P0s1s1S0T",
        "2P1P0P0s1s1T0S",
        "2P1P0P0s1s1T0T",
        "2P1P0P0s1t1S0S",
        "2P1P0P0s1t1S0T",
        "2P1P0P0s1t1T0S",
        "2P1P0P0s1t1T0T",
        "2P1P0P0t1s1S0S",
        "2P1P0P0t1s1S0T",
        "2P1P0P0t1s1T0S",
        "2P1P0P0t1s1T0T",
        "2P1P0P0t1t1S0S",
        "2P1P0P0t1t1S0T",
        "2P1P0P0t1t1T0S",
        "2P1P0P0t1t1T0T",
    ]

    phl = MODEL.get_reflections(2, ["pP"])
    assert_phaselist(phl, pph + Pp2P)

    phl = MODEL.get_reflections(2, ["pS"])
    assert_phaselist(phl, pph + Pp2S)

    phl = MODEL.get_reflections(2, ["sS"])
    assert_phaselist(phl, pph + Ps2S)

    phl = MODEL.get_reflections(2, ["sP"])
    assert_phaselist(phl, pph)

    phl = MODEL.get_reflections(2, ["pP", "sS"])
    assert_phaselist(phl, pph + Pp2P + Ps2S)

    phl = MODEL.get_reflections(modes=["pP", "sS"])
    assert_phaselist(phl, pph + Pp2P + Ps2S + Pp1P + Ps1S)

    phl = MODEL.get_reflections()
    assert_phaselist(phl, pph + Pp2P + Pp1P + Pp1S + Pp2S + Ps2S + Ps1S)
