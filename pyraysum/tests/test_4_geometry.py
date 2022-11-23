from pyraysum import prs, Geometry
import tempfile
from os import remove


def test_def_geometry():
    geom = Geometry(90, [0.04, 0.08], 10, 20)
    return geom

def test_write_read_geometry():
    geom1 = test_def_geometry()
    tmpf = tempfile.NamedTemporaryFile(delete=False)
    geom1.write(tmpf.name)
    geom2 = prs.read_geometry(tmpf.name)
    remove(tmpf.name)
    assert geom1 == geom2
    assert isinstance(geom1, Geometry)
    assert isinstance(geom2, Geometry)

def test_plot_geometry():
    geom = test_def_geometry()
    geom.plot()

def test_getitem_setitem_add():
    # __getitem__
    geom = test_def_geometry()
    assert geom[0] == (90, 0.04, 10, 20)
    assert geom[-1] == (90, 0.08, 10, 20)
    assert geom[0, "baz"] == 90
    assert geom[0, "slow"] == 0.04
    assert geom[0, "dn"] == 10
    assert geom[0, "de"] == 20
    assert geom["baz"] == [90, 90]
    assert geom["slow"] == [0.04, 0.08]

    # __setitem__
    geom[1] = (270, 0.06, 5, 10)
    assert geom[1] == (270, 0.06, 5, 10)

    # __add__
    geom += (360, 0.06, 15, 15)
    assert geom[2] == (360, 0.06, 15, 15)
    assert len(geom) == 3