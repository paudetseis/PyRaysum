from pyraysum import Control

def test_def_control():
    ctrl = Control()
    assert ctrl.verbose == 0
    ctrl.verbose = "1"
    assert ctrl.verbose == 1
    assert ctrl.parameters[7] == 1