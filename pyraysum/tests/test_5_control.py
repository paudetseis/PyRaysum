from pyraysum import Control

def test_def_control():
    ctrl = Control(align="SV", rot="PVH")

    assert ctrl.verbose == 0

    ctrl.wvtype = "SH"
    assert ctrl._iphase == 3
    assert ctrl.parameters[0] == 3

    ctrl.mults = 1
    assert ctrl.mults == 1
    assert ctrl.parameters[1] == 1

    ctrl.npts = 99
    assert ctrl.npts == 99
    assert ctrl.parameters[2] == 99
    
    ctrl.dt = 0.33
    assert ctrl.dt == 0.33
    assert ctrl.parameters[3] == 0.33

    ctrl.align = 2
    assert ctrl.align == 2
    assert ctrl.parameters[4] == 2
    
    ctrl.align = "none"
    assert ctrl.align == 0
    assert ctrl.parameters[4] == 0
    
    ctrl.align = "P"
    assert ctrl.align == 1
    assert ctrl.parameters[4] == 1
    
    ctrl.align = "SV"
    assert ctrl.align == 2
    assert ctrl.parameters[4] == 2
    
    ctrl.align = "SH"
    assert ctrl.align == 3
    assert ctrl.parameters[4] == 3
    
    ctrl.shift = 5
    assert ctrl.shift == 5
    assert ctrl.parameters[5] == 5
    
    ctrl.rot = 2
    assert ctrl.rot == 2
    assert ctrl.parameters[6] == 2
    
    ctrl.rot = "ZNE"
    assert ctrl.rot == 0
    assert ctrl.parameters[6] == 0
    
    ctrl.rot = "RTZ"
    assert ctrl.rot == 1
    assert ctrl.parameters[6] == 1
    
    ctrl.verbose = "1"
    assert ctrl.verbose == 1
    assert ctrl.parameters[7] == 1