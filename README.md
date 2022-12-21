![](./pyraysum/examples/picture/PyRaysum_logo.png)
## Software for modeling ray-theoretical body-wave propagation

This program generates sets of ray-theoretical seismograms for incident plane waves 
(teleseismic approximation) for seismic velocity models consisting of a stack of layers 
with planar but non-parallel (dipping) interfaces, allowing the possibility of anisotropy 
in the layers. Incident P and S waves are supported.

`PyRaysum` is a Python wrapper around the Fortran software `Raysum`, originally developed by 
[Andrew Frederiksen](https://umanitoba.ca/faculties/environment/departments/geo_sciences/research_facilities/AndrewFrederiksen.html). 
A trimmed down version of the Fortran code is supplied with `PyRaysum`. You can find the 
original version [here](https://home.cc.umanitoba.ca/~frederik/Software/).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6095749.svg)](https://doi.org/10.5281/zenodo.6095749)
[![build](https://github.com/paudetseis/PyRaysum/workflows/Build/badge.svg)](https://github.com/paudetseis/PyRaysum/actions)
[![codecov](https://codecov.io/gh/paudetseis/PyRaysum/branch/main/graph/badge.svg?token=59F1SWLM9Q)](https://codecov.io/gh/paudetseis/PyRaysum)
![GitHub](https://img.shields.io/github/license/paudetseis/pyraysum)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Authors: [Wasja Bloch](https://www.eoas.ubc.ca/people/wasjabloch), [Pascal Audet](https://www.uogeophysics.com/authors/admin/) (Developers and Maintainers of `PyRaysum`) & [Andrew Frederiksen](https://umanitoba.ca/faculties/environment/departments/geo_sciences/research_facilities/AndrewFrederiksen.html) & (Developer of original `Fortran` version)


#### Installation

*PyRaysum* can be installed from PyPI or from source.

To avoid conflicts with other programs, it is recommended to install *PyRaysum* inside a designated `conda` environment (called here `prs`) alongside its dependecies.
```
conda create -n prs python fortran-compiler obspy -c conda-forge
conda activate prs
```

##### Installing from PyPI

```
pip install pyraysum
```

##### Installing from source

The source code of *PyRaysum* can also be downloaded from *GitHub* and installed via `pip`.

```
git clone https://github.com/paudetseis/PyRaysum.git
cd PyRaysum
pip install .
```

#### Getting Started

To compute receiver functions for a range of back-azimuths considering a simple 1-layer 
over a half-space subsurface seismic velocity model, in a *Python* or *iPython* console, execute:
```
from pyraysum import Model, Geometry, Control, run

# Build a 1-layer-over-half-space subsurface model
model = Model(
    thickn=[32000, 0],  # m; half-space thickness is irrelevant
    rho=[2800, 3600],  # kg/m^3
    vp=[6400, 8100],  # m/s
    vs=[3600, 4650],  # m/s
)
model.plot()

# Define back-azimuth range and single horizontal slowness value
geometry = Geometry(baz=range(0, 360, 30), slow=0.07)
geometry.plot()

# Set sampling interval, number of points, alignment and ray-coordinate rotation
control = Control(dt=1./20., npts=800, align="P", rot="PVH")

# Run the simulation
result = run(model, geometry, control, rf=True)

# Filter below 2s period
result.filter("rfs", "lowpass", freq=0.5)

# Plot the results
result.plot("rfs")
```

#### Documentation
The complete API documentation, scripts and tutorials are described at https://paudetseis.github.io/PyRaysum/

#### Citing

If you use `PyRaysum` in your work, please cite the following references:

- Frederiksen, A.W., and Bostock, M.G. (1999) Modelling teleseismic waves in dipping anisotropic structures. Geophysical Journal International 141: 401-412. https://doi.org/10.1046/j.1365-246x.2000.00090.x

- Bloch, W., and Audet, P. (in press). PyRaysum: Software for modeling ray-theoretical plane body-wave propagation in dipping anisotropic media. Seismica.

- Audet, P., and Bloch, W. (2022). PyRaysum: Software for modeling ray-theoretical body-wave propagation. Zenodo. https://doi.org/10.5281/zenodo.6095749

#### Contributing

All constructive contributions are welcome, e.g. bug reports, discussions or suggestions for new features. You can either [open an issue on GitHub](https://github.com/paudetseis/PyRaysum/issues) or make a pull request with your proposed changes. Before making a pull request, check if there is a corresponding issue opened and reference it in the pull request. If there isn't one, it is recommended to open one with your rationale for the change. New functionality or significant changes to the code that alter its behavior should come with corresponding tests and documentation. If you are new to contributing, you can open a work-in-progress pull request and have it iteratively reviewed. 

Other examples of contributions include notebooks that describe published examples of `PyRaysum` usage and processing. Suggestions for improvements (speed, accuracy, plotting, etc.) are also welcome.

