![](./pyraysum/examples/picture/PyRaysum_logo.png)
## Software for modeling ray-theoretical body-wave propagation

This program generates sets of ray-theoretical seismograms for an
incident plane wave (teleseismic approximation) for models consisting
of a stack of layers with planar but nonparallel (dipping) interfaces,
allowing the possibility of anisotropy in the layers. Incident P and S
waves are supported.

`PyRaysum` is a Python wrapper around the Fortran software `Raysum`, originally developed by [Andrew Frederiksen](https://umanitoba.ca/faculties/environment/departments/geo_sciences/research_facilities/AndrewFrederiksen.html). Although a trimmed down version of the Fortran code is supplied with `PyRaysum`, you can find the original version [here](https://home.cc.umanitoba.ca/~frederik/Software/).

[![build](https://github.com/paudetseis/PyRaysum/workflows/Build/badge.svg)](https://github.com/paudetseis/PyRaysum/actions)
[![codecov](https://codecov.io/gh/paudetseis/PyRaysum/branch/main/graph/badge.svg?token=59F1SWLM9Q)](https://codecov.io/gh/paudetseis/PyRaysum)
![GitHub](https://img.shields.io/github/license/paudetseis/pyraysum)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6095749.svg)](https://doi.org/10.5281/zenodo.6095749)


Installation, API documentation, scripts and tutorials are described at https://paudetseis.github.io/PyRaysum/

Authors: [`Pascal Audet`](https://www.uogeophysics.com/authors/admin/), [`Wasja Bloch`](https://www.eoas.ubc.ca/people/wasjabloch) (Developers and Maintainers of `PyRaysum`) & [`Andrew Frederiksen`](https://umanitoba.ca/faculties/environment/departments/geo_sciences/research_facilities/AndrewFrederiksen.html) & (Developer of original `Fortran` version)

#### Contributing

All constructive contributions are welcome, e.g. bug reports, discussions or suggestions for new features. You can either [open an issue on GitHub](https://github.com/paudetseis/PyRaysum/issues) or make a pull request with your proposed changes. Before making a pull request, check if there is a corresponding issue opened and reference it in the pull request. If there isn't one, it is recommended to open one with your rationale for the change. New functionality or significant changes to the code that alter its behavior should come with corresponding tests and documentation. If you are new to contributing, you can open a work-in-progress pull request and have it iteratively reviewed. 

Examples of contributions include notebooks that describe published examples of `PyRaysum` usage and processing. Suggestions for improvements (speed, accuracy, plotting, etc.) are also welcome.

#### Citing

If you use `PyRaysum` in your work, please cite the [`Zenodo DOI`](https://doi.org/10.5281/zenodo.6095749) and the following paper:

- Frederiksen, A.W., and Bostock, M.G. (1999) Modelling teleseismic waves in dipping anisotropic structures. Geophysical Journal International 141: 401-412.
