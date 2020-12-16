![](./pyraysum/examples/picture/PyRaysum_logo.png)
## Software for modeling ray-theoretical body-wave propagation

This program generates sets of ray-theoretical seismograms for an
incident plane wave (teleseismic approximation) for models consisting
of a stack of layers with planar but nonparallel (dipping) interfaces,
allowing the possibility of anisotropy in the layers. Incident P and S
waves are supported.

`PyRaysum` is a Python wrapper around the Fortran software `Raysum`, originally developed by [Andrew Frederiksen](https://umanitoba.ca/faculties/environment/departments/geo_sciences/research_facilities/AndrewFrederiksen.html). Although a trimmed down version of the Fortran code is supplied with `PyRaysum`, you can find the original version [here](https://home.cc.umanitoba.ca/~frederik/Software/).

![build](https://github.com/paudetseis/PyRaysum/workflows/Build/badge.svg)
[![codecov](https://codecov.io/gh/paudetseis/PyRaysum/branch/main/graph/badge.svg?token=59F1SWLM9Q)](https://codecov.io/gh/paudetseis/PyRaysum)

Installation, API documentation, scripts and tutorials are described at https://paudetseis.github.io/PyRaysum/

Authors: [`Pascal Audet`](https://www.uogeophysics.com/authors/admin/) (Developer and Maintainer of `PyRaysum`) & [`Andrew Frederiksen`](https://umanitoba.ca/faculties/environment/departments/geo_sciences/research_facilities/AndrewFrederiksen.html) & (Developer of original `Fortran` version)

#### Contributing

All constructive contributions are welcome, e.g. bug reports, discussions or suggestions for new features. You can either [open an issue on GitHub](https://github.com/paudetseis/PyRaysum/issues) or make a pull request with your proposed changes. Before making a pull request, check if there is a corresponding issue opened and reference it in the pull request. If there isn't one, it is recommended to open one with your rationale for the change. New functionality or significant changes to the code that alter its behavior should come with corresponding tests and documentation. If you are new to contributing, you can open a work-in-progress pull request and have it iteratively reviewed. 

Examples of contributions include notebooks that describe published examples of `PyRaysum` usage and processing. Suggestions for improvements (speed, accuracy, plotting, etc.) are also welcome.

#### Citing

If you use `PyRaysum` in your work, please cite the Zenodo DOI (not yet available) and the following paper:

- Frederiksen, A.W., and Bostock, M.G. (1999) Modelling teleseismic waves in dipping anisotropic structures. Geophysical Journal International 141: 401-412.

<!-- ### Dependencies

`PyRaysum` requires a modern Fortran compiler (e.g., gfortran, ifort). In addition, the following packages are required:

- [pandas](https://pandas.pydata.org)
- [obspy](https://docs.obspy.org)

### Installing

Documentation for `PyRaysum` is in development and this Readme file provides minimal guidance on how to install and run the code. To install `PyRaysum` using `conda`, you can first create a `conda` environment and install the dependencies at the same time:

```
conda create -n prs python=3.8 obspy pandas fortran-compiler -c conda-forge
conda activate prs
```

Then clone and install `PyRaysum` using `pip`:

```
git clone https://github.com/paudetseis/PyRaysum.git
cd PyRaysum
pip install .
```

After `PyRaysum` is installed, you need to compile and install the Raysum binaries by running the provided `install-raysum` script. Check out the helper page for the available options:

```
install-raysum -h
```

Running the script without arguments will attempt to install Raysum to default paths and fortran compiler (these are contained in the $PATH and $FC environment variables). If you installed `PyRaysum` using `conda` (and the `fortran-compiler` package), you can safely ignore the options and type `install-raysum` in the terminal. To call the script with options (which take precedence over any default `conda` environment variables), you would specify:

```
install-raysum --path=/path/to/binaries/bin --fcompiler=gfortran
```

Note that the `--path` option requires the full, absolute path. If you specify a path in one of the root directories, you may have to run `install-raysum` with super-user privileges, for example:

```
sudo install-raysym --path=/usr/bin
```

### Notebooks

There is currently one notebook bundled with `PyRaysum`, which reproduces some published examples in [Porter et al., 2011](https://doi.org/10.1130/L126.1). To download it locally and run it, first install `jupyter`:

```
conda install jupyter
```

Then open a `Python` window and type:

```
from pyraysum import doc
doc.install_doc(path='Notebooks')
```

which will install the notebook locally in the folder `'Notebooks'`. Then you simply open a `jupyter notebook` and navigate to the `Notebook` folder and open the corresponding file. 
 -->
  

