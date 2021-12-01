
.. figure:: ../pyraysum/examples/picture/PyRaysum_logo.png
   :align: center

Documentation
=============

``PyRaysum`` is a Python wrapper around the Fortran software ``Raysum``, originally developed by `Andrew Frederiksen <https://umanitoba.ca/faculties/environment/departments/geo_sciences/research_facilities/AndrewFrederiksen.html>`_ in collaboration with `Michael Bostock <https://www.eoas.ubc.ca/people/michaelbostock>`_. This program generates sets of ray-theoretical seismograms for an incident plane wave (teleseismic approximation) for models consisting of a stack of layers with planar (but possibly non-parallel and dipping) interfaces, allowing the flexibility of adding seismic anisotropy in the layers. Incident P and S waves are supported. The code include methods to process the synthetic seismograms for receiver function calculation and for plotting. 

``PyRaysum`` is bundled with a trimmed version of the original Fortran software and provides a module to interact with Raysum, as well as functions to post-process the data and generate plots of models and seismograms. Common computational workflows are covered in the ``Jupyter`` notebooks bundled with this package.

.. image:: https://github.com/paudetseis/PyRaysum/workflows/Build/badge.svg
    :target: https://github.com/paudetseis/PyRaysum/actions
.. image:: https://codecov.io/gh/paudetseis/PyRaysum/branch/main/graph/badge.svg
    :target: https://codecov.io/gh/paudetseis/PyRaysum

.. toctree::
   :maxdepth: 1
   :caption: Quick Links

   links

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   init

.. toctree::
   :maxdepth: 1
   :caption: API

   api

.. toctree::
   :maxdepth: 1
   :caption: Jupyter Notebooks

   Example 1: Reproducing Figure 2 in Porter et al. (2011) <https://nbviewer.jupyter.org/github/paudetseis/PyRaysum/blob/main/pyraysum/examples/notebooks/sim_Porter2011.ipynb>

