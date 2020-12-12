
.. figure:: ../pyraysum/examples/picture/PyRaysum_logo.png
   :align: center

Documentation
=============

``PyRaysum`` is a Python wrapper around the Fortran software ``Raysum``, originally developed by `Andrew Frederiksen <https://umanitoba.ca/faculties/environment/departments/geo_sciences/research_facilities/AndrewFrederiksen.html>`_ in collaboration with `Michael Bostock <https://www.eoas.ubc.ca/people/michaelbostock>`_. This program generates sets of ray-theoretical seismograms for an incident plane wave (teleseismic approximation) for models consisting of a stack of layers with planar but nonparallel (dipping) interfaces, allowing the possibility of anisotropy in the layers. Incident P and S waves are supported.

``PyRaysum`` is bundled with a trimmed version of the original Fortran software and provides a script to compile and install Raysum, as well as functions to interact with the software and generate plots of models and seismograms. Common computational workflows are covered in the ``Jupyter`` notebooks bundled with this package.

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

   Example 1: A simple example for OBS <https://nbviewer.jupyter.org/github/paudetseis/PyRaysum/blob/main/pyraysum/examples/notebooks/sim_Porter2011.ipynb>

