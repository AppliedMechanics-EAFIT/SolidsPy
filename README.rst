SolidsPy: 2D-Finite Element Analysis with Python
================================================

.. figure:: https://raw.githubusercontent.com/AppliedMechanics-EAFIT/SolidsPy/master/docs/img/wrench.png
   :alt: Wrench under bending.



.. image:: https://img.shields.io/pypi/v/solidspy.svg
   :target: https://pypi.python.org/pypi/continuum_mechanics
   :alt: PyPI download

.. image:: https://readthedocs.org/projects/solidspy/badge/?version=latest
   :target: https://solidspy.readthedocs.io/en/latest/
   :alt: Documentation Status

.. image:: https://img.shields.io/pypi/dm/solidspy
   :target: https://pypistats.org/packages/solidspy
   :alt: Downloads frequency



A simple finite element analysis code for 2D elasticity problems.
The code uses as input simple-to-create text files
defining a model in terms of nodal, element, material and load data.

- Documentation: http://solidspy.readthedocs.io
- GitHub: https://github.com/AppliedMechanics-EAFIT/SolidsPy
- PyPI: https://pypi.org/project/solidspy/
- Free and open source software: `MIT license <http://en.wikipedia.org/wiki/MIT_License>`__



Features
--------

* It is based on an open-source environment.

* It is easy to use.

* The code allows to find displacement, strain and stress solutions for
  arbitrary two-dimensional domains discretized into finite elements and
  subject to point loads.

* The code is organized in independent modules for pre-processing, assembly and
  post-processing allowing the user to easily modify it or add features like
  new elements or analyses pipelines.

* It was created with academic and research purposes.

* It has been used to tech the following courses:

  - Computational Modeling.
  - Introduction to the Finite Element Methods.


Installation
------------

The code is written in Python and it depends on ``numpy``, ``scipy`` and
``sympy``. It has been tested under Windows, Mac, Linux and Android.

To install *SolidsPy* open a terminal and type:

::

    pip install solidspy

To specify through a GUI the folder where the input
files are stored  you will need to install `easygui <http://easygui.readthedocs.org/en/master/>`__.

To easily generate the required SolidsPy text files out of a
`Gmsh <http://gmsh.info/>`__ model you will need
`meshio <https://github.com/nschloe/meshio>`__.

These two can be installed with:

::

    pip install easygui
    pip install meshio


How to run a simple model
-------------------------

For further explanation check the `docs <http://solidspy.readthedocs.io/en/latest/>`__.

Let's suppose that we have a simple model represented by the following
files (see `tutorials/square example <http://solidspy.readthedocs.io/en/latest/tutorials/square_example.html>`__
for further explanation).

.. code:: python

    import numpy as np
    from solidspy.solids_GUI import solids_auto

    ### Define the data
    nodes = np.array([
        [0, 0.00, 0.00],
        [1, 2.00, 0.00],
        [2, 2.00, 2.00],
        [3, 0.00, 2.00],
        [4, 1.00, 0.00],
        [5, 2.00, 1.00],
        [6, 1.00, 2.00],
        [7, 0.00, 1.00],
        [8, 1.00, 1.00]])

    cons = np.array([
        [0, -1],
        [0, -1],
        [0,  0],
        [0,  0],
        [-1, -1],
        [0,  0],
        [0,  0],
        [0,  0],
        [0,  0]])

    elements = np.array([
        [0, 1, 0, 0, 4, 8, 7],
        [1, 1, 0, 4, 1, 5, 8],
        [2, 1, 0, 7, 8, 6, 3],
        [3, 1, 0, 8, 5, 2, 6]])

    mats = np.array([[1.0, 0.3]])

    loads = np.array([
        [2, 0.0, 1.0],
        [3, 0.0, 1.0],
        [6, 0.0, 2.0]])

    data = {"nodes": nodes,
            "cons": cons,
            "elements": elements,
            "mats": mats,
            "loads": loads}

    ### Run the simulation
    disp = solids_auto(data)
    plt.show()


Save it as ``example_solidspy.py`` and run it in your terminal:

.. code:: bash

    python example_solidspy.py


License
-------

This project is licensed under the `MIT
license <http://en.wikipedia.org/wiki/MIT_License>`__. The documents are
licensed under `Creative Commons Attribution
License <http://creativecommons.org/licenses/by/4.0/>`__.

Citation
--------

To cite SolidsPy in publications use

    Juan Gómez, Nicolás Guarín-Zapata (2018). SolidsPy: 2D-Finite
    Element Analysis with Python, <https://github.com/AppliedMechanics-EAFIT/SolidsPy>.

A BibTeX entry for LaTeX users is

.. code:: bibtex

    @software{solidspy,
     title = {SolidsPy: 2D-Finite Element Analysis with Python},
     author = {Gómez, Juan and Guarín-Zapata, Nicolás},
     year = 2018,
     keywords = {Python, Finite elements, Scientific computing, Computational mechanics},
     abstract = {SolidsPy is a simple finite element analysis code for
       2D elasticity problems. The code uses as input simple-to-create text
       files defining a model in terms of nodal, element, material and
       load data.},
     url = {https://github.com/AppliedMechanics-EAFIT/SolidsPy}
    }
