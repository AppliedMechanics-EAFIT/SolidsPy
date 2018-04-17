SolidsPy: 2D-Finite Element Analysis with Python
================================================

.. figure:: https://raw.githubusercontent.com/AppliedMechanics-EAFIT/SolidsPy/master/docs/img/wrench.png
   :alt: Wrench under bending.

   Wrench under bending.

A simple finite element analysis code for 2D elasticity problems.
The code uses as input simple-to-create text files 
defining a model in terms of nodal, element, material and load data.

It has been created for academic purposes and it is part of the
teaching material developed for the courses IC0602 Introduction to
the Finite Element Methods and IC0285 Computational Modeling at
Universidad EAFIT.

Features
--------

The code allows the user to find the displacement, strain and stress
solution for an arbitrary two-dimensional domain discretized into finite
elements and subjected to point loads. 

The code is organized in independent modules for pre-processing, assembly
and post-processing allowing the user to easily modify it or add features
like new elements or analyses.


Installation
------------

The code is written in Python and it depends on ``numpy``, ``scipy`` and
``sympy``.

To install *SolidsPy* open a terminal and type:

::

    pip install solidspy

To run the examples with specification of the folder stoing the input
files through a GUI you will need to install 
`easygui <http://easygui.readthedocs.org/en/master/>`__.

To easily generate the required SolidsPy text files out of a
`Gmsh <http://gmsh.info/>`__ model you will need
`meshio <https://github.com/nschloe/meshio>`__.

These two can be installed with:

::

    pip install easygui
    pip install meshio

How to run a simple model
-------------------------

After installation, you can run an analysis in 3 easy steps (see
`template <./docs/template/README.md>`__): - Create the model (i.e.,
geometry and mesh) using `Gmsh <http://gmsh.info/>`__. Several meshes
are available in the repo
`SOLIDSPy-meshes <https://github.com/AppliedMechanics-EAFIT/SolidsPy-meshes>`__
- Generate the text files (eles.txt, nodes.txt, mater.txt and loads.txt)
required by *SolidsPy* using a python script based on
`meshio <https://github.com/nschloe/meshio>`__. - Run it in Python
as follows:

.. code:: python

    import matplotlib.pyplot as plt  # load matplotlib
    from solidspy import solids_GUI  # import our package
    disp = solids_GUI()  # run the Finite Element Analysis
    plt.show()    # plot contours

This would not work properly in Anaconda for Mac OS. In that case is
suggested to use an IPython console to run the example.

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

.. code-block::

    @software{solidspy,
     title = {SolidsPy: 2D-Finite Element Analysis with Python},
     author = {Gómez, Juan and Guarín-Zapata, Nicolás},
     year = 2017,
     keywords = {Python, Computer algebra system, Symbolics},
     abstract = {SolidsPy is a simple finite element analysis code for
       2D elasticity problems. The code uses as input simple-to-create text
       files defining a model in terms of nodal, element, material and
       load data.},
     url = {https://github.com/AppliedMechanics-EAFIT/SolidsPy}
    }
