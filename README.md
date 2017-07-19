# SolidsPy: 2D-Finite Element Analysis with Python


![Wrench under bending.](./docs/img/wrench.png)

This _repo_ contains a simple finite element analysis code for 2D elasticity
problems. The code uses as input data simple-to-create text files containing
nodal, element, material and load data.

The _repo_ contains 3 main folders:

1. `solidspy/` stores the package routines:

  - `solids_GUI.py`: run the program;
  - `preprocesor.py` model input subroutines);
  - `assemutil.py` (assembly subroutines);
  - `femutil.py` (general finite element method subroutines);
  - `uelutil.py` (local matrix subroutines for different elements; and
  - `postprocesor.py` (results handling subroutines)

2. `meshes/` contains input files and meshes in `.msh` format corresponding to
    different examples.

3. `docs/` contains the documentation files.

4. `examples/` contains some examples of the use of the package.

5. `tests/` contains unit testing.

## Authors
- [Juan David Gómez Cataño](http://www.eafit.edu.co/docentes-investigadores/Paginas/juan-gomez.aspx),
    Professor at Universidad EAFIT.
- [Nicolás Guarín-Zapata](https://github.com/nicoguaro), PhD Student at
    Purdue University.

## Installation
The code is written in Python and it depends on `numpy`, `scipy` and `sympy`.

To run the examples with GUI input you will need to install
[`easygui`](http://easygui.readthedocs.org/en/master/). And, you will
need [`meshio`](https://github.com/nschloe/meshio) to automatically read
[Gmsh](http://gmsh.info/) mesh files. These two can be installed with

    pip install easygui
    pip install meshio

To install _SolidsPy_ use

    pip install solidspy

## Run a simple model
After installation, you can run an analysis in 3 easy steps (see [template](./docs/template/README.md)):
- Create the mesh using [Gmsh](http://gmsh.info/).
- Generate the model files (eles.txt, nodes.txt, mater.txt and loads.txt) using
  a python script with the aid of [`meshio`](https://github.com/nschloe/meshio).
- Load the function `solids_GUI` with (and also matplotlib)

    import matplotlib.pyplot as plt
    from solidspy import solids_GUI

  run it

    solids_GUI()

  and plot the contours with

    plt.show()

## License
This project is licensed under the
[MIT license](http://en.wikipedia.org/wiki/MIT_License). The documents are
licensed under
[Creative Commons Attribution License](http://creativecommons.org/licenses/by/4.0/).

Since this project is used to teach Finite Element Methods and Computational
Mechanics, we have included some examples and documents and code snippets
developed by students. Their license is specified in each particular directory.
