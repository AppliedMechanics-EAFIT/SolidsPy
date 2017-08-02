# SolidsPy: 2D-Finite Element Analysis with Python


![Wrench under bending.](./docs/img/wrench.png)

This _repo_ contains a simple finite element analysis code for 2D elasticity
problems. The code uses as input simple-to-create text files defining a model in terms of
nodal, element, material and load data.

The _repo_ contains 4 main folders:

1. `solidspy/` Which stores the source code in the following modules:

    - `solids_GUI.py`: The main program;
    - `preprocesor.py` Pre-processing subrotuines including Gmsh convertion functions using meshio;
    - `assemutil.py` Assembly of elemental stiffnesss matrices ;
    - `femutil.py` Shape functions, iots derivatives and general finite element method subroutines;
    - `uelutil.py` Elemental or local matrix subroutines for different elements; and
    - `postprocesor.py` Several results handling subroutines

2. `meshes/` Complete models including its gmsh representation and a Python script to produce the required
    (nodes, elements, materials and load) text files ready for input.

3. `docs/` Contains the documentation files like easy-to-follow tutorials
     showing how to define a SolidsPy model in trms of text files and model
     creation with gmsh.

4. `examples/` Specific applications using SolidsPy functions to conduct analysis.

5. `tests/` contains unit testing.

## Authors
- [Juan Gomez](http://www.eafit.edu.co/docentes-investigadores/Paginas/juan-gomez.aspx),
    Professor at Universidad EAFIT.
- [Nicolás Guarín-Zapata](https://github.com/nicoguaro), PhD Student at
    Purdue University.

## Installation
The code is written in Python and it depends on `numpy`, `scipy` and `sympy`.

To install _SolidsPy_ use:

    pip install solidspy

To run the examples with GUI input you will need to install
[`easygui`](http://easygui.readthedocs.org/en/master/).

To easily generate the required SolidsPy text files out of a [Gmsh](http://gmsh.info/) model
you will need [`meshio`](https://github.com/nschloe/meshio).These two can be installed with:

    pip install easygui
    pip install meshio

## Run a simple model
After installation, you can run an analysis in 3 easy steps (see [template](./docs/template/README.md)):
- Create the mesh using [Gmsh](http://gmsh.info/).
- Generate the model files (eles.txt, nodes.txt, mater.txt and loads.txt) using
  a python script with the aid of [`meshio`](https://github.com/nschloe/meshio).
- Run it in Python

```python
import matplotlib.pyplot as plt  # load matplotlib
from solidspy import solids_GUI  # import our package
solids_GUI()  # run the Finite Element Analysis
plt.show()    # plot contours
```

## License
This project is licensed under the
[MIT license](http://en.wikipedia.org/wiki/MIT_License). The documents are
licensed under
[Creative Commons Attribution License](http://creativecommons.org/licenses/by/4.0/).

Since this project is used to teach Finite Element Methods and Computational
Mechanics, we have included some examples and documents and code snippets
developed by students. Their license is specified in each particular directory.
