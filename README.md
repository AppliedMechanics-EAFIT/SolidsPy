# SolidsPy: 2D/3D-Finite Element Analysis with Python

![Wrench under bending](https://raw.githubusercontent.com/AppliedMechanics-EAFIT/SolidsPy/master/docs/img/wrench.png)

[![PyPI version](https://img.shields.io/pypi/v/solidspy.svg)](https://pypi.python.org/pypi/continuum_mechanics)
[![Documentation Status](https://readthedocs.org/projects/solidspy/badge/?version=latest)](https://solidspy.readthedocs.io/en/latest/)
[![Downloads frequency](https://img.shields.io/pypi/dm/solidspy)](https://pypistats.org/packages/solidspy)
[![DOI](https://zenodo.org/badge/48294591.svg)](https://zenodo.org/badge/latestdoi/48294591)

A simple finite element analysis code for 2D/3D elasticity problems, written in Python. The code uses easily created text files defining a model in terms of nodal, element, material, and load data.

- **Documentation**: http://solidspy.readthedocs.io
- **GitHub**: https://github.com/AppliedMechanics-EAFIT/SolidsPy
- **PyPI**: https://pypi.org/project/solidspy/
- **License**: [MIT license](http://en.wikipedia.org/wiki/MIT_License)
- **Year**: 2025

## Features

* **Open-Source Environment**: Entirely written in Python, leveraging popular libraries like NumPy and SciPy.

* **Easy to Use**: Simple text input files; minimal overhead to set up new simulations.

* **2D/3D Elasticity**: Solve displacement, strain, and stress for arbitrary 2D/3D domains using finite elements.

* **Modular Code**: Independent modules for pre-processing, assembly, and post-processing. Extend or modify as needed.

* **Academic & Research**: Ideal for teaching courses such as:
  - Computational Modeling
  - Introduction to the Finite Element Methods
  
  and for rapid prototyping of new FEM elements.

* **Wide Compatibility**: Tested under Windows, macOS, Linux, and even Android (via Termux).

## Installation

SolidsPy runs on Python 3.11+ and depends on numpy and scipy.

Install via PyPI:
```bash
pip install solidspy
```

For a GUI file selector, install:
```bash
pip install easygui
```

To convert Gmsh models ([gmsh.info](http://gmsh.info/)) to SolidsPy input files, install:
```bash
pip install meshio
```

Install via Conda:
```bash
conda install solidspy
```

## How to run a simple model

For further explanation check the [docs](http://solidspy.readthedocs.io/en/latest/).

Let's suppose that we have a simple model represented by the following files (see [tutorials/square example](http://solidspy.readthedocs.io/en/latest/tutorials/square_example.html) for further explanation).

- nodes.txt
```
0  0.00  0.00   0  -1
1  2.00  0.00   0  -1
2  2.00  2.00   0   0
3  0.00  2.00   0   0
4  1.00  0.00  -1  -1
5  2.00  1.00   0   0
6  1.00  2.00   0   0
7  0.00  1.00   0   0
8  1.00  1.00   0   0
```

- eles.txt
```
0   1   0   0   4   8   7
1   1   0   4   1   5   8
2   1   0   7   8   6   3
3   1   0   8   5   2   6
```

- mater.txt
```
1.0  0.3
```

- loads.txt
```
3  0.0  1.0
6  0.0  2.0
2  0.0  1.0
```

Run it in Python as follows:
```python
import matplotlib.pyplot as plt  # load matplotlib
from solidspy import solids_GUI  # import our package
disp = solids_GUI()  # run the Finite Element Analysis
plt.show()    # plot contours
```

For Mac users it is suggested to use an IPython console to run the example.

## License

This project is licensed under the [MIT license](http://en.wikipedia.org/wiki/MIT_License). The documents are licensed under [Creative Commons Attribution License](http://creativecommons.org/licenses/by/4.0/).

## Citation

To cite SolidsPy in publications use:

> Nicolás Guarín-Zapata, Juan Gomez (2025). SolidsPy: Version 2.0.0 (Version v2.0.0). Zenodo. http://doi.org/10.5281/zenodo.4029270

A BibTeX entry for LaTeX users is:

```bibtex
@software{solidspy,
 title = {SolidsPy: 2D/3D-Finite Element Analysis with Python},
 version = {2.0.0},
 author = {Guarín-Zapata, Nicolás and Gómez, Juan},
 year = 2020,
 keywords = {Python, Finite elements, Scientific computing, Computational mechanics},
 abstract = {SolidsPy is a simple finite element analysis code for
   2D/3D elasticity problems. The code uses as input simple-to-create text
   files defining a model in terms of nodal, element, material and
   load data.},
 url = {https://github.com/AppliedMechanics-EAFIT/SolidsPy},
 doi = {http://doi.org/10.5281/zenodo.4029270}
}
```