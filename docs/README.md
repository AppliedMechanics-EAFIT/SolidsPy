# SolidsPy

Here we show through simple examples how to conduct an analysis using SolidsPy.
First we describe the structure of the text files for the problem of small 2x2 square
plate under axial loading.

In the second part of the documents we describe the creation of a SolidsPy model
with the aid of [Gmsh](http://gmsh.info/). This is necessary when conducting analysis in
large and complex geometries. In this case the Gmsh files need to be converted into text
files using subroutines based upon meshio [meshio](https://pypi.python.org/pypi/meshio).

## Contents

- [2×2 square with axial load](square_example.md): describes the use of SolidsPy, through a simple
example corresponding to a square plate under point loads.
- [Creation of a simple geometry using Gmsh](geometry_gmsh/README.md) show the user how to create a model using
the third-party software [Gmsh](http://gmsh.info/).
