# SolidsPy

The aim of this project is to find the displacement, strain and stress
solution for an arbitrary two-dimensional domain discretized into finite
elements and subjected to point loads. It has been  created for
academic purposes and it is part of the teaching material developed for
the courses IC0602 Introduction to the Finite Element Methods and
IC0285 Computational Modeling at Universidad EAFIT. The code is
written in Python 2 dialect (although, most of it works as well
in Python 3) and it is organized in independent modules for
pre-processing, assembly and post-processing allowing the user to
easily modify it or add features like new elements.

## Contents

- [2×2 square with axial load](square_example.md): describes the use of SolidsPy, through a simple
example corresponding to a square plate under point loads.
- [Creation of a simple geometry using Gmsh](geometry_gmsh/README.md) show the user how to create a model using
the third-party software [Gmsh](http://gmsh.info/).
