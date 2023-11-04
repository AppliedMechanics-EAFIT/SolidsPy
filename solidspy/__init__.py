# -*- coding: utf-8 -*-
from __future__ import absolute_import
from solidspy.solids_GUI import solids_GUI

__all__ = ["assemutil",
           "femutil",
           "gaussutil",
           "postprocesor",
           "preprocesor",
           "solutil",
           "uelutil",
           "solids_GUI"]

__version__ = "1.1.0"


__citation__ = """@software{solidspy,
 title = {SolidsPy: 2D-Finite Element Analysis with Python},
 author = {Guarín-Zapata, Nicolás and Gómez, Juan},
 year = 2023,
 keywords = {Python, Finite elements, Scientific computing, Computational mechanics},
 abstract = {SolidsPy is a simple finite element analysis code for
   2D elasticity problems. The code uses as input simple-to-create text
   files defining a model in terms of nodal, element, material and
   load data.},
 url = {https://github.com/AppliedMechanics-EAFIT/SolidsPy}
}"""
