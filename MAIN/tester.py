# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 18:52:09 2015

@author: eafit
"""

from __future__ import division
import numpy as np
from sympy import *
import hello
from sympy import init_printing
init_printing()
#
print hello.__doc__
print hello.foo.__doc__
print hello.suma.__doc__
hello.foo(4)
c=hello.suma(4.0,5.0)
