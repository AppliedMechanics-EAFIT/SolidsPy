# 2D-Finite Element Analysis with Python

This _repo_ contains a simple finite element analysis code for 2D elasticity problems. The code uses as input data simple-to-create text files containing nodal, element, material and load data.

The _repo_ contains 3 folders:

1. `MAIN/` stores the python scripts divided in

  - `SOLIDS_ISO.py` (the main program),
  - `PREPROCESOR.py` (model input subroutines),
  - `ASSEMUTILS.py` (assembly subroutines), FEMUTIL.py(general finite element method subroutines),
  - `UELUTIL.py` (local matrix subroutines for different elements),
  - `POSTPROCESOR.py` (results handling subroutines);

2. `MESHES/` contains input files and meshes in `.msh` format corresponding to different examples and

3. `MESHUTILS/` stores programs to conduct pre-processing and post-processing from and to GMESH. Since the code has been created for academic purposes it gives as main results the displacement vector corresponding to the nodal points of the domain. These vectors are used, together with the Python function griddata to interpolate and plot the displacement solution all over the domain. Post-processing to compute strain and stress values is left for students projects.


## Authors
- [Juan David Gómez Cataño](http://www.eafit.edu.co/docentes-investigadores/Paginas/juan-gomez.aspx), Professor at Universidad EAFIT.
- [Juan Carlos Vergara](https://github.com/jvergar2), PhD Student at Universidad EAFIT.
- [Nicolás Guarín-Zapata](https://github.com/nicoguaro), PhD Student at Purdue University.

## Instructionsto 
The code is written in Python 2 dialect (we believe that it will work in Python 3 but we have not tested yet) and it depends on `numpy`, `scipy` and `sympy`. To use it clone the repo with 

   git clone https://github.com/jgomezc1/FEM_PYTHON.git
   
uncompress the zip folder an run the main file in the Python console of your preference.


To run the files in `MESHUTILS` you will need to have `gfortran` and run the following `make` instructions

    make --file="CONTOUR Makefile"
    make --file="MESHER Makefile"


## License
This project is licensed under the [MIT license](http://en.wikipedia.org/wiki/MIT_License).
