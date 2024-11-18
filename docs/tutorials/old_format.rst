Using the old format format for SolidsPy
========================================

:Author: Nicolás Guarín-Zapata
:Date: September 2019

Until version 1.2.8 the data for a simulation was stored in 4 files

- ``nodes.txt``

- ``eles.txt``

- ``mater.txt``

- ``loads.txt``

As shown below.

- ``nodes.txt``

::

    0  0.00  0.00   0  -1
    1  2.00  0.00   0  -1
    2  2.00  2.00   0   0
    3  0.00  2.00   0   0
    4  1.00  0.00  -1  -1
    5  2.00  1.00   0   0
    6  1.00  2.00   0   0
    7  0.00  1.00   0   0
    8  1.00  1.00   0   0

- ``eles.txt``

::

    0   1   0   0   4   8   7
    1   1   0   4   1   5   8
    2   1   0   7   8   6   3
    3   1   0   8   5   2   6

- ``mater.txt``

::

    1.0  0.3

- ``loads.txt``

::

    3  0.0  1.0
    6  0.0  2.0
    2  0.0  1.0

Run it in Python as follows:

.. code:: python

    import matplotlib.pyplot as plt  # load matplotlib
    from solidspy import solids_GUI  # import our package
    disp = solids_GUI()  # run the Finite Element Analysis
    plt.show()    # plot contours