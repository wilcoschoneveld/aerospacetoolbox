Aerospace Toolbox
================

The Aerospace Toolbox for Python contains a collection of functions intended for aerospace engineering applications. The toolbox is developed (and currently tested only) in a Python 2.7 environment and requires a recent version of NumPy and SciPy. The toolbox is designed to allow quick and robust evaluation of different functions and to properly accept and handle array-like data structures.

### Features

* Functions for evaluation of the standard atmosphere and custom defined atmospheres. Isentropic relations, normal shock wave relations and Prandtl-Meyer expansion waves. And more!
* The theoretical formulas are evaluated by either direct substition or iteration with the Newton-Raphson method if the former is not possible.
* Each function has a full docstring explaining it's syntax and usage. Docstrings follow the NumPy docstring standard.
* All functions accept input data in any form that can be converted to an array. This includes scalars, lists, tuples, matrices and ndarrays. Each functions returns it's results in original form.
* Includes a module which wraps most of the functions in [MATLAB Aerospace Toolbox](http://www.mathworks.nl/help/aerotbx/index.html) syntax.

### Installation

The easiest way to obtain the most recent version of the Aerospace Toolbox is by downloading the source as a [ZIP-file](https://github.com/wilcoschoneveld/aerospacetoolbox/zipball/master). Extract the contents to a folder and install the package by running the following command:

````
python setup.py install
````

A Python 2.7 distribution with [NumPy](http://www.numpy.org/) and [SciPy](http://www.scipy.org/) is needed in order to make proper use of the Aerospace Toolbox. A recent version of [matplotlib](http://matplotlib.org/) is required to run the example scripts.

### Usage

If you want to use the Aerospace Toolbox in your scripts, you can write the following command to import all the functions:

````python
from aerotbx import *
````

This will import the main functions and modules into your namespace. From here you can start using the functions provided by the toolbox. In the case you come from MATLAB and are used to the syntax from it's aerospace toolbox, you can use the matlab module that comes with this package. You can import all the MATLAB-style functions with the following:

````python
from aerotbx.matlab import *
````

### Examples

*to be added*
