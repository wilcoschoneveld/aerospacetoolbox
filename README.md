aerospacetoolbox
================

Functions for aerospace analysis to develop and evaluate your designs. Currently contains the following:

Environment (accepts multidimensional arrays as input)
- atmosisa: Evaluate the international standard atmosphere (ISA) at a given altitude. The function assumes a continued troposphere below 0 meters and an infinite mesosphere above 71 kilometers geopotential height.
- geoidheight: Calculates the geoid height using the EGM96 Geopotential Model.

Gas Dynamics (can convert multidimensional arrays from any property to another):
- flowisentropic: Isentropic relations with a given set of specific heat ratios and any one of the isentropic flow variables.
- flownormalshock: Normal shock relations with a given set of specific heat ratios and any one of the normal shock variables.
- flowprandtlmeyer: Prandtl-Meyer function for expansion waves.

Unit Conversion (can convert multidimensional arrays from any unit to another):
- convert: conversion of units, such as mass, pressure and density.

Created to (partially) look like the equivalent set of functions in the matlab
aerospace toolbox: http://www.mathworks.nl/help/aerotbx/index.html