.. pyMOE documentation master file, created by
   sphinx-quickstart on Wed Aug 23 11:35:11 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

	
.. toctree::
   :maxdepth: 4


Welcome to pyMOE's documentation!
=================================

pyMOE is a Python library for designing and producing masks for Micro Optical Elements. 

Features
********

Mask design features: 

* Designing multi-layer (grayscale) masks from analytical and/or numerical (e.g. Gerchberg–Saxton) methods  

* Designing single layer (binary) dithered masks from grayscale masks 

* Designing rudimentary metasurfaces masks 

* Numerical propagation methods (Rayleigh-Sommerfeld, Fresnel, Fraunhofer) of electric fields from (grayscale) masks 

Mask production features: 

* Image file (binary/grayscale) <-> CAD files input-output conversion  

* Operations with/within mask files (e.g. change layer numbers, make intances, CAD file import, merge shapes, ...)  

Dependencies
************

A list of dependencies (with versions) is at https://github.com/INLnano/pyMOE/blob/main/requirements.txt . To use the package straighforwardly please make sure you have those dependencies (and versions) installed. 

Modular architechture 
*********************

pyMOE is based on a modular architechture, where each module provides specific functionalities. The use of each module and its inter-modular connections is exemplified in notebooks (see Notebooks, there is one notebook for each module). 

Example Notebooks 
*****************

Check the notebooks in the sidebar menu.

Full function/class listings 
****************************

Please check the function/class Indices and Tables in the sidebar menu or use Search functionality.

.. toctree::
   notebooks
   listings