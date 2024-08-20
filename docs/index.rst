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

* Numerical scalar propagation methods (Rayleigh-Sommerfeld, Fresnel, Fraunhofer) of electric fields from (grayscale) masks 

Mask production features: 

* Image file (binary/grayscale) <-> CAD files input-output conversion  

* Operations with/within mask files (e.g. change layer numbers, make intances, CAD file import, merge shapes, ...)  


Dependencies
************

A list of dependencies (with versions) is at https://github.com/INLnano/pyMOE/blob/main/requirements.txt . To use the package straighforwardly please make sure you have these dependencies (and versions) installed. 


Modular architecture 
*********************

pyMOE is based on a modular architecture, where each module provides specific functionalities. The use of each module and its inter-modular connections is exemplified at  https://www.sciencedirect.com/science/article/pii/S0010465524002546#fg0060 and practically in the example notebooks (there is one notebook for each module). 


Literature  
**********

For background, theory, package structure, practical results and discussion of package's performance, please refer to the publication: 

J. Cunha, J. Queiroz, C. Silva, F. Gentile and D. E. Aguiam, pyMOE: Mask Design and Modelling for Micro Optical Elements and Flat Optics, Computer Physics Communications 109331 (2024). Link: https://doi.org/10.1016/j.cpc.2024.109331


Example Notebooks 
*****************

Check the notebooks in the sidebar menu.

.. toctree::
   notebooks

Full function/class listings 
****************************

Please check the function/class Indices and Tables in the sidebar menu or use Search functionality.

.. toctree::
   listings