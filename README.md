# pyMOE
Python library for designing and producing masks for Micro Optical Elements. 


# Features
Mask design features: 
* Designing multi-layer (grayscale) masks from analytical and/or numerical (e.g. Gerchberg–Saxton) methods  
* Designing single layer (binary) dithered masks from grayscale masks 
* Designing rudimentary metasurfaces masks 
* Numerical propagation methods (Rayleigh-Sommerfeld, Fresnel, Fraunhofer) of electric fields from (grayscale) masks 

Mask production features: 
* Image file (binary/grayscale) <-> CAD files input-output conversion  
* Operations with/within mask files (e.g. change layer numbers, make intances, CAD file import, merge shapes, ...)  


# Getting started

Follow the notebooks inside `notebooks\` directory


# Documentation

For documentation on classes and functions within package's modules, please visit https://pymoe-doc.readthedocs.io 

# Dependencies

A list of dependencies (with versions) is at https://github.com/INLnano/pyMOE/blob/main/requirements.txt . To use the package straighforwardly please make sure you have those dependencies (and versions) installed. 

# Test

From root folder run
```
pytest
```    

