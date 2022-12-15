# pyMOE
Python library for designing and producing masks for Micro Optical Elements. 

# Features
Mask design features: 
* Designing multi-layer (grayscale) masks from analytical and/or numerical (e.g. Gerchberg–Saxton) methods  
* Designing single layer (binary) dithered masks from grayscale masks 
* Designing rudimentary metasurfaces masks  

Mask production features: 
* Image file (binary/grayscale) <-> CAD files input-output conversion  
* Operations with/within mask files (e.g. change layer numbers, make intances, CAD file import, merge shapes, ...) 
* STL-to-GDS, convert 3D file into multi-layer GDS 


# Documentation

Follow the notebooks inside `notebooks\` directory


# Test

From root folder run
```
pytest
```    

