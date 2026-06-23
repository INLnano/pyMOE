# pyMOE
Python library for designing and producing masks for Micro Optical Elements. 


# Features
Mask design features: 
* Designing multi-layer (grayscale) masks from analytical and/or numerical (e.g. Gerchberg–Saxton) methods  
* Designing single layer (binary) dithered masks from grayscale masks 
* Designing rudimentary metasurfaces masks 
* Numerical scalar propagation methods (Rayleigh-Sommerfeld, Fresnel, Fraunhofer, Angular Spectrum Method, Bluestein) of electric fields from (grayscale) masks
* Optimization and inverse design based on gradient-descent methods

Mask production features: 
* Image file (binary/grayscale) <-> CAD files input-output conversion  
* Operations with/within mask files (e.g. change layer numbers, make intances, CAD file import, merge shapes, ...)  


# Getting started

For background, theory, package structure, practical results and discussion of package's performance, please refer to the publication: https://doi.org/10.1016/j.cpc.2024.109331.

For optimization and inverse design using pyMOE, please refer to the proceedings paper https://doi.org/10.1117/12.3099849.

For a practical walk through on how to use the package, besides the aforementioned publications, please follow the example notebooks inside `notebooks\` directory. For running we recommend having Jupyter Notebook or Lab with a Python distribution >3.7, as well as required dependencies. 


# Documentation

For documentation please visit https://pymoe-doc.readthedocs.io 


# Dependencies

A list of dependencies (with versions) is at https://github.com/INLnano/pyMOE/blob/main/requirements.txt . To use the package straighforwardly please make sure you have these dependencies (and versions) installed. 


# Citing pyMOE  

If pyMOE was useful for you, please cite the companion paper published in Computer Physics Communications (https://doi.org/10.1016/j.cpc.2024.109331), bibtex entry:

@article{CUNHA2024109331,\
author = {Joao Cunha and José Queiroz and Carlos Silva and Fabio Gentile and Diogo E. Aguiam},\
title = {pyMOE: Mask Design and Modelling for Micro Optical Elements and Flat Optics},\
journal = {Computer Physics Communications},\
pages = {109331},\
year = {2024},\
issn = {0010-4655},\
doi = {https://doi.org/10.1016/j.cpc.2024.109331 }}

And/or the proceedings paper dedicated to inverse design and optimization using pyMOE (https://doi.org/10.1117/12.3099849), bibtex entry:

@inproceedings{10.1117/12.3099849,\
author = {Joao Cunha and Filipe Camarneiro and Diogo E. Aguiam},\
title = {{Inverse design framework for large-area micro-optical elements enabling beam shaping and beam steering applications}},\
volume = {14088},\
booktitle = {Optical and Photonic Microsystems: Technology and Applications},\
editor = {Hans Zappe and Niels Quack},\
organization = {International Society for Optics and Photonics},\
publisher = {SPIE},\
pages = {140880B},\
year = {2026},\
doi = {https://doi.org/10.1117/12.3099849 }}


# Zenodo repository 

pyMOE releases are also mirrored in Zenodo at https://doi.org/10.5281/zenodo.7326572. 


# Interactions and Contributions 

To interact/contribute please use the following: 
* **Email contact**: Please feel free to contact via email - for email addresses, please use the email(s) provided at https://doi.org/10.1016/j.cpc.2024.109331.
* **Issues page**: Please feel free to post suggestions, bug reports, code malfunctions, or requests in the Issues page (https://github.com/INLnano/pyMOE/issues). Those will be cathegorized and addressed in due time. 
* **Code contributions**: To contribute, please (1) fork the repository, (2) commit your code changes in your forked repository, and (3) make a pull request *to a dev-(...) branch* of this pyMOE repository (*pull requests directly to pyMOE's main branch will not be accepted*).  


# Test

From root folder run
```
pytest
```    
