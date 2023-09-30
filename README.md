# Generalized Mie Theory - Electric Dipole Source

[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

This MATLAB code calculates electromagnetic properties for a spherical scatterer illuminated by an electric dipole source based on the [Mie theory](https://en.wikipedia.org/wiki/Mie_scattering).

This repository contains:

1. [Main.m](https://github.com/CoFFeeSooDa/MieDipole/blob/main/Main.m): main program to initiate a calculation.
2. [InputFiles](https://github.com/CoFFeeSooDa/MieDipole/tree/main/InputFiles): input files, which are read by main.m (demonstration).
3. [Functions](https://github.com/CoFFeeSooDa/MieDipole/tree/main/Functions): function called by Main.m.


![](https://i.imgur.com/V6hKZPp.png)


## Table of Contents

- [Installation and Usage](#Installation-and-Usage)
	- [About Main.m](#About-Main.m)
- [How to Set an Input File](#How-to-Set-an-Input-File)
	- [First Part](#First-Part)
	- [Second Part](#Second-Part)
	- [Third Part](#Third-Part)
- [Citation](#Citation)


## Installation and Usage

1. To use this code, you can either download the zipped file or clone the repository (you need to install Git on your PC first).
```sh
git clone https://github.com/CoFFeeSooDa/MieDipole.git
```
2. After the download, you can start a calculation in MATLAB (2016a and after is recommended).


### About Main.m
Main.m starts a calculation by reading a JSON file in lines 22-23,
```
FilePath = './InputFiles/'; % Folder Path of Input Files
FileName = 'Demo_WavelengthMode_CF_PCRET'; % File Name
```
![](https://i.imgur.com/ogu2RwZ.png)



##  How to Set an Input File

To customize your own calculation, you can establish your settings in a JSON file. Here are some examples to demonstrate what you can do in this program.

A standard JSON input file comprises three main objects, "Settings", "tmp_set" and "fplot".

### **First Part**
**The first object defines the geometry of a system.**
```json
"Settings" : {
	  "ModeName"   : "wavelength",  
	  "Quantity"   : "CF",          
	  "DPos"	   : {
		  "Cart"   : [0, 0, 80e-9]
	  },
	  "APos"	   : {
		  "Cart"   : [0, 0, -80e-9]
	  },
	  "DOri"	   : {
		  "Cart"   : [0, 0, 1]
	  },
	  "AOri"	   : {
		  "Cart"   : [0, 0, 1]
	  },
	  "nmax"	   : 70,
	  "BC"		   : "coreshell",
	  "rbc"		   : [70e-9,60e-9],
	  "Dpstrength" : 1
    }
```
1. ```ModeName```: sets the calculation type.
Currently, the code supports the following types: ```wavelength``` , ```angle``` and ```mapping```. In the ```wavelength``` mode, the program sweeps the wavelength (wavenumber) of the donor dipole; in the ```angle``` mode, the program sweeps the position of the acceptor dipole, which is described by $\theta$; in the ```mapping``` mode, the program plot the mapping of a specific plane (see also ```tmp_set```).
2. ```Quantity```: specifies the quantity presented in the figures. 
Currently, ```CF``` (coupling factor, the definition can be found in [here](https://pubs.acs.org/doi/full/10.1021/acs.jpclett.0c01989)) and ```Purcell```(Purcell factor, the enhancement of spontaneous emission) are available.
3. ```DPos```(```APos```): sets the position of the donor (acceptor) dipole.(unit: meter)
4. ```DOri```(```AOri```): set the orientation of the donor (acceptor) dipole. (dimensionless, be sure that it is normalized)
5. ```nmax```: the highest expansion order in a calculation.
6. ```BC```: the boundary condition --- ```single``` or ```coreshell```.
7. ```rbc```: the radius of a single (core/shell) sphere (descending order)
8. ```Dpstrength```: the strength of the donor dipole (Unit: cgs)

### **Second Part**
**The second object defines dielectric properties of the space,**
```json
"tmp_set"	 : {
	  "mode"       : "Auto",
	  "lambda_s"   : 300e-9,
	  "lambda_e"   : 700e-9,
	  "epsi0"      : 1,
	  "epsi1"      : 4,
	  "epsi2"      : ".\\InputFiles\\DielectricFunctions\\Ag_JPCL.csv"
    //(for angle mode only)-----------------------------------------------
      "ThetaResol" : 1,
	  "Theta_i"    : 5,
	  "Theta_f"    : 180,
	  "Ar"         : 80e-9,
	  "Phi"        : 0
    //(for mapping mode only)---------------------------------------------
      "x_points"   : 101,
	  "x_start"    : -300e-9,
	  "x_end"      :  300e-9,
	  "y_points"   : 101,
	  "y_start"    : -300e-9,
	  "y_end"      :  300e-9,
	  "plane"      : "xz",
	  "third_coord": 0
    }
```
1. ```mode```: control of the fitting of dielectric functions.
2. ```lambda_s```: the starting point of wavelength (unit: m)
3. ```lambda_e```: the end point of wavelength (unit:m)
4. ```epsi0```: the dielectric function of the outermost region.
5. ```epsi1```: the dielectric function of the middle region.
6. ```epsi2```: the dielectric function of the innermost region. (for ```"BC"="coreshell"``` only)
7. ```ThetaResol```: the resolution of $\theta$ (the position of the acceptor dipole, unit: degree)
8. ```Theta_i```: the initial value of $\theta$. (unit: degree)
9. ```Theta_f```: the final value of $\theta$. (unit: degree)
10. ```Ar```: the $r$ of the acceptor dipole (unit: m)
11. ```Phi```: the $\phi$ of the acceptor dipole (unit: degree)
12. ```x_points```: the number of grid points along x axis.
13. ```x_start```: the starting point of the x-axis (unit: m)
14. ```x_end```: the endpoint of the x-axis (unit: m)
15. ```y_points```: the number of grid points along the y-axis.
16. ```y_start```: the starting point of the y-axis (unit: m)
17. ```y_end```: the endpoint of the y-axis (unit: m)


### **Third Part**
**The third object defines the format of the output figures,**
```json
"fplot"	 : {
	  "colorstyle" : "-k",
	  "range" 	   : [null,null,null,null],
	  "yscale"     : "log",
	  "xlabel"	   : "$\\mathrm{Wavenumber}~(\\mathrm{cm}^{-1})$",
	  "ylabel"	   : "$\\mathrm{Coupling~Factor}~(\\mathrm{cm}^{-6})$",
	  "subaxis"    : 1,
	  "subrange"   : [null,null,null,null],
	  "subxlabel"  : "$\\mathrm{Wavelength}~(\\mathrm{nm})$"
    }	  
```
It is recommended that the readers should be familiar with the plot setting in MATLAB before changing the settings. Otherwise, the setting should be kept in default.


## Citation
If you use MieDipole in your work, please cite:

> **Ming-Wei Lee** and Liang-Yan Hsu,
> "Controllable Frequency Dependence of Resonance Energy Transfer Coupled with Localized Surface Plasmon Polaritons",
> *J. Phys. Chem. Lett.* 2020, **11**, 16, 6796–6804


```bib
@article{Lee2020,
    author = {Lee, Ming-Wei and Hsu, Liang-Yan},
    doi = {10.1021/acs.jpclett.0c01989},
    issn = {1948-7185},
    journal = {J. Phys. Chem. Lett.},
    number = {16},
    pages = {6796--6804},
    pmid = {32787214},
    url = {https://pubs.acs.org/doi/10.1021/acs.jpclett.0c01989},
    volume = {11},
    year = {2020}
}
```

