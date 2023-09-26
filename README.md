# Generalized Mie Theory - Electric Dipole Source

[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

This MATLAB code calculates electromagnetic properties for a spherical scatterer illuminated by an electric dipole source based on the [Mie theory](https://en.wikipedia.org/wiki/Mie_scattering).

This repository contains:

1. [Main.m](https://github.com/CoFFeeSooDa/MieDipole/blob/main/Main.m): main program to initiate a calculation.
2. [InputFiles](https://github.com/CoFFeeSooDa/MieDipole/tree/main/InputFiles): input files, which are read by main.m (demonstration).
3. [Functions](https://github.com/CoFFeeSooDa/MieDipole/tree/main/Functions): function called by Main.m.

## Table of Contents

- [Installation and Usage](#Installation-and-Usage)
	- [About Main.m](#About-Main.m)
- [How to Set a Input File](#How-to-Set-a-Input-File)
	- [Generate Trimmed Star Table](#Generate-Trimmed-Star-Table)
	- [Generate Database](#Generate-Database)
- [Others](#Others)


## Installation and Usage

1. To use this code, you can either download the zipped file or clone the repository (you need to install Git in your PC first).
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
![](https://hackmd.io/_uploads/r1IwsGgga.png)



##  How to Set a Input File

To customize your own calculation, you can establish your settings in a JSON file. Here are some examples to demostrate what you can do in this program.

A standard JSON input file comprises three main objects, "Settings", "tmp_set" and "fplot". The first object defines the geometry of a system.
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

The second object defines dielectric properties of the space,
```json
"tmp_set"	 : {
	  "mode"	   : "Auto",
	  "lambda_s"   : 300e-9,
	  "lambda_e"   : 700e-9,
	  "epsi0"	   : 1,
	  "epsi1"	   : 4,
	  "epsi2"	   : ".\\InputFiles\\DielectricFunctions\\Ag_JPCL.csv"
    }
```
The third object defines the format of the output figures,
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

### Demo_WavelengthMode_CF_PCRET

In this JSON file, we calculate the coupling factor of a Ag core-shell sphere sandwiched by a donor and acceptor dipole.

The JSON file reads
```json
{
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
    },
  "tmp_set"	 : {
	  "mode"	   : "Auto",
	  "lambda_s"   : 300e-9,
	  "lambda_e"   : 700e-9,
	  "epsi0"	   : 1,
	  "epsi1"	   : 4,
	  "epsi2"	   : ".\\InputFiles\\DielectricFunctions\\Ag_JPCL.csv"
    },
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
}
```



## Others

