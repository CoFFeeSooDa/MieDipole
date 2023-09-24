# Mie Scattering with an Electric Dipole Source

[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

This is a MATLAB code to calculate the electromagnetic properties for a spherical scatterer illuminated by an electric dipole source based on the [Mie theory](https://en.wikipedia.org/wiki/Mie_scattering) in 1908.

This repository contains:

1. [MyTest_Tetra3_xyls.py](https://github.com/CoFFeeSooDa/StarTrackerTest/blob/main/MyTest_Tetra3_xyls.py): main program to test the robustness of tetra3.
2. [Centroid_Gen.py](https://github.com/CoFFeeSooDa/StarTrackerTest/blob/main/Centroid_Gen.py): generation of (star) centroids in a specified boresight.
3. [Database_Gen.py](https://github.com/CoFFeeSooDa/StarTrackerTest/blob/main/Database_Gen.py): generation of tetra3 database.
4. [Star_Catalog_to_mat.py](https://github.com/CoFFeeSooDa/StarTrackerTest/blob/main/Star_Catalog_to_mat.py): generation of .mat files from star catalog (Tycho-2) for [Centroid_Gen.py](https://github.com/CoFFeeSooDa/StarTrackerTest/blob/main/Centroid_Gen.py)
5. [My_Tetra3.py](https://github.com/CoFFeeSooDa/StarTrackerTest/blob/main/My_Tetra3.py): modified tetra3 solving engine of the compatibility of Tycho-2 catalog 

## Table of Contents

- [Installation and Usage](#Installation-and-Usage)
	- [About MyTest_Tetra3_xyls](#About-MyTest_Tetra3_xyls)
- [Customize Your Own Database](#Customize-Your-Own-Database)
	- [Generate Trimmed Star Table](#Generate-Trimmed-Star-Table)
	- [Generate Database](#Generate-Database)
- [Others](#Others)


## Installation and Usage

1. To use this repository, you can either download the zipped file or clone the repository (you need to install git).
```sh
git clone https://github.com/CoFFeeSooDa/StarTrackerTest.git
```
2. After the download, the script in ./catalogs/tycho2/combine.py should be run frist if you choose Tycho-2 as your star catalog.
Also, the Hipparcos catalog is also available in the catalog directory. (pending)

3. In the current version, the required files have been saved in the corresponding directory. In other words, you can run [MyTest_Tetra3_xyls.py](https://github.com/CoFFeeSooDa/StarTrackerTest/blob/main/MyTest_Tetra3_xyls.py) directly to test tetra3.

### About MyTest_Tetra3_xyls

The script provides four different tests:
1. Continuous test for ideal star centroids in a boresight with a specific FOV. The boresight is determined randomly.
2. Continuous test for ideal star centroids plus random pixel deviations, where we randomly assigned the pixel deviation and the target stars. The boresight is determined randomly.
3. Continuous test for ideal star centroids with additional distractors in a boresight with a specific FOV. The boresight is determined randomly.
4. Continuous test for ideal star centroids and delete stars randomly in a boresight with a specific FOV. The boresight is determined randomly.

```python
mode = 'continue_table' # 'continue_table', 'continue_table_devi', 'continue_table_addstar', and 'continue_table_delstar'
```

Moreover, the ideal star centroids are calculated by the trimmed star table which can be found in ./trimmed_table/*.mat.
To specify a trimmed star table, please change the following variable.
```python
read_mat_file = './trimmed_table/tycho2-Vmag06.mat'
```


##  Customize Your Own Database

To customize your own database for a specific use, you need to load a star catalog first. However, a star catalog is usually huge. 
Here I provide a script to trim the star catalog (Tycho-2) and extract the required information for the built database later.

### Generate Trimmed Star Table

You can create your customized trimmed star table (.mat files) using [Star_Catalog_to_mat.py](https://github.com/CoFFeeSooDa/StarTrackerTest/blob/main/Star_Catalog_to_mat.py). This script reads a star catalog (in ./catalogs) and trim the stars by setting the maximum Vmag. In the current version, the Tycho-2 catalog is available. 
In the script, please feel free to change the maximal magnitude (Johnson V filter) and the name of saved .mat file.
```python
star_max_magnitude = 6
s_name = f'tycho2-Vmag{star_max_magnitude:02d}.mat'
```
Further information about the [Tycho-2 catalog](https://cdsarc.u-strasbg.fr/ftp/cats/I/259/) can be found in ./catalog/tycho2/README.md.

### Generate Database

The Database for tetra3 can be generated by [Database_Gen.py](https://github.com/CoFFeeSooDa/StarTrackerTest/blob/main/Database_Gen.py).
Here are the best parameters obtained from the optimizations during the several months of my internship,
```python
t3.generate_database(max_fov=16.7, 
                     min_fov=None, # For single FOV, please set None
                     save_as='./database/my_tycho2_pspf17_pme005_vspf40_mag6',
                     path_catalog='catalogs/tycho2',
                     star_catalog='tyc2', # Use Tycho-2 catalog
                     pattern_stars_per_fov=17, # Crucial !!
                     verification_stars_per_fov=40, # Crucial!!
                     star_max_magnitude=6,
                     pattern_max_error=.005, # Crucial!! Associated with the bin number round(1/4/pattern_max_error)
                     simplify_pattern=False, # Do not change to True
                     range_ra=None,
                     range_dec=None,
                     presort_patterns=True, # True for faster matching speed
                     save_largest_edge=False, # True for faster matching speed
                     multiscale_step=1.5) # default 1.5
```
Further explanation for each parameter can be found in [tetra3](https://github.com/esa/tetra3).

## Others

In this section, the instructions for other files are put in here.

### Generate Centroids

Centroid_Gen.py is to compute the centroids of stars from a specific catalog (default is Tycho-2). You can import this function by using 
```python
import Centroid_Gen as gen
xyls = gen.generate_xyls(star_table_path, boresight, debug=False, fov=16.7, width=1124, height=1124)
```
Centroid_Gen.py can also run independently (usually served as debugging).
Cetnroid_Gen.py contains two functions, generate_xyls and calculate_vectors.
The latter function is called by generate_xyls.
Note that xyls is a numpy array of star-centroid pixels, star_table_path is the path of trimmed star table (.mat files, placed in ./trimmed table), and boresight is a list [RA, Dec] in degrees. 

### Tetra3 Solving Engine

My_Tetra3.py is the solving engine developed by ESA and adapted by me for my personal use. To correct run my script, please to do not use the python file in [https://github.com/esa/tetra3](https://github.com/esa/tetra3).

