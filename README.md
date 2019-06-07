[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/JaGeo/IR/blob/master/LICENSE) [![DOI](https://zenodo.org/badge/101991065.svg)](https://zenodo.org/badge/latestdoi/101991065) [![Build Status](https://travis-ci.org/JaGeo/IR.svg?branch=master)](https://travis-ci.org/JaGeo/IR)

# IR
This python package can calculate infrared intensities based on the dipole approximation. To do so, you need [```VASP```](https://www.vasp.at/) and [```Phonopy```](https://github.com/atztogo/phonopy). 
<hr></hr>

What to cite
------------
It is based on the following two publications: 

1. P. Giannozzi, S. Baroni, *J. Chem. Phys.*, **1994**, *100*, 8537. 

2. D. Karhánek, T. Bučko, J. Hafner, *J. Phys.: Condens. Matter.*, **2010**, *22*, 265006.


 
They should be cited if you use the program. 

Moreover, the following should be cited:
1.  [A. Görne, J. George, J. van Leusen, R. Dronskowski, *Inorganics*, **2017**, *5*, 10.](https://doi.org/10.3390/inorganics5010010) 

2. J. George, & R. Dronskowski. (2018, February 7). IR Version 1.0.2 (Version v1.0.2). Zenodo. [http://doi.org/10.5281/zenodo.1168027](http://doi.org/10.5281/zenodo.1168027)  ([Bibtex](https://zenodo.org/record/1168027/export/hx)). 

Of course, also [```VASP```](https://www.vasp.at/) and [```Phonopy```](https://github.com/atztogo/phonopy).

Intallation
-----------
This package can be installed with ```pip install IR-JaGeo```.
Alternative Installation:
To use this package you need to install [```Phonopy```](https://github.com/atztogo/phonopy) correctly. Furthermore, ```numpy``` and ```matplotlib``` are required. Also, the python path should be exported correctly.

How to
--------
1. Perform a phonon calculation with Phonopy and VASP (finite displacements or DFPT) ([More information on this procedure](https://atztogo.github.io/phonopy/procedure.html))
2. Generate the ```FORCE_SETS``` or ```FORCE_CONSTANTS``` file
3. Calculate BORN charges ([More information on this procedure and the ```BORN``` file](https://atztogo.github.io/phonopy/vasp.html))
4. Download this repository, export the Python path correctly
5. Copy an example script, adapt the names of the files and the supercell size (the one you used for the phonon calculation!)
6. Run the script


Todo
--------
1. Other functionalities
2. Include tests

Information about the Author
--------

- J. George (Université catholique de Louvain, before: RWTH Aachen University)
- PI during the development of the code: R. Dronskowski, RWTH Aachen University

