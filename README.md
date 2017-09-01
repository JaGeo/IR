# IR

This python package can calculate infrared intensities based on the dipole approximation. To do so, you need [```VASP```](https://www.vasp.at/) and [```Phonopy```](https://github.com/atztogo/phonopy). 
<hr></hr>
It is based on the following two publications: 

1. P. Giannozzi, S. Baroni, *J. Chem. Phys.*, **1994**, *100*, 8537. 

2. D. Karhánek, T. Bučko, J. Hafner, *J. Phys.: Condens. Matter.*, **2010**, *22*, 265006.

They should be cited if you use the program. Of course, also [```VASP```](https://www.vasp.at/) and [```Phonopy```](https://github.com/atztogo/phonopy).

Intallation
-----------
To use this package you need to install [```Phonopy```](https://github.com/atztogo/phonopy) correctly. Furthermore, ```numpy``` and ```matplotlib``` are required. Also, the python path should be exported correctly.

How to
--------
1. Perform a phonon calculation with Phonopy and VASP (finite displacements or DFTP) ([More information on this procedure](https://atztogo.github.io/phonopy/procedure.html))
2. Generate the ```FORCE_SETS``` or ```FORCE_CONSTANTS``` file
3. Calculate BORN charges ([More information on this procedure](https://atztogo.github.io/phonopy/procedure.html)) and the ```BORN``` file
4. Download this repository, export the Python path correctly
5. Copy an example script, adapt the names of the files and the supercell size (the one you used for the phonon calculation!)
6. Run the script



