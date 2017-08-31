from IR import IR


supercell=[[1, 0, 0],[0, 1, 0], [0, 0, 1]]
masses=[12.011, 1.000, 1.000, 1.000, 1.000,  16.000]
new=IR(PoscarName='Testfiles/POSCAR',BornFileName='Testfiles/BORN',ForceConstants=True,ForceFileName='Testfiles/FORCE_CONSTANTS',supercell=supercell,nac=False,masses=masses)


#print(new.get_intensities())
#print(new.get_frequencies())
#print(new.get_spectrum())

#print(new.get_gaussiansmearedspectrum(10))
new.write_spectrum('Outputs/NormalSpectrum2.txt')
new.write_gaussiansmearedspectrum('Outputs/Gaussian2.txt',2)

