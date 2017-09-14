from IR import IR


supercell=[[1, 0, 0],[0, 1, 0], [0, 0, 1]]

#IF you would like to reproduce the frequencies calculated by DFPT in VASP and relied on the masses given in the POTCAR files:
masses=[12.011, 1.000, 1.000, 1.000, 1.000,  16.000]


myIR=IR(PoscarName='POSCAR',BornFileName='BORN',ForceConstants=True,ForceFileName='FORCE_CONSTANTS',supercell=supercell,nac=False,masses=masses)







myIR.write_spectrum('OscillatorStrengths.txt')
myIR.write_gaussiansmearedspectrum('SmearedOscillatorStrengths.txt',4)
myIR.plot_spectrum('OscillatorStrengths.eps')
myIR.plot_gaussiansmearedspectrum('SmearedOscillatorStrengths.eps',4)
