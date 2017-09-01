from IR import IR


supercell=[[3, 0, 0],[0, 3, 0], [0, 0, 4]]


myIR=IR(PoscarName='POSCAR',BornFileName='BORN',ForceConstants=False,ForceFileName='FORCE_SETS',supercell=supercell,nac=False,masses=[])







myIR.write_spectrum('OscillatorStrengths.txt')
myIR.write_gaussiansmearedspectrum('SmearedOscillatorStrengths.txt',2)

