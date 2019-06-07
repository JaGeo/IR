from IRDipoleApprox import IR


supercell=[[3, 0, 0],[0, 3, 0], [0, 0, 4]]


myIR=IR(PoscarName='POSCAR',BornFileName='BORN',ForceConstants=False,ForceFileName='FORCE_SETS',supercell=supercell,nac=False,masses=[])







myIR.write_spectrum('OscillatorStrengths.txt')
myIR.write_gaussiansmearedspectrum('SmearedOscillatorStrengths.txt',4)
myIR.plot_spectrum('OscillatorStrengths.eps')
myIR.plot_gaussiansmearedspectrum('SmearedOscillatorStrengths.eps',4)
