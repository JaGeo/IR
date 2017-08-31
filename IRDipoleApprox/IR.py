import numpy as np
from phonopy import Phonopy
from phonopy.interface.vasp import read_vasp
from phonopy.file_IO import parse_FORCE_SETS, parse_BORN, write_FORCE_CONSTANTS, parse_FORCE_CONSTANTS
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.units import VaspToCm, VaspToTHz
import os



class IR:
	
	def __init__(self,PoscarName='POSCAR',BornFileName='BORN',ForceConstants=False,ForceFileName='FORCE_SETS',supercell=[[1, 0, 0],[0, 1, 0], [0, 0, 1]]
,nac=False,masses=[]):
		"""
		Class for calculating the IR spectra in the dipole approximation according to:
		P. Giannozzi and S. Baroni, J. Chem. Phys. 100, 8537 (1994). 
		
		and
	
		D. Karhánek, T. Bučko, J. Hafner, J. Phys.: Condens. Matter 22 265006 (2010).  
		(http://homepage.univie.ac.at/david.karhanek/downloads.html )

		This class was also carefully tested against the script by D. Karhánek. 
	
		Args:
			PoscarNamse (str): name of the POSCAR that was used for the phonon calculation
			BornFileName (str): name of the file with BORN charges (formatted with outcar-born)
			ForceConstants (boolean): If True, ForceConstants are read in. If False, forces are read in.
			ForceFileName (str): name of the file including force constants or forces
			supercell (list of lists): reads in supercell
			nac (boolean): If true, NAC is applied.
			masses (list): Masses in this list are used instead of the ones prepared in Phonopy. Useful for isotopes.
		"""

		
		self.__unitcell =read_vasp(PoscarName)
		self.__supercell=supercell
		self.__phonon= Phonopy(self.__unitcell,supercell_matrix=self.__supercell,factor=VaspToCm,symprec=1e-4)
	        self.__natoms=self.__phonon.get_primitive().get_number_of_atoms()	

		#If different masses are supplied
		if masses: 
			self.__phonon.set_masses(masses)
                self.__masses=self.__phonon.get_primitive().get_masses() 	
		#Forces or Force Constants
		if not ForceConstants:
			self.__force_sets = parse_FORCE_SETS(filename=ForceFileName)
			self.__phonon.set_displacement_dataset(self.__force_sets)
			self.__phonon.produce_force_constants()
		
		if ForceConstants:
			force_constants = parse_FORCE_CONSTANTS(filename=ForceFileName)
			self.__phonon.set_force_constants(force_constants)
	
		#Read in BORN file
		BORN_file = parse_BORN(self.__phonon.get_primitive(),filename=BornFileName)
                self.__BORN_CHARGES=BORN_file['born']

		#Apply NAC Correction
                if nac:
                        self.__phonon.set_nac_params(BORN_file)
		self.__frequencies,self.__eigvecs=self.__phonon.get_frequencies_with_eigenvectors([0, 0, 0])


		self.__NumberOfBands=len(self.__frequencies)
		#Nicer format of the eigenvector file
		self.__FormatEigenvectors()
		#Get dipole approximation of the intensitiess
		self.__set_intensities()
		
	


	def __FormatEigenvectors(self):
		"""
		Formats eigenvectors to a dictionary: the first argument is the number of bands, the second the number of atoms, the third the Cartesian coordinate
		"""

		self.__EigFormat = {}
		for alpha in  range(self.__NumberOfBands):
			laufer=0
    	 	   	for beta in range(self.__natoms):
        			for xyz in range(0,3):
      	    	     	 		self.__EigFormat[beta,alpha,xyz]=self.__eigvecs[laufer][alpha]
                       			laufer=laufer+1

	def __Eigenvector(self, atom, band, xoryorz ):
		"""
		Gives a certain eigenvector corresponding to one specific atom, band and Cartesian coordinate

		args:
			
			
		"""

		return np.real(self.__EigFormat[atom,band,xoryorz])

	def __massEig(self,atom,band,xoryorz):
		"""
		Gives a certain eigenvector corresponding to one specific atom, band and Cartesian coordinate

			

		"""

		return self.__Eigenvector(atom,band,xoryorz)/np.sqrt(self.__masses[atom])

	def __set_intensities(self):
		Intensity={}
		for freq in range(len(self.__frequencies)):
			Intensity[freq]=0
			for alpha in range(3):
				sum=0
				for l in range(self.__natoms):
                			for beta in range(3):
						sum=sum+self.__BORN_CHARGES[l,alpha,beta]*self.__massEig(l,freq,beta)
				Intensity[freq]=Intensity[freq]+np.power(np.absolute(sum),2)

		ReformatIntensity=[]		
		for i in Intensity:
			ReformatIntensity.append(Intensity[i])	
			
		self.__Intensity=np.array(ReformatIntensity)
	
	#moeglicherweise 3 Translationen und Imaginaermoden auslassen!	
	def get_intensities(self):
		return self.__Intensity
	#moeglicherweise 3 Translationen und Imaginaermoden auslassen!	
	def get_frequencies(self):
		return self.__frequencies

	#hier wird bisher keine Entartung beachtet
	def get_spectrum(self):
		spectrum = {'Frequencies': self.__frequencies, 'Intensities': self.__Intensity }
		return spectrum

	def get_gaussiansmearedspectrum(self,sigma):
		unsmearedspectrum=self.get_spectrum()
		frequencies=self.get_frequencies()
		Intensity=self.get_intensities()
		rangex=np.linspace(0,np.nanmax(frequencies),num=int(np.nanmax(frequencies))*100)
		y=np.zeros(int(np.nanmax(frequencies))*100)
		for i in range(len(frequencies)):
			y=y+self.__gaussiansmearing(rangex,frequencies[i],Intensity[i],sigma,np.nanmax(frequencies))
		smearedspectrum={'Frequencies': rangex, 'Intensities': y}
		return smearedspectrum



	def __gaussiansmearing(self,rangex,frequency,Intensity,sigma,numberofFrequencies):
		y=np.zeros(int(numberofFrequencies)*100)		
		y=Intensity*np.exp(-np.power((rangex-frequency),2)/(2*np.power(sigma,2)))*np.power(np.sqrt(2*np.pi)*sigma,-1)
		return y
	

	def write_spectrum(self,filename):
		spectrum=self.get_spectrum()
		self.__write_file(filename,spectrum)	

	def write_gaussiansmearedspectrum(self,filename,sigma):
		spectrum=self.get_gaussiansmearedspectrum(sigma)
		self.__write_file(filename,spectrum)

	def __write_file(self,filename,spectrum):
		Freq=np.array(spectrum['Frequencies'].tolist())  
		Intens=np.array(spectrum['Intensities'].tolist())
		file  = open(filename, 'w')
		for i in range(len(Freq)):
			file.write('%s %s \n' % (Freq[i], Intens[i]))
		file.close()				

