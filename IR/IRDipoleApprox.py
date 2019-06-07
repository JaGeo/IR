import numpy as np
import matplotlib.pyplot as plt
from phonopy import Phonopy
from phonopy.interface.vasp import read_vasp
from phonopy.file_IO import parse_FORCE_SETS, parse_BORN, parse_FORCE_CONSTANTS
from phonopy.units import VaspToCm
from phonopy.phonon.degeneracy import degenerate_sets as get_degenerate_sets


class IR:
    def __init__(self, PoscarName='POSCAR', BornFileName='BORN', ForceConstants=False, ForceFileName='FORCE_SETS',
                 supercell=[[1, 0, 0], [0, 1, 0], [0, 0, 1]]
                 , nac=False, symprec=1e-5, masses=[], primitive=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                 degeneracy_tolerance=1e-5):
        """
        Class for calculating the IR spectra in the dipole approximation according to:
        P. Giannozzi and S. Baroni, J. Chem. Phys. 100, 8537 (1994).
        and
        D. Karhanek, T. Bucko, J. Hafner, J. Phys.: Condens. Matter 22 265006 (2010).


        This class was also carefully compared to the results of the script by D. Karhanek available at http://homepage.univie.ac.at/david.karhanek/downloads.html

        args:
            PoscarNamse (str): name of the POSCAR that was used for the phonon calculation
            BornFileName (str): name of the file with BORN charges (formatted with outcar-born)
            ForceConstants (boolean): If True, ForceConstants are read in. If False, forces are read in.
            ForceFileName (str): name of the file including force constants or forces
            supercell (list of lists): reads in supercell
            nac (boolean): If true, NAC is applied.
            symprec (float): contains symprec tag as used in Phonopy
            masses (list): Masses in this list are used instead of the ones prepared in Phonopy. Useful for isotopes.
            primitive (list of lists): contains rotational matrix to arrive at primitive cell
            degeneracy_tolerance (float): tolerance for degenerate modes
        """

        self.__unitcell = read_vasp(PoscarName)
        self.__supercell = supercell
        self.__phonon = Phonopy(self.__unitcell, supercell_matrix=self.__supercell, primitive_matrix=primitive,
                                factor=VaspToCm, symprec=symprec)
        self.__natoms = self.__phonon.get_primitive().get_number_of_atoms()
        self._degeneracy_tolerance = degeneracy_tolerance
        # If different masses are supplied
        if not ForceConstants:
            self.__force_sets = parse_FORCE_SETS(filename=ForceFileName)
            self.__phonon.set_displacement_dataset(self.__force_sets)
            self.__phonon.produce_force_constants()

        if ForceConstants:
            force_constants = parse_FORCE_CONSTANTS(filename=ForceFileName)
            self.__phonon.set_force_constants(force_constants)


        if masses:
            self.__phonon._build_supercell()
            self.__phonon._build_primitive_cell()
            # if ForceConstants:
            #     force_constants = parse_FORCE_CONSTANTS(filename=ForceFileName)
            #     self.__phonon.set_force_constants(force_constants)
            self.__phonon.set_masses(masses)

        self.__masses = self.__phonon.get_primitive().get_masses()
        # Forces or Force Constants

        # Read in BORN file
        BORN_file = parse_BORN(self.__phonon.get_primitive(), filename=BornFileName)
       
        self.__BORN_CHARGES = BORN_file['born']

        # Apply NAC Correction
        if nac:
            self.__phonon.set_nac_params(BORN_file)
        self._frequencies, self._eigvecs = self.__phonon.get_frequencies_with_eigenvectors([0, 0, 0])

        self.__NumberOfBands = len(self._frequencies)

        # Nicer format of the eigenvector file
        self.__FormatEigenvectors()

        # Get dipole approximation of the intensitiess
        self.__set_intensities()

    def __FormatEigenvectors(self):
        """
        Formats eigenvectors to a dictionary: the first argument is the number of bands, the second the number of atoms, the third the Cartesian coordinate
        """

        self._EigFormat = {}
        for alpha in range(self.__NumberOfBands):
            laufer = 0
            for beta in range(self.__natoms):
                for xyz in range(0, 3):
                    self._EigFormat[beta, alpha, xyz] = self._eigvecs[laufer][alpha]
                    laufer = laufer + 1

    def _Eigenvector(self, atom, band, xoryorz):
        """
        Gives a certain eigenvector corresponding to one specific atom, band and Cartesian coordinate

        args:
            atom (int) : number of the atoms (same order as in POSCAR)
            band (int) : number of the frequency (ordered by energy)
            xoryorz (int): Cartesian coordinate of the eigenvector


        """

        return np.real(self._EigFormat[atom, band, xoryorz])

    def __massEig(self, atom, band, xoryorz):
        """
        Gives a certain eigenvector divided by sqrt(mass of the atom) corresponding to one specific atom, band and Cartesian coordinate

        args:
            atom (int) : number of the atoms (same order as in POSCAR)
            band (int) : number of the frequency (ordered by energy)
            xoryorz (int): Cartesian coordinate of the eigenvector


        """

        return self._Eigenvector(atom, band, xoryorz) / np.sqrt(self.__masses[atom])

    def __set_intensities(self):
        """
        Calculates the oscillator strenghts according to "P. Giannozzi and S. Baroni, J. Chem. Phys. 100, 8537 (1994)."
        """

        Intensity = {}
        for freq in range(len(self._frequencies)):
            Intensity[freq] = 0
            for alpha in range(3):
                sum = 0
                for l in range(self.__natoms):
                    for beta in range(3):
                        sum = sum + self.__BORN_CHARGES[l, alpha, beta] * self.__massEig(l, freq, beta)
                Intensity[freq] = Intensity[freq] + np.power(np.absolute(sum), 2)

        # get degenerate modes
        freqlist_deg = get_degenerate_sets(self._frequencies, cutoff=self._degeneracy_tolerance)

        ReformatIntensity = []
        for i in Intensity:
            ReformatIntensity.append(Intensity[i])

        # if degenerate modes exist:
        if (len(freqlist_deg) < len(self._frequencies)):

            Intensity_deg = {}
            for sets in range(len(freqlist_deg)):
                Intensity_deg[sets] = 0
                for band in range(len(freqlist_deg[sets])):
                    Intensity_deg[sets] = Intensity_deg[sets] + ReformatIntensity[freqlist_deg[sets][band]]

            ReformatIntensity = []
            for i in range(len(Intensity_deg)):
                ReformatIntensity.append(Intensity_deg[i])

            Freq = []
            for band in range(len(freqlist_deg)):
                Freq.append(self._frequencies[freqlist_deg[band][0]])

            self.__frequencies_deg = np.array(Freq)

        else:
            self.__frequencies_deg = self._frequencies

        self.__Intensity = np.array(ReformatIntensity)

    def get_intensities(self):
        """
        returns calculated oscillator strengths as a numpy array
        """

        return self.__Intensity

    def get_frequencies(self):
        """
        returns frequencies as a numpy array
        """

        return self.__frequencies_deg

    def get_spectrum(self):
        """
        returns spectrum as a dict of numpy arrays
        """

        """
        Degeneracy should be treated for the Oscillator strengths

        """

        spectrum = {'Frequencies': self.get_frequencies(), 'Intensities': self.get_intensities()}
        return spectrum

    # only gaussian broadening so far
    def get_gaussiansmearedspectrum(self, sigma):
        """
        returns a spectrum with gaussian-smeared intensities

        args:
            sigma (float): smearing
        """

        #unsmearedspectrum = self.get_spectrum()
        frequencies = self.get_frequencies()
        Intensity = self.get_intensities()
        rangex = np.linspace(0, np.nanmax(frequencies) + 50, num=int(np.nanmax(frequencies) + 50) * 100)
        y = np.zeros(int(np.nanmax(frequencies) + 50) * 100)
        for i in range(len(frequencies)):
            y = y + self.__gaussiansmearing(rangex, frequencies[i], Intensity[i], sigma)
        smearedspectrum = {'Frequencies': rangex, 'Intensities': y}
        return smearedspectrum

    def __gaussiansmearing(self, rangex, frequency, Intensity, sigma):
        """
        applies gaussian smearing to a range of x values, a certain frequency, a given intensity

        args:
            rangex (ndarray): Which values are in your spectrum
            frequency (float): frequency corresponding to the intensity that will be smeared
            Intensity (float): Intensity that will be smeared
            sigma (float): value for the smearing

        """

        y = np.zeros(rangex.size)
        y = Intensity * np.exp(-np.power((rangex - frequency), 2) / (2 * np.power(sigma, 2))) * np.power(
            np.sqrt(2 * np.pi) * sigma, -1)
        return y

    def write_spectrum(self, filename, type='yaml'):
        """
        writes oscillator strenghts to file

        args:
            filename(str): Filename
            type(str): either txt or yaml
        """
        #TODO: csv
        spectrum = self.get_spectrum()
        if type == 'txt':
            self.__write_file(filename, spectrum)
        elif type == 'yaml':
            self.__write_file_yaml(filename, spectrum)

    def write_gaussiansmearedspectrum(self, filename, sigma, type='txt'):
        """
        writes smeared oscillator strenghts to file

        args:
            filename(str): Filename
            sigma(float): smearing of the spectrum
            type(str): either txt or yaml
        """
        #TODO csv
        spectrum = self.get_gaussiansmearedspectrum(sigma)
        if type == 'txt':
            self.__write_file(filename, spectrum)
        elif type == 'yaml':
            self.__write_file_yaml(filename, spectrum)

    def __write_file(self, filename, spectrum):
        """
        writes dict for any spectrum into txt file

        args:
            filename(str): Filename
            spectrum (dict): Includes nparray for 'Frequencies'
             and 'Intensities'

        """

        Freq = np.array(spectrum['Frequencies'].tolist())
        Intens = np.array(spectrum['Intensities'].tolist())
        file = open(filename, 'w')
        file.write('Frequency (cm-1) Oscillator Strengths \n')
        for i in range(len(Freq)):
            file.write('%s %s \n' % (Freq[i], Intens[i]))
        file.close()

    def __write_file_yaml(self, filename, spectrum):
        """
        writes dict for any spectrum into yaml file

        args:
                filename(str): Filename
                spectrum (dict): Includes nparray for 'Frequencies'
                and 'Intensities'

        """
        Freq = np.array(spectrum['Frequencies'].tolist())
        Intens = np.array(spectrum['Intensities'].tolist())
        file = open(filename, 'w')
        file.write('Frequency: \n')
        for i in range(len(Freq)):
            file.write('- %s \n' % (Freq[i]))
        file.write('Oscillator Strengths: \n')
        for i in range(len(Intens)):
            file.write('- %s \n' % (Intens[i]))
        file.close()

    def plot_spectrum(self, filename):
        """
        Plots frequencies in cm-1 and oscillator strengths
        args:
            filename(str): name of the file
        """
        spectrum = self.get_spectrum()
        plt.stem(spectrum['Frequencies'].tolist(), spectrum['Intensities'].tolist(), markerfmt=' ')
        plt.xlabel('Wave number (cm$^{-1}$)')
        plt.ylabel('Oscillator Strengths')
        plt.savefig(filename)
        plt.show()

    def plot_gaussiansmearedspectrum(self, filename, sigma):
        """
        Plots frequencies in cm-1 and smeared oscillator strengths
        args:
            filename(str): name of the file
            sigma(float): smearing
        """
        spectrum = self.get_gaussiansmearedspectrum(sigma)
        plt.plot(spectrum['Frequencies'].tolist(), spectrum['Intensities'].tolist())
        plt.xlabel('Wave number (cm$^{-1}$)')
        plt.ylabel('Oscillator Strengths')
        plt.savefig(filename)
        plt.show()

# def __plot_xydata(self,filename,spectrum):
#	spectrum=spectrum
