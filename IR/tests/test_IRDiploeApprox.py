from IR.IRDipoleApprox import IR
import unittest
import numpy as np
import os

path_here = os.path.dirname(__file__)


class IRTest(unittest.TestCase):
    def setUp(self):
        self.IR_calc = IR(PoscarName=os.path.join(path_here, 'POSCAR'),
                          ForceConstants=False,
                          ForceFileName=os.path.join(path_here, 'FORCE_SETS'),
                          supercell=[[3, 0, 0], [0, 3, 0], [0, 0, 4]],
                          primitive=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], BornFileName=os.path.join(path_here, 'BORN'))
        self.IR_calc_masses = IR(PoscarName=os.path.join(path_here, 'POSCAR'),
                                 ForceConstants=False,
                                 ForceFileName=os.path.join(path_here, 'FORCE_SETS'),
                                 supercell=[[3, 0, 0], [0, 3, 0], [0, 0, 4]],
                                 primitive=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                 masses=[12.010700, 12.010700, 15.999400, 15.999400,
                                         14.006700, 14.006700, 14.006700, 14.006700, 2,
                                         2,
                                         2, 2, 2, 2, 2, 2],
                                 BornFileName=os.path.join(path_here, 'BORN'))

        self.IR_calc2 = IR(PoscarName=os.path.join(path_here, 'POSCAR.NaCl'),
                           ForceConstants=False,
                           ForceFileName=os.path.join(path_here, 'FORCE_SETS.NaCl'),
                           supercell=[[2, 0, 0], [0, 2, 0], [0, 0, 2]], nac=True,
                           BornFileName=os.path.join(path_here, 'BORN.NaCl'),
                           primitive=[[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])
        self.IRFC = IR(PoscarName=os.path.join(path_here, 'POSCAR_Methanol'),
                       ForceConstants=True,
                       ForceFileName=os.path.join(path_here, 'FORCE_CONSTANTS_Methanol'),
                       BornFileName=os.path.join(path_here, 'BORN.Methanol'),
                       supercell=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], nac=False)

    def test_attributes(self):
        # test calculation of frequencies
        self.assertAlmostEqual(self.IR_calc._frequencies[47], 3490.6434922723, places=1)
        # test calculation of eigenvectors
        self.assertAlmostEqual(abs(self.IR_calc._EigFormat[15, 47, 0]), 0.00084433323436)
        self.assertAlmostEqual(abs(self.IR_calc._EigFormat[15, 47, 1]), 0.00084433323436)
        self.assertAlmostEqual(abs(self.IR_calc._EigFormat[15, 47, 2]), 0.37170414232138)
        # check if sign of eigenvectors is consistent!!
        self.assertEqual(np.sign(self.IR_calc._EigFormat[14, 47, 2]),
                         np.sign(self.IR_calc._EigFormat[15, 47, 0]))
        self.assertEqual(np.sign(self.IR_calc._EigFormat[14, 47, 2]),
                         np.sign(self.IR_calc._EigFormat[15, 47, 1]))
        self.assertEqual(np.sign(self.IR_calc._EigFormat[14, 47, 2]),
                         np.sign(self.IR_calc._EigFormat[15, 47, 2]))

        # test intensities
        self.assertAlmostEqual(self.IR_calc.get_intensities()[-1], 2.403735277960254)

        # TODO: test NAC
        self.assertAlmostEqual(self.IR_calc2._frequencies[-1], 153.7212069157, places=2)

        # TODO: set masses externally [e.g., use D mass]
        self.assertAlmostEqual(self.IR_calc_masses._frequencies[47], 2598.2875793589, places=1)
        # test calculation of eigenvectors
        self.assertAlmostEqual(abs(self.IR_calc_masses._EigFormat[15, 47, 0]), 0.00378948635566)
        self.assertAlmostEqual(abs(self.IR_calc_masses._EigFormat[15, 47, 1]), 0.00378948635566)
        self.assertAlmostEqual(abs(self.IR_calc_masses._EigFormat[15, 47, 2]), 0.33223420830758)
        # check if sign of eigenvectors is consistent
        self.assertEqual(np.sign(self.IR_calc_masses._EigFormat[14, 47, 2]),
                         np.sign(self.IR_calc_masses._EigFormat[15, 47, 0]))
        self.assertEqual(np.sign(self.IR_calc_masses._EigFormat[14, 47, 2]),
                         np.sign(self.IR_calc_masses._EigFormat[15, 47, 1]))
        self.assertEqual(np.sign(self.IR_calc_masses._EigFormat[14, 47, 2]),
                         np.sign(self.IR_calc_masses._EigFormat[15, 47, 2]))

        # test intensities
        self.assertAlmostEqual(self.IR_calc_masses.get_intensities()[-1], 1.1658569108600518)

        # start from FORCE constants instead
        self.assertAlmostEqual(self.IRFC._frequencies[-1], 3741.4132865293, places=1)

        # test intensities
        self.assertAlmostEqual(self.IRFC.get_intensities()[-1], 0.02218877246084381)

    def test_printing_files(self):
        # TODO: implement test of writing of files
        pass

    def test_smearing(self):
        # TODO: test smearing

        pass


if __name__ == '__main__':
    unittest.main()
