import unittest
from src.Molecule import Molecule
import pandas as pd
from rdkit import Chem

scaffold = Chem.MolFromSmiles('O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]')
try:
    csv_file = pd.read_csv('src/r_group_decomp.csv')
except FileNotFoundError:
    csv_file = pd.read_csv('tests/r_group_decomp.csv')


class TestMolecule(unittest.TestCase):
    """Tests for getting descriptors for molecules and ranking them"""
    def test_descriptors(self):
        """Can get all specified descriptors for single moiety and
        return as dict"""
        mol = Molecule('Oc1ccc(C[*:1])cc1')
        true_case = {'mol': 'Oc1ccc(C[*:1])cc1',
                     'MW': 107.0497,
                     'logP': 1.4415,
                     'TPSA': 20.23,
                     'HA': 8,
                     'h_acc': 1,
                     'h_don': 1,
                     'rings': 1
                     }
        test_case = mol.descriptors()
        for key in true_case.keys():
            self.assertAlmostEqual(true_case[key], test_case[key], places=4)

    def test_filter_passes(self):
        """Tests a molecule correctly passes the filters"""
        message = "PASSED FILTERS"
        passing_mol = Molecule('CC')
        self.assertEqual(passing_mol.filter_properties(), message)

    def test_filter_fails(self):
        """Tests a molecule containing a thiocarbonyl fails the filter"""
        message = "FAIL FILTERS"
        failing_mol = Molecule('O=C(NC(C1=C(C)C=C2N1C=CC=C2)=S)C3=CC=CC=C3')
        self.assertEqual(failing_mol.filter_properties(), message)

    def test_lipinski_no_violations(self):
        """Tests water correctly passes the Lipsinki rule of 5."""
        message = 0, "passes"
        no_violation_molecule = Molecule('O')
        descriptors = no_violation_molecule.descriptors()
        self.assertEqual(no_violation_molecule.lipinski(descriptors), message)

    def test_lipinski_no_violations2(self):
        """Tests glucose (an edge case with h_don = 5) correctly
        passes the Lipsinki rule of 5."""
        message = 0, "passes"
        no_violation_molecule2 = Molecule('C(C1C(C(C(C(O1)O)O)O)O)O')
        descriptors = no_violation_molecule2.descriptors()
        self.assertEqual(no_violation_molecule2.lipinski(descriptors), message)

    def test_lipinski_one_violation(self):
        """Tests that cholestrol fails one Lipinski rule and
        passses the rest"""
        message = 1, "passes"
        one_violation_molecule = Molecule('CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3('
                                          'CCC(C4)O)C)C')
        descriptors = one_violation_molecule.descriptors()
        self.assertEqual(one_violation_molecule.lipinski(descriptors), message)

    def test_lipinski_fails(self):
        """Tests that Dextran fails the Lipinski rule of 5 and
        violates 3 of the 4 rules."""
        message = 3, "fails"
        one_violation_molecule = Molecule('C(C1C(C(C(C(O1)OCC2C(C(C(C(O2)OCC('
                                          'C(C(C(C=O)O)O)O)O)O)O)O)O)O)O)O')
        descriptors = one_violation_molecule.descriptors()
        self.assertEqual(one_violation_molecule.lipinski(descriptors), message)

# if __name__ == '__main__':
#     unittest.main()


# Further tests are needed - relevant tests from last year that could be
# adapted are test_apply_filter(),  test_molchoose_correct()
