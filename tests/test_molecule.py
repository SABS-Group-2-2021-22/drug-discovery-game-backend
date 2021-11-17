import unittest
from Molecule import Molecule
import pandas as pd
from rdkit import Chem


scaffold = Chem.MolFromSmiles('O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]')
try:
    csv_file = pd.read_csv('./drug-discovery-game-backend/r_group_decomp.csv')
except FileNotFoundError:
    csv_file = pd.read_csv('r_group_decomp.csv')


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

    # def test_lipinski_one_violation(self):
    #     "Tests that ... fails one Lipinski rule and passses the rest"

# one_violation_molecule = Molecule(...)
# descriptors = one_violation_molecule.descriptors()
# print(descriptors)
# print(one_violation_molecule.lipinski(descriptors))

    def test_lipinski_no_violations2(self):
        """Tests glucose (an edge case with h_don = 5) correctly
        passes the Lipsinki rule of 5."""
        message = 0, "passes"
        no_violation_molecule2 = Molecule('C(C1C(C(C(C(O1)O)O)O)O)O')
        descriptors = no_violation_molecule2.descriptors()
        self.assertEqual(no_violation_molecule2.lipinski(descriptors), message)

# if __name__ == '__main__':
#     unittest.main()


# Further tests needed - relevant tests from last year are below

# def test_apply_filter():
#     """Tests that a dataframe is filtered correctly"""
#     from game_scripts.filters import apply_filter
#     test_df = pd.DataFrame({'mols': [
#         'C=CCN1C(CCC(O)=O)=CC=C1C2=CC=C(F)C=C2',
#         'OC1=C(C(O)=O)C=CC=C1O',
#         'CC'
#         ]})
#     correct_df = pd.DataFrame({'mols': ['CC']})
#     correct_df.index = correct_df.index + 2
#     pd.testing.assert_frame_equal(correct_df, apply_filter(test_df, 'mols'))

# def test_molchoose_correct():
#     """Tests that the correct data is returned from a trial csv document"""

#     test = [[1,'A01B01','A01','B01','6.5','OC(=O)[C@H](Cc1ccc(O)cc1)NS(=O)
# (=O)c2ccc(cc2)c3ccccc3',7,12,'>98','Gen-5']]
#     test_frame = pd.DataFrame(test, columns=['Index','Tag','atag','btag',
# 'pIC50_MMP12','Smiles','A_SortMax','B_SortMax','Final QC Purity',
# 'Generation-No'])

#     pd.testing.assert_frame_equal(combine.MolChoose('A01', 'B01',
# DataSource='tests/test_input.csv'), test_frame)

# test_molchoose_correct()
