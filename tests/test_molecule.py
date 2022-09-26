import unittest
from src.Molecule import Molecule
import pandas as pd
from rdkit import Chem

scaffold = Chem.MolFromSmiles('O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]')
try:
    csv_file = pd.read_csv('r_group_decomp.csv')
except FileNotFoundError:
    pass


class TestMolecule(unittest.TestCase):
    """Tests for Molecule class"""
    def test_get_smile_string(self):
        mol = Molecule('Oc1ccc(C[*:1])cc1')
        true_case = 'Oc1ccc(C[*:1])cc1'
        self.assertEqual(mol.get_smile_string, true_case)

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
        message = {
            'PAINS': 'passing',
            'ZINC': 'passing',
            'BRENK': 'passing',
            'NIH': 'passing'
        }
        passing_mol = Molecule('CC')
        self.assertEqual(passing_mol.filter_properties(), message)

    def test_filter_fails(self):
        """Tests a molecule containing a thiocarbonyl fails the filter"""
        message = {
            'PAINS': ['thio_keto_het(2)'],
            'ZINC': 'passing',
            'BRENK': ['Thiocarbonyl_group'],
            'NIH': 'passing'}
        failing_mol = Molecule('O=C(NC(C1=C(C)C=C2N1C=CC=C2)=S)C3=CC=CC=C3')
        self.assertEqual(failing_mol.filter_properties(), message)

    def test_lipinski_no_violations(self):
        """Tests water correctly passes the Lipsinki rule of 5."""
        message = {
            'MW': True,
            'h_acc': True,
            'h_don': True,
            'logP': True
            }
        no_violation_molecule = Molecule('O')
        descriptors = no_violation_molecule.descriptors()
        self.assertEqual(no_violation_molecule.lipinski(descriptors), message)

    def test_lipinski_no_violations2(self):
        """Tests glucose (an edge case with h_don = 5) correctly
        passes the Lipsinki rule of 5."""
        message = {
                'MW': True,
                'h_acc': True,
                'h_don': True,
                'logP': True
            }
        no_violation_molecule2 = Molecule('C(C1C(C(C(C(O1)O)O)O)O)O')
        descriptors = no_violation_molecule2.descriptors()
        self.assertEqual(no_violation_molecule2.lipinski(descriptors), message)

    def test_lipinski_one_violation(self):
        """Tests that cholestrol fails one Lipinski rule and
        passses the rest"""
        message = {
            'MW': True,
            'h_acc': True,
            'h_don': True,
            'logP': False
            }
        one_violation_molecule = Molecule('CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3('
                                          'CCC(C4)O)C)C')
        descriptors = one_violation_molecule.descriptors()
        self.assertEqual(one_violation_molecule.lipinski(descriptors), message)

    def test_lipinski_fails(self):
        """Tests that Dextran fails the Lipinski rule of 5 and
        violates 3 of the 4 rules."""
        message = {
            'MW': False,
            'h_acc': False,
            'h_don': False,
            'logP': True
            }
        one_violation_molecule = Molecule('C(C1C(C(C(C(O1)OCC2C(C(C(C(O2)OCC('
                                          'C(C(C(C=O)O)O)O)O)O)O)O)O)O)O)O')
        descriptors = one_violation_molecule.descriptors()
        self.assertEqual(one_violation_molecule.lipinski(descriptors), message)

    def test_drawMoleculeAsByteStream(self):
        """Tests function by comparing bytestreams produced
        """
        test_smiles = ["O=C(O)C(NS(=O)(=O)c1ccc9cc1)8"
                       ".Oc1ccc(C8)cc1.c1ccc9cc1",
                       "O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]"]
        for i in test_smiles:
            passing_mol = Molecule(i)
            test_byte_stream = passing_mol.drawMoleculeAsByteStream(
                orient_with_scaffold=True
                )
            self.assertIsInstance(test_byte_stream, str)
            self.assertGreater(len(test_byte_stream), 1000)
