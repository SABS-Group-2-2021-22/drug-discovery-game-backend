import unittest
from src.Molecule import Molecule, R_group, FinalMolecule
import pandas as pd
from rdkit import Chem

scaffold = Chem.MolFromSmiles('O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]')
try:
    csv_file = pd.read_csv('src/r_group_decomp.csv')
except FileNotFoundError:
    csv_file = pd.read_csv('tests/r_group_decomp.csv')


class TestMolecule(unittest.TestCase):
    """Tests for getting descriptors for molecules and ranking them"""
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
        passing_mol = Molecule("O=C(O)C(NS(=O)(=O)c1ccc9cc1)8"
                               ".Oc1ccc(C8)cc1.c1ccc9cc1")
        test_byte_stream = passing_mol.drawMoleculeAsByteStream(
            orient_with_scaffold=True
            )
        f = open("tests/test_draw_molecule_bytes.txt", "r")
        true_byte_stream = f.read()
        print(test_byte_stream)
        self.assertEqual(test_byte_stream, true_byte_stream)

    def test_init_R_group(self):
        """Tests that when a non-valid R-group id is included,
        an Exception is raised"""
        with self.assertRaises(ValueError):
            R_group("A65")
        with self.assertRaises(ValueError):
            R_group("C01")

    def test_extract_smilefromcsv(self):
        """Tests R_group function by comparing values extracted from the csv
        """
        passing_mol = R_group("A01")
        test_smiles = passing_mol.extract_smilefromcsv()
        true_smiles = "Oc1ccc(C[*:1])cc1"
        self.assertEqual(test_smiles, true_smiles)

    def test_add_r_group(self):
        """Test R_group function by comparing smiles string generated when
        adding R group to a molecule
        """
        r_group_mol_1 = R_group("A01")
        intermediate_mol = r_group_mol_1.add_r_group(
            Molecule('O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]')
            )
        true_smiles = 'O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)8.Oc1ccc(C8)cc1'
        test_smiles = intermediate_mol.get_smile_string
        self.assertEqual(test_smiles, true_smiles)

    def test_build_final_smiles(self):
        """Tests FinalMolecule function by comparing smiles strings
        """
        passing_mol = FinalMolecule("A01", "B01")
        test_smiles = passing_mol.build_final_smiles()
        true_smiles = "O=C(O)C(NS(=O)(=O)c1ccc9cc1)8.Oc1ccc(C8)cc1.c1ccc9cc1"
        self.assertEqual(test_smiles, true_smiles)

    def test_drug_properties(self):
        """Tests that correct drug properties are returned for a molecule
        """
        passing_mol = FinalMolecule("A04", "B08")
        test_dict = passing_mol.drug_properties()
        true_dict = {
            'pic50': '5.2',
            'clearance_mouse': 'medium (5.6-30.5)',
            'clearance_human': 'low (< 12)',
            'logd': 1.51,
            'pampa': 'med2high'
            }
        self.assertEqual(test_dict, true_dict)
