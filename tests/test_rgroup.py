import unittest
from src.Molecule import Molecule, R_group
import pandas as pd
from rdkit import Chem

scaffold = Chem.MolFromSmiles('O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]')
try:
    csv_file = pd.read_csv('r_group_decomp.csv')
except FileNotFoundError:
    pass


class TestRGroup(unittest.TestCase):
    """Tests for RGroup class"""
    def test_init_R_group(self):
        """Tests that when a non-valid R-group id is included,
        an Exception is raised"""
        with self.assertRaises(ValueError):
            R_group("A65")
        with self.assertRaises(ValueError):
            R_group("B65")
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
