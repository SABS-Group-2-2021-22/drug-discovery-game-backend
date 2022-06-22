import unittest
from src.Molecule import FinalMolecule
import pandas as pd
from rdkit import Chem

scaffold = Chem.MolFromSmiles('O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]')
try:
    csv_file = pd.read_csv('r_group_decomp.csv')
except FileNotFoundError:
    pass


class TestFinalMolecule(unittest.TestCase):
    """Tests for FinalMolecule class"""
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

    # def test_indices(self):
    #     """Tests that ligand lipophilicity efficiency index (LLE), ligand
    #     efficiency index (LEI) and ligand efficiency (LE) are calculated
    #     accurately.
    #     """
    #     passing_mol = FinalMolecule("A01", "B01")
    #     desc_dict = passing_mol.descriptors()
    #     print(desc_dict)
    #     drug_dict = passing_mol.drug_properties()
    #     print(drug_dict)
    #     test_dict = passing_mol.indices(desc_dict, drug_dict)
    #     true_dict = {'LLE': 3.466599999999999,
    #                  'LEI': 0.23214285714285715,
    #                  'LE': 0.3180357142857143}
    #     self.assertEqual(test_dict, true_dict)

    # def test_indices_no_pIC50(self):
    #     """Tests that ligand lipophilicity efficiency index (LLE), ligand
    #     efficiency index (LEI) and ligand efficiency (LE) are returned None
    #     if pIC50 has no value
    #     """
    #     passing_mol = FinalMolecule("A10", "B01")
    #     desc_dict = passing_mol.descriptors()
    #     drug_dict = passing_mol.drug_properties()
    #     test_dict = passing_mol.indices(desc_dict, drug_dict)
    #     true_dict = {
    #         'LLE': None,
    #         'LEI': None,
    #         'LE': None
    #     }
    #     self.assertEqual(test_dict, true_dict)

    def test_astrazeneca(self):
        """Tests that LLE>5 calculatator is true
        """
