from ctypes.util import find_library
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

    def test_indices(self):
        """Tests that ligand lipophilicity efficiency index (LLE), ligand
        efficiency index (LEI) and ligand efficiency (LE) are calculated
        accurately.
        """
        passing_mol = FinalMolecule("A02", "B07")
        desc_dict = {'logP': 5, 'HA': 10}
        drug_dict = {'pic50': 10}
        test_dict = passing_mol.indices(desc_dict, drug_dict)
        true_dict = {'LE': 1.37, 'LEI': 1.0, 'LLE': 5.0}
        self.assertEqual(test_dict, true_dict)

    def test_indices_no_pIC50(self):
        """Tests that ligand lipophilicity efficiency index (LLE), ligand
        efficiency index (LEI) and ligand efficiency (LE) are returned None
        if pIC50 has no value
        """
        passing_mol = FinalMolecule("A10", "B01")
        desc_dict = {'logP': 5, 'HA': 10}
        drug_dict = {'pic50': 'Not Made'}
        test_dict = passing_mol.indices(desc_dict, drug_dict)
        true_dict = {
            'LLE': None,
            'LEI': None,
            'LE': None
        }
        self.assertEqual(test_dict, true_dict)

    def test_astrazeneca(self):
        """Tests that LLE>5 calculatator is true
        """
        passing_mol = FinalMolecule("A09", "B22")
        indices = {'LE': 1.37, 'LEI': 1.0, 'LLE': 5.0}
        astrazeneca_rule = passing_mol.astrazeneca(indices)
        self.assertEqual(astrazeneca_rule, False)

    def test_astrazeneca_no_pIC50(self):
        """Tests that LLE>5 calculatator returns None if has no pIC50 value
        """
        passing_mol = FinalMolecule("A10", "B01")
        indices = {'LE': None, 'LEI': None, 'LLE': None}
        astrazeneca_rule = passing_mol.astrazeneca(indices)
        self.assertEqual(astrazeneca_rule, None)

    def test_pfizer(self):
        """Tests that clogP < 3 and TPSA > 75 is calculated correctly
        """
        passing_mol = FinalMolecule("A03", "B02")
        desc_dict_1 = {'logP': 1, 'TPSA': 80}
        pfizer_rule = passing_mol.pfizer(desc_dict_1)
        self.assertEqual(pfizer_rule, {'clogP': True, 'TPSA': True})
        desc_dict_2 = {'logP': 5, 'TPSA': 80}
        pfizer_rule = passing_mol.pfizer(desc_dict_2)
        self.assertEqual(pfizer_rule, {'clogP': False, 'TPSA': True})

    def test_gsk(self):
        """Tests that MW < 400, cLogP < 4 is calculated correctly
        """
        passing_mol = FinalMolecule("A10", "B40")
        desc_dict_1 = {'MW': 100, 'logP': 1}
        gsk_rule = passing_mol.gsk(desc_dict_1)
        self.assertEqual(gsk_rule, {'clogP': True, 'MW': True})
        desc_dict_2 = {'MW': 500, 'logP': 6}
        gsk_rule = passing_mol.gsk(desc_dict_2)
        self.assertEqual(gsk_rule, {'clogP': False, 'MW': False})
