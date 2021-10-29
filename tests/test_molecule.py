import pytest
import sys
# sys.path.append('path')
# sys.path.insert(0, '/Molecule.py')
# import os.path
# cwd = os.path.dirname(__file__)  # get current working directory
from Molecule import Molecule
# import .r_group_decomp

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Crippen
from rdkit.Chem.FilterCatalog import *
from rdkit.Chem import AllChem


import numpy.testing as npt
# import pandas as pd
# from game_scripts import filters
# from rdkit import Chem
# from io import StringIO

#Make into class
scaffold = Chem.MolFromSmiles('O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]')
csv_file = pd.read_csv('r_group_decomp.csv')


"""Tests for getting descriptors for molecules and ranking them"""
def test_descriptors():
    """Can get all specified descriptors for single moiety and return as dict"""
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
    assert true_case == test_case

def test_filter_passes():
    """Tests a molecule correctly passes the filters"""
    message = "PASSED FILTERS"
    passing_mol = Molecule('CC')
    testmessage = passing_mol.filter_properties()
    assert passing_mol.filter_properties() == message


def test_filter_fails():    
    """Tests a molecule containing a thiocarbonyl fails the filter"""
    message = "FAIL FILTERS"
    failing_mol = Molecule('O=C(NC(C1=C(C)C=C2N1C=CC=C2)=S)C3=CC=CC=C3')
    assert failing_mol.filter_properties() == message

# Further tests needed

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

#     test = [[1,'A01B01','A01','B01','6.5','OC(=O)[C@H](Cc1ccc(O)cc1)NS(=O)(=O)c2ccc(cc2)c3ccccc3',7,12,'>98','Gen-5']]
#     test_frame = pd.DataFrame(test, columns=['Index','Tag','atag','btag','pIC50_MMP12','Smiles','A_SortMax','B_SortMax','Final QC Purity','Generation-No'])

#     pd.testing.assert_frame_equal(combine.MolChoose('A01', 'B01', DataSource='tests/test_input.csv'), test_frame)
 
# test_molchoose_correct()