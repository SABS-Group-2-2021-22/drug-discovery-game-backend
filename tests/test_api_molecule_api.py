import unittest
from src.Molecule import Molecule
import pandas as pd
from rdkit import Chem

from flask import Flask

from api.backend.molecule_api import *

app = Flask(__name__)       # emtpy app to run tests within flask app context

class TestMoleculeAPI(unittest.TestCase):

    def test_all_comparison_texts(self):
        with app.app_context():
            for i in range(1, 26):
                if i < 10:
                    mol_i = 'A0'+ str(i)
                else: 
                    mol_i = 'A' + str(i) 
                for j in range(1, 26):
                    if j < 10:
                        mol_j = 'B0'+ str(j)
                    else: 
                        mol_j = 'B' + str(j) 
                
                    comparison_txt((mol_i, mol_j))

