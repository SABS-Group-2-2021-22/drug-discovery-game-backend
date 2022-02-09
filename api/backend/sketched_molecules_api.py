from api.backend.backend_api import Backend
from src.SketchedMolecule import sketchedMolecule
from rdkit import Chem

class sketched_molecule(Backend):
    def __init__(self):
        

    def run_lipinski(mol):
        drug_mol = sketchedMolecule(mol)
        return drug_mol.lipinski()
    
    def 

        
