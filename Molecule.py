import pandas as pd
from rdkit import Chem
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Crippen
# from rdkit.Chem import FilterCatalog
# from rdkit.Chem.FilterCatalog import FilterCatalogParams
from rdkit.Chem import AllChem

#can't delete the following!!! Steph: last year they couldn't 
# get rid of the warning message
from rdkit.Chem.FilterCatalog import *

# Didn't work:
# from rdkit.Chem import FilterCatalog
# from rdkit.Chem.FilterCatalog import FilterCatalogParams

# Build scaffold and read in csv
scaffold = Chem.MolFromSmiles('O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]')
csv_file = pd.read_csv('r_group_decomp.csv')


class Molecule:
    """A molecule. In particular, either the R1 or R2 group, or the scaffold
    and one or two groups.
     There are methods which tell you the properties of the molecule and if it
     passes the Lipsinki test."""

    def __init__(self, mol_smiles):
        self.mol_smiles = mol_smiles

    def descriptors(self):
        """Calculate descriptors"""
        mol = Chem.MolFromSmiles(self.mol_smiles)
        mw = round(Descriptors.ExactMolWt(mol), 4) #round so easier to test
        log_p = Crippen.MolLogP(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol)  # topological polar surface area
        ha = Lipinski.HeavyAtomCount(mol)  # heavy atom count
        h_acceptors = Lipinski.NumHAcceptors(mol)
        h_donors = Lipinski.NumHDonors(mol)
        rings = Lipinski.RingCount(mol)

        desc_dict = {'mol': self.mol_smiles,
                     'MW': mw,
                     'logP': log_p,
                     'TPSA': tpsa,
                     'HA': ha,
                     'h_acc': h_acceptors,
                     'h_don': h_donors,
                     'rings': rings
                     }
        return desc_dict

    def lipinski(self, desc_dict):
        """Calculate Lipinski from the descriptor dictionary.
        Return the number of rules broken and whether the molecule passes.
        """
        violations = [desc_dict['MW'] >= 500.0,
                      desc_dict['h_acc'] >= 10,
                      desc_dict['h_don'] >= 5,
                      desc_dict['logP'] >= 5
                      ].count(True)
        if violations > 1:
            result = 'fails'
        else:
            result = 'passes'
        return violations, result

    def draw_molecule(self, drawn_file_name, orient_with_scaffold):
        """Draws the molecule. The variable 'orient_with_scaffold' can take
        the values True and False. It should take the value True if and
        only if the molecule contains the scaffold"""
        drawn_mol = Chem.MolFromSmiles(self.mol_smiles)
        # Align molecule with scaffold if the molecule contains the scaffold.
        if orient_with_scaffold is True:
            AllChem.Compute2DCoords(scaffold)
            _ = AllChem.GenerateDepictionMatching2DStructure(drawn_mol,
                                                             scaffold)
        d = rdMolDraw2D.MolDraw2DCairo(250, 200)
        d.drawOptions().addStereoAnnotation = True
        d.drawOptions().clearBackground = False
        d.DrawMolecule(drawn_mol)
        d.FinishDrawing()
        d.WriteDrawingText(f'{drawn_file_name}.png')

    def filter_properties(self, drawn_mol_final):
        """See whether molecule passes or fails FILTERS"""
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.ZINC)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
        catalog = FilterCatalog(params)
        mol = Chem.MolFromSmiles(self.mol_smiles) 
        if catalog.HasMatch(mol):
            return "FAIL FILTERS"
        else:
            return "PASSED FILTERS"


class R_group(Molecule):
    """Name of R group is of the form 'Axy' or 'Bxy' e.g. A01 etc.
    Number corresponds to whether it is an R1 or R2 group"""
    def ___init___(self, name, number):
        self.name = name
        self.number = number

    def extract_smilefromcsv(self):
        """Extracts the SMILE for the R group """
        if self.number == 1:
            rgroup_smiles = csv_file[csv_file['atag'] == self.name]['R1'][0]
        if self.number == 2:
            rgroup_smiles = csv_file[csv_file['btag'] == self.name]['R2'][0]
        print(rgroup_smiles)


class Scaffold_and_Rgroups(Molecule):
    # Can I also make this have R-group???
    """Add a R group to the old molecule (either the scaffold with R1 attached
    or just the scaffold)."""
    def ___init___(self, old_molecule, rgroup, number):
        self.old_molecule = old_molecule
        self.rgroup = rgroup
        self.number = number
        self.new_molecule = None

    def add_r_group(self, new_mol):
        """ Add R group to molecule"""
        new_mol = Chem.MolToSmiles(self.old_molecule) + '.' + self.rgroup
        new_mol = new_mol.replace(f'[*:{self.number}]', '9')
        new_mol = new_mol.replace('(9)', '9')
        self.new_molecule = new_mol


class FinalMolecule(Molecule):
    """Final molecule with scaffold and two R groups."""
    def ___init___(self, rgroup1, rgroup2):
        # Name of R groups should be in the form 'Axy' or 'Bxy' e.g. A01 etc.
        self.rgroup1 = rgroup1
        self.rgroup2 = rgroup2

    def drug_properties(self):
        """Selects properties of the final drug from the data."""
        drug_properties = [
                            'pic50',
                            'clearance_mouse',
                            'clearance_human',
                            'logd',
                            'pampa'
                            ]
        for d in drug_properties:
            value = csv_file[csv_file['atag'] == self.rgroup1]
            value_two = value[csv_file['btag'] == self.rgroup2][d][0]
        print(d, value_two)
