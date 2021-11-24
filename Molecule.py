import pandas as pd
from rdkit import Chem
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Crippen
from rdkit.Chem import AllChem
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

import io
import base64

# Build scaffold and read in csv
scaffold = Chem.MolFromSmiles('O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]')
orient_scaffold = Chem.MolFromSmiles('O=C(O)C(NS(=O)(=O)c1ccccc1)')
try:
    csv_file = pd.read_csv('./drug-discovery-game-backend/r_group_decomp.csv')
except FileNotFoundError:
    csv_file = pd.read_csv('r_group_decomp.csv')


class Molecule:
    """ A molecule. In particular, either the R1 or R2 group, or the scaffold
    and one or two groups.
    There are methods which tell you the properties of the molecule and if it
    passes the Lipsinki test
    """

    def __init__(self, mol_smiles=None):
        """Constructor for Molecule class. Initialises Molecule instance from
        smile string.

        :param mol_smiles: smile string of molecule
        :type mol_smiles: String
        """
        self.__mol_smiles = mol_smiles

    @property
    def get_smile_string(self):
        """Returns molecule's smile string

        :return: smile string of molecule
        :rtype: String
        """
        return self.__mol_smiles

    def descriptors(self):
        """Calculate molecule descriptor metrics as dict:
        | mol - smile string
        | MW - molecular weight
        | logP - logP
        | TPSA - topological polar surface area
        | HA - heavy atom count
        | h_acc - H acceptor count
        | h_don - H donator count
        | rings - ring count

        :return: molecule descriptor metrics
        :rtype: dict
        """
        mol = Chem.MolFromSmiles(self.__mol_smiles)
        mw = Descriptors.ExactMolWt(mol)
        log_p = Crippen.MolLogP(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol)  # topological polar surface area
        ha = Lipinski.HeavyAtomCount(mol)  # heavy atom count
        h_acceptors = Lipinski.NumHAcceptors(mol)
        h_donors = Lipinski.NumHDonors(mol)
        rings = Lipinski.RingCount(mol)
        desc_dict = {'mol': self.__mol_smiles,
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

        :param desc_dict: molecule descriptor metrics
        :type desc_dict: dict
        :return: violations, result
        | Number of violations and 'fails' or passes'
        :rtype: int, String
        """
        violations = [desc_dict['MW'] >= 500.0,
                      desc_dict['h_acc'] > 10,
                      desc_dict['h_don'] > 5,
                      desc_dict['logP'] >= 5
                      ].count(True)
        if violations > 1:
            result = 'fails'
        else:
            result = 'passes'
        return violations, result

# csv_file = pd.read_csv('r_group_decomp.csv')
    def draw_molecule(self, drawn_file_name, orient_with_scaffold):
        """Draws the molecule.

        :param drawn_file_name: filename to save drawn molecule with
        only if the molecule contains the scaffold.
        :type orient_with_scaffold: bool
        """
        drawn_mol = Chem.MolFromSmiles(self.__mol_smiles)
        # Align molecule with scaffold if the molecule contains the scaffold.
        if orient_with_scaffold is True:
            AllChem.Compute2DCoords(orient_scaffold)
            # TODO: test if _ assignment is needed or if fn call wo./
            # assignment is sufficient
            _ = AllChem.GenerateDepictionMatching2DStructure(drawn_mol,
                                                             orient_scaffold)
        d = rdMolDraw2D.MolDraw2DCairo(250, 200)
        d.FinishDrawing()
        d.WriteDrawingText(f'{drawn_file_name}.png')

    def drawMoleculeAsByteStream(self, orient_with_scaffold=False):
        """Returns png image of molecule as bytestream

        :return: base64 png image bytestream
        :rtype: String
        """
        drawn_mol = Chem.MolFromSmiles(self.__mol_smiles)
        if orient_with_scaffold is True:
            AllChem.Compute2DCoords(orient_scaffold)
            AllChem.Compute2DCoords(drawn_mol)
            # TODO: test if _ assignment is needed or if fn call wo./
            # assignment is sufficient
            _ = AllChem.GenerateDepictionMatching2DStructure(drawn_mol,
                                                             orient_scaffold)
        img = Chem.Draw.MolToImage(drawn_mol)
        imgByteArray = io.BytesIO()
        img.save(imgByteArray, format='png')
        imgByteArray = imgByteArray.getvalue()
        imgByteArray = base64.b64encode(imgByteArray).decode("utf-8")
        return imgByteArray

    def filter_properties(self):
        """See whether molecule passes or fails FILTERS"""
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.
                          FilterCatalogs.PAINS_A)
        params.AddCatalog(FilterCatalogParams.
                          FilterCatalogs.PAINS_B)
        params.AddCatalog(FilterCatalogParams.
                          FilterCatalogs.PAINS_C)
        params.AddCatalog(FilterCatalogParams.
                          FilterCatalogs.ZINC)
        params.AddCatalog(FilterCatalogParams.
                          FilterCatalogs.BRENK)
        params.AddCatalog(FilterCatalogParams.
                          FilterCatalogs.NIH)
        catalog = FilterCatalog(params)
        mol = Chem.MolFromSmiles(self.__mol_smiles)
        if catalog.HasMatch(mol):
            return "FAIL FILTERS"
        else:
            return "PASSED FILTERS"


class R_group(Molecule):
    """Name of R group is of the form 'Axy' or 'Bxy' e.g. A01 etc.
    Number corresponds to whether it is an R1 or R2 group"""
    def __init__(self, name, number):
        self.name = name
        self.number = number
        mol_smiles = self.extract_smilefromcsv()
        super().__init__(mol_smiles)

    def extract_smilefromcsv(self):
        """Extracts the SMILE for the R group """

        if self.number == 1:
            rgroup_smiles = csv_file[csv_file['atag'] ==
                                     self.name]['R1'].iloc[0]
        if self.number == 2:
            rgroup_smiles = csv_file[csv_file['btag'] ==
                                     self.name]['R2'].iloc[0]
        return(rgroup_smiles)


class Scaffold_and_Rgroups(Molecule):
    # Can I also make this have R-group???
    """Add a R group to the old molecule (either the scaffold with R1 attached
    or just the scaffold)."""
    def ___init___(self, old_molecule, rgroup, number):

        self.old_molecule = old_molecule
        self.rgroup = rgroup
        self.number = number
        self.new_molecule = None
        self.add_r_group()

    # def add_r_group(self, new_mol):
    def add_r_group(self):
        """ Add R group to molecule"""
        new_mol = Chem.MolToSmiles(self.old_molecule) + '.' + self.rgroup
        new_mol = new_mol.replace(f'[*:{self.number}]', '9')
        new_mol = new_mol.replace('(9)', '9')
        self.new_molecule = new_mol


class FinalMolecule(Molecule):
    """Final molecule with scaffold and two R groups."""
    def __init__(self, rgroup1, rgroup2):
        # Name of R groups should be in the form 'Axy' or 'Bxy' e.g. A01 etc.
        super().__init__()
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
        drug_property_dict = {}
        for d in drug_properties:
            value = csv_file[csv_file['atag'] == self.rgroup1]
            value_two = value[csv_file['btag'] == self.rgroup2][d].iloc[0]
            drug_property_dict[d] = value_two
        return drug_property_dict
