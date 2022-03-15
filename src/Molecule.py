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
    csv_file = pd.read_csv('src/r_group_decomp.csv')
except FileNotFoundError:
    csv_file = pd.read_csv('r_group_decomp.csv')


class Molecule:
    """ A molecule. In particular, either the R1 or R2 group, or the scaffold
    and one or two groups.
    There are methods which tell you the properties of the molecule and if it
    passes the Lipsinki test
    """

    def __init__(self, mol_smiles):
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
        """Calculate Lipinski from the descriptor dictionary. Returns the
        number of rules broken and whether the molecule passes.

        :param desc_dict: molecule descriptor metrics
        :type desc_dict: dict
        :return: violations
        :rtype: dict
        """
        violations = {'MW': desc_dict['MW'] < 500.0,
                      'h_acc': desc_dict['h_acc'] <= 10,
                      'h_don': desc_dict['h_don'] <= 5,
                      'logP': desc_dict['logP'] < 5}
        return violations

    def drawMoleculeAsByteStream(self, orient_with_scaffold=False, size=None):
        """Returns png image of molecule as bytestream

        :return: base64 png image bytestream
        :rtype: String
        """

        # Determine R1, R2
        if '[*:1]' in self.__mol_smiles:
            R1 = True
        else:
            R1 = False

        if '[*:2]' in self.__mol_smiles:
            R2 = True
        else:
            R2 = False

        # drawn_mol = Chem.MolFromSmiles(self.__mol_smiles)
        smiles = self.__mol_smiles
        colours = [(0.8, 0.0, 0.8), (0.8, 0.8, 0)]

        if R1 is True:
            # Replace the R1 group tag with Xe
            smiles = smiles.replace('*:1', 'Xe')

        if R2 is True:
            # Replace the R2 group tag with Ne
            smiles = smiles.replace('*:2', 'Ne')

        # if R1 != True and R2 != True:
        #     raise ValueError('Neither R1 nor R2 is True')

        mol = Chem.MolFromSmiles(smiles)

        patt1 = Chem.MolFromSmarts('[Xe]')
        patt2 = Chem.MolFromSmarts('[Ne]')

        # This will be populated with the (former) locations of the R1 and R2
        # group tags, and corresponding colours
        atom_cols = {}

        if R1:
            # Finds the location of the R1 group
            hit_ats1 = [mol.GetSubstructMatch(patt1)[0]]
            # Assigns colour[0] to this R1 group for labelling
            atom_cols[hit_ats1[0]] = colours[0]

        else:
            hit_ats1 = []

        if R2:                      # Similarly to R1
            hit_ats2 = [mol.GetSubstructMatch(patt2)[0]]
            atom_cols[hit_ats2[0]] = colours[1]

        else:
            hit_ats2 = []

        hit_ats = hit_ats1 + hit_ats2  # List of R group locations

        # Remove the Xe and Ne atoms marking the R group locations
        # and replace with C atoms
        smi_draw = smiles.replace('[Xe]', 'C')
        smi_draw = smi_draw.replace('[Ne]', 'C')

        # Ensuring the case A11 = '[H][*:1]' diplays as H not CH4
        if smi_draw == '[H]C':
            smi_draw = '[H]'
        else:
            pass

        mol_draw = Chem.MolFromSmiles(smi_draw)

        if orient_with_scaffold is True:
            AllChem.Compute2DCoords(orient_scaffold)
            AllChem.Compute2DCoords(mol_draw)
            # TODO: test if _ assignment is needed or if fn call wo./
            # assignment is sufficient
            _ = AllChem.GenerateDepictionMatching2DStructure(mol_draw,
                                                             orient_scaffold)

        if size is not None:
            d = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
        else:
            size = d = rdMolDraw2D.MolDraw2DCairo(500, 500)

        # Add the highlighting
        d.DrawMolecule(mol_draw, highlightAtoms=hit_ats,
                       highlightAtomColors=atom_cols)

        # Add the highlighting
        d.FinishDrawing()
        b = io.BytesIO(d.GetDrawingText())
        imgByteArray = base64.b64encode(b.read()).decode("utf-8")
        return imgByteArray

    def filter_properties(self):
        """See whether molecule passes or fails FILTERS
        """
        pains_params = FilterCatalogParams()
        pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
        pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
        pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)

        zinc_params = FilterCatalogParams()
        zinc_params.AddCatalog(FilterCatalogParams.FilterCatalogs.ZINC)

        brenk_params = FilterCatalogParams()
        brenk_params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)

        nih_params = FilterCatalogParams()
        nih_params.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)

        params = [pains_params, zinc_params, brenk_params, nih_params]
        filt_names = ['PAINS', 'ZINC', 'BRENK', 'NIH']
        mol = Chem.MolFromSmiles(self.__mol_smiles)
        all_warnings = {}

        for param, name in zip(params, filt_names):
            catalog = FilterCatalog(param)
            entries = catalog.GetMatches(mol)
            warnings = []
            for e in entries:
                warnings.append(e.GetDescription())
            if len(warnings) != 0:
                all_warnings[name] = warnings
            else:
                all_warnings[name] = 'passing'
        return all_warnings


class R_group(Molecule):
    """Name of R group is of the form 'Axy' or 'Bxy' e.g. A01 etc.
    Number corresponds to whether it is an R1 or R2 group
    """

    def __init__(self, name, number=None):
        self.name = name
        if number is None:
            number_key = int(name[1]+name[2])
            if number_key <= 50:
                if name[0] == 'A':
                    number = 1
                elif name[0] == 'B':
                    number = 2
                else:
                    raise ValueError('R group ID not valid')
            else:
                raise ValueError('R group ID not valid')
        self.number = number
        mol_smiles = self.extract_smilefromcsv()
        super().__init__(mol_smiles)

    def extract_smilefromcsv(self):
        """Extracts the SMILE for the R group
        """
        if self.number == 1:
            rgroup_smiles = csv_file[csv_file['atag'] ==
                                     self.name]['R1'].iloc[0]
        if self.number == 2:
            rgroup_smiles = csv_file[csv_file['btag'] ==
                                     self.name]['R2'].iloc[0]
        return(rgroup_smiles)

    def add_r_group(self, base_molecule):
        new_mol = base_molecule.get_smile_string + '.' + self.get_smile_string
        connection_site = str(7+self.number)
        new_mol = new_mol.replace(f'[*:{self.number}]', connection_site)
        new_mol = new_mol.replace('('+connection_site+')', connection_site)
        return Molecule(new_mol)


class FinalMolecule(Molecule):
    """Final molecule with scaffold and two R groups.
    """

    def __init__(self, rgroup1, rgroup2):
        # Name of R groups should be in the form 'Axy' or 'Bxy' e.g. A01 etc.
        self.rgroup1 = rgroup1
        self.rgroup2 = rgroup2
        mol_smiles = self.build_final_smiles()
        super().__init__(mol_smiles)

    def build_final_smiles(self):
        r_group_mol_1 = R_group(self.rgroup1, 1)
        r_group_mol_2 = R_group(self.rgroup2, 2)
        intermediate_mol = r_group_mol_1.add_r_group(
            Molecule('O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]')
            )
        final_mol = r_group_mol_2.add_r_group(intermediate_mol)
        final_mol_smiles = final_mol.get_smile_string
        return(final_mol_smiles)

    def drug_properties(self):
        """Selects properties of the final drug from the data.
        """
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
