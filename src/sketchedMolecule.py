from rdkit import Chem
from rdkit.Chem import Lipinski,\
    Descriptors, rdMolDescriptors, Crippen, rdDepictor
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

import io
import base64


class sketchedMolecule:
    def __init__(self, mol_block):
        self.mol_block = mol_block
        self.mol = Chem.rdmolfiles.MolFromMolBlock(self.mol_block)
        old_mol = Chem.rdmolfiles.MolFromMolBlock(self.mol_block)
        rdDepictor.Compute2DCoords(self.mol, clearConfs=True)
        Chem.rdMolAlign.AlignMol(self.mol, old_mol)

    @property
    def smiles(self):
        '''Returns smile string for the molecule

        :return: smiles string
        :type: str
        '''
        return Chem.MolToSmiles(self.mol)

    def drawMoleculeAsByteStream(self, size=None):
        """Returns png image of molecule as bytestream

        :return: base64 png image bytestream
        :rtype: String
        """
        if size is not None:
            img = Chem.Draw.MolToImage(self.mol, size=size)
        else:
            img = Chem.Draw.MolToImage(self.mol)
        imgByteArray = io.BytesIO()
        img.save(imgByteArray, format='png')
        imgByteArray = imgByteArray.getvalue()
        imgByteArray = base64.b64encode(imgByteArray).decode("utf-8")
        return "data:;base64, " + imgByteArray

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
        all_warnings = {}

        for param, name in zip(params, filt_names):
            catalog = FilterCatalog(param)
            entries = catalog.GetMatches(self.mol)
            warnings = []
            for e in entries:
                warnings.append(e.GetDescription())
            if len(warnings) != 0:
                all_warnings[name] = warnings
            else:
                all_warnings[name] = 'passing'
        return all_warnings

    def descriptors(self):
        """Calculate molecule descriptor metrics as dict:
        | mol - RDKit molecule
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
        mw = Descriptors.ExactMolWt(self.mol)
        log_p = Crippen.MolLogP(self.mol)
        tpsa = rdMolDescriptors.CalcTPSA(self.mol)
        ha = Lipinski.HeavyAtomCount(self.mol)
        h_acceptors = Lipinski.NumHAcceptors(self.mol)
        h_donors = Lipinski.NumHDonors(self.mol)
        rings = Lipinski.RingCount(self.mol)
        desc_dict = {'mol': self.smiles,
                     'MW': mw,
                     'logP': log_p,
                     'TPSA': tpsa,
                     'HA': ha,
                     'h_acc': h_acceptors,
                     'h_don': h_donors,
                     'rings': rings
                     }
        return desc_dict

    def lipinski(self):
        """Calculate Lipinski from the descriptor dictionary. Returns the
        number of rules broken and whether the molecule passes.

        :return: violations
        :rtype: dict
        """
        desc_dict = self.descriptors()
        violations = {'MW': desc_dict['MW'] < 500.0,
                      'h_acc': desc_dict['h_acc'] <= 10,
                      'h_don': desc_dict['h_don'] <= 5,
                      'logP': desc_dict['logP'] < 5}
        return violations
