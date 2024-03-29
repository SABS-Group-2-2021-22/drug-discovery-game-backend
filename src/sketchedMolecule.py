from rdkit import Chem
from rdkit.Chem import Lipinski,\
    Descriptors, rdMolDescriptors, Crippen, rdDepictor
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
import pandas as pd
import io
import base64
import pickle


class sketchedMolecule:
    def __init__(self, mol_block):
        """Initialises sketchedMolecule class using Mol block. If Mol block is
        not valid, will flag error.

        :param mol_block: Mol Block of molecule
        :type mol_block: str
        """
        self.mol_block = mol_block
        try:
            self.mol = Chem.rdmolfiles.MolFromMolBlock(self.mol_block)
            old_mol = Chem.rdmolfiles.MolFromMolBlock(self.mol_block)
            rdDepictor.Compute2DCoords(self.mol, clearConfs=True)
            Chem.rdMolAlign.AlignMol(self.mol, old_mol)
        except:  # noqa: E722
            print('The Mol block is either incorrectly structured or the compound cannot be \
                processed by RDKit')

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
        :rtype: str
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
        """Calculates whether molecule passes or fails PAINS filters (PAINS,
        ZINC, BRENK, NIH).

        :returns: Dictionary of if it passes each type of PAINS filters for
        each case.
        :rtype: dict
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

        :param desc_dict: Molecule descriptor metrics calculated by
        descriptors().
        :type desc_dict: dict
        :return: Dictionary of it passes the Lipinski Rules for each case.
        :rtype: dict
        """
        desc_dict = self.descriptors()
        violations = {'MW': desc_dict['MW'] < 500.0,
                      'h_acc': desc_dict['h_acc'] <= 10,
                      'h_don': desc_dict['h_don'] <= 5,
                      'logP': desc_dict['logP'] < 5}
        return violations

    def drug_properties(self):
        """Predicts properties of the final drug from the models trained on the
        MMP12 dataset (assay_ml_models). See assays_models.ipynb.

        :return: dictionary with each property as a key and the predicted value
        from the ML models as the value
        :rtype: dict
        """
        drug_property_dict = {}
        drug_properties = [
            'pic50',
            'clearance_mouse',
            # 'clearance_human', # No ML model needed for this
            'logd',
            'pampa'
        ]
        descriptors = {d[0]: d[1] for d in Descriptors.descList}
        rdkit_features = pd.read_csv('r_group_rdkit_features.csv', index_col=0)
        rdkit_desc_dict = {}
        rdkit_desc_dict[0] = \
            {r: descriptors[r](self.mol) for r in rdkit_features.columns}
        desc_df = pd.DataFrame.from_dict(rdkit_desc_dict)
        desc_df = desc_df.T
        for d in drug_properties:
            model = pickle.load(
                open(f'assay_ml_models/automl_{d}_model.pkl', 'rb'))
            value = model.predict(desc_df)[0]
            if isinstance(value, float):
                value = "{:.2f}".format(value)
            drug_property_dict[d] = value
        # All values in the csv file are this so just set molecule as the same
        drug_property_dict['clearance_human'] = 'low (< 12)'
        return drug_property_dict
