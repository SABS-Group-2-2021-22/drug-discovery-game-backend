from flask import jsonify, request
import rdkit
from src.sketchedMolecule import sketchedMolecule
from src.Molecule import FinalMolecule
from .utils import numerise_params


def sketcher_save_molecule(mol_block, mol_dict):
    '''Receives Mol block of sketched molecules and returns descriptors,
    filters, lipinski rules,  assay information, smiles, tanimoto similarity
    to Roche's drug and Lipinski rules and image of the molecule

    :returns: A json dictionary of the bytestream, descriptors, filters and
    Lipinksi rules for the molecule
    :rtypes: json dict
    '''
    try:
        new_mol = sketchedMolecule(mol_block)
        molecule_key = new_mol.smiles
        new_dict = {'smiles': new_mol.smiles,
                    'img_html': new_mol.drawMoleculeAsByteStream(),
                    'descriptors': new_mol.descriptors(),
                    'filters': new_mol.filter_properties(),
                    'lipinski': new_mol.lipinski(),
                    'assays_run': {},
                    'drug_props': new_mol.drug_properties(),
                    }
        if molecule_key not in mol_dict:
            mol_dict[molecule_key] = new_dict

        return jsonify(new_dict), mol_dict
    except:  # noqa: E722
        return jsonify('failure'), mol_dict


def sketcher_choose_molecule():
    """Saves the final choice of molecule for the end of the game as a tuple of
    the final molecule's R Group IDs.
    Pass R group IDs as queries: /choose?r1=A01&r2=B10

    :returns: Tuple of the R Group IDs as a json dict. Access tuple with
    'chosen_mol' key.
    :rtype: json dict
    """
    chosen_mol = [request.args.get('id'), request.args.get('smiles')]
    return jsonify({'chosen_mol': chosen_mol}), chosen_mol


def sketcher_comparison_txt(chosen_mol):
    """ Returns comparison text depending on pic50, logd, and
    clearance_human of chosen molecule

    returns: json dict with text in value depending on metric
    rtype: json dict
    """
    assay_list = ['pic50', 'clearance_mouse', 'clearance_human',
                  'logd', 'pampa']
    mol = rdkit.Chem.MolFromSmiles(chosen_mol[1])
    mol_block = rdkit.Chem.rdmolfiles.MolToMolBlock(mol)
    drug_mol = sketchedMolecule(mol_block)
    drug_properties = {
        label: drug_mol.drug_properties()[label] for label in assay_list
    }
    # some properties for logd are strings ('best' or 'good') (no idea why)
    drug_properties = numerise_params(drug_properties)

    comp_dict = {}
    with open('./src/comparison.txt', 'r') as f:
        lines = f.readlines()
        if float(drug_properties['pic50']) < 6.5:
            comp_dict['pic50'] = str(lines[0])
        else:
            comp_dict['pic50'] = str(lines[1])
        if float(drug_properties['logd']) < 0.95:
            comp_dict['logd'] = str(lines[2])
        elif 0.95 < float(drug_properties['logd']) < 1.15:
            comp_dict['logd'] = str(lines[3])
        else:
            comp_dict['logd'] = str(lines[4])
        if drug_properties['clearance_human'] != str(1):
            comp_dict['clearance_human'] = str(lines[5])
        else:
            comp_dict['clearance_human'] = str(lines[6])
    return jsonify({'comparison': comp_dict})


def sketcher_return_spider_data(chosen_mol):
    """Takes final chosen molecule (chosen_mol) and calls FinalMolecule()
     from Molecule.py to return quantitative assay parameters of the chosen
     molecule and the reference drug
     Calls /getspiderdata + FinalMolecule() + numerise_params()

     :returns: A json dictionary containing a list of 2 dictionaries, one
     containing chosen mol.parameters and the other containing reference
     drug parameters
     rtype: json dict
     """

    assay_list = ['pic50', 'clearance_mouse', 'clearance_human',
                  'logd', 'pampa']
    mol = rdkit.Chem.MolFromSmiles(chosen_mol[1])
    mol_block = rdkit.Chem.rdmolfiles.MolToMolBlock(mol)
    drug_mol = sketchedMolecule(mol_block)

    drug_properties = {
        label: drug_mol.drug_properties()[label] for label in assay_list
    }
    ref_mol = FinalMolecule('A05', 'B07')
    ref_properties = {
        label: ref_mol.drug_properties()[label] for label in assay_list
    }
    drug_properties = numerise_params(drug_properties)
    ref_properties = numerise_params(ref_properties)
    property_arr = [drug_properties, ref_properties]
    return jsonify({'param_dict': property_arr})
