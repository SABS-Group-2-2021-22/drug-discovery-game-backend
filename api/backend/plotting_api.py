from src.Molecule import FinalMolecule
from .utils import numerise_params

from flask import jsonify


def return_spider_data(chosen_mol):
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

    if chosen_mol[0] is not None and chosen_mol[1] is not None:
        r_group_1_id = chosen_mol[0]
        r_group_2_id = chosen_mol[1]
    else:
        r_group_1_id = 'A01'
        r_group_2_id = 'B01'

    drug_mol = FinalMolecule(r_group_1_id, r_group_2_id)

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