
from flask import jsonify, request
from .utils import tuple2str


def choose_molecule():
    """Saves the final choice of molecule for the end of the game as a list of
    the final molecule's R Group IDs.

    :returns: Json dict contained list of the R Group IDs and the list on its
    own.
    :rtype: json dict, list
    """
    chosen_mol = [request.args.get('r1'), request.args.get('r2')]
    return jsonify({'chosen_mol': chosen_mol}), chosen_mol


def save_molecule(mol_dict):
    """Saves a molecule representing it as a pair of R Group IDs. Appends it to
    mol_dict and returns updated mol_dict.

    :param mol_dict: All the information about the current state of the game
    :type mol_dict: dict
    :return: List of tuples, containing the R Group IDs as a json dict and
    updated mol_dict.
    :rtype: json dict, dict
    """
    r_group_1_id = request.args.get('r1')
    r_group_2_id = request.args.get('r2')
    molecule_key = tuple2str((r_group_1_id, r_group_2_id))
    if molecule_key not in mol_dict:
        mol_dict[molecule_key] = {}
        mol_dict[molecule_key]["keys"] = (r_group_1_id, r_group_2_id)
    list_saved_mols = [mol_dict[m]["keys"]
                       for m in mol_dict.keys()]
    return jsonify({'saved_mols': list_saved_mols}), mol_dict


def reset(mol_dict):
    """ Resets values and so restarts the game

    :param mol_dict: All the information about the current state of the game
    :type mol_dict: dict
    :return: The new empty molecule_info dictionary as a json dict
    :rtype: json dict
    """
    return jsonify({"new_info": mol_dict})
