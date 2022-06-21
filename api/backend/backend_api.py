
from flask import jsonify, request
from .utils import tuple2str

# Temporary storage of data
global molecule_info
molecule_info = {}

global chosen_mol
chosen_mol = [None, None]

global money
money = [100000.0]

global time
time = [30.0]


def get_all_mol_info(molecule_info):
    """Returns all assay, lipinski, filters and descriptor values for those run,
    for all molecules saved.
    Call just /get_all_mol_info

    :returns: Json dictionary of all saved information during a game
    :rtype: json dict
    """
    return jsonify(molecule_info)


def choose_molecule():
    """Saves the final choice of molecule for the end of the game as a tuple of
    the final molecule's R Group IDs.
    Pass R group IDs as queries: /choose?r1=A01&r2=B10

    :returns: Tuple of the R Group IDs as a json dict. Access tuple with
    'chosen_mol' key.
    :rtype: json dict
    """
    chosen_mol = [request.args.get('r1'), request.args.get('r2')]
    return jsonify({'chosen_mol': chosen_mol}), chosen_mol


def save_molecule(mol_dict):
    """Saves a molecule from the front-end, storing in the back-end,
    representing it as a pair of R Group IDs. Uses POST API call.
    Pass R group IDs as queries: /save?r1=A01&r2=B10

    :return: List of tuples, containing the R Group IDs as a json dict. Access
    list with 'saved_mols' key
    :rtype: json dict
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


def reset():
    """ Resets values and so restarts the game

    returns: The new empty molecule_info dictionary as a json dict
    rtype: json dict
    """
    molecule_info.clear()
    chosen_mol.clear()
    chosen_mol.append(None)
    chosen_mol.append(None)
    money.clear()
    money.append(100000.0)
    time.clear()
    time.append(30.0)
    return jsonify({"new_info": molecule_info})
