
from flask import jsonify, request

from .utils import tuple2str, numerise_params

# Temporary storage of data
global molecule_info
molecule_info = {}

global chosen_mol
chosen_mol = [None, None]

global money
money = [100000.0]

global time
time = [30.0]


def hello_world():
    return "<p>Hello, World!</p>"


def get_all_mol_info(molecule_info):
    """Returns all assay, lipinski, filters and descriptor values for those run,
    for all molecules saved.
    Call just /get_all_mol_info

    :returns: Json dictionary of all saved information during a game
    :rtype: json dict
    """
    return jsonify(molecule_info)


def update_time_and_money(user):
    """Updates the time and money values.
    Call just /update_time_money.

    :returns: A json tuple of the latest money and time value.
    :rtype: json tuple
    """
    user.update_time_and_money()
    return jsonify((money[-1], time[-1]))


def choose_molecule():
    """Saves the final choice of molecule for the end of the game as a tuple of
    the final molecule's R Group IDs.
    Pass R group IDs as queries: /choose?r1=A01&r2=B10

    :returns: Tuple of the R Group IDs as a json dict. Access tuple with
    'chosen_mol' key.
    :rtype: json dict
    """
    # Prevents overwrite if a molecule has already been chosen
    # if len(chosen_mol) > 0:
    #     return jsonify({'chosen_mol': chosen_mol})
    # else:
    #     r_group_1_id = request.args.get('r1')
    #     r_group_2_id = request.args.get('r2')
    #     chosen_mol.append((r_group_1_id, r_group_2_id))
    # r_group_1_id = request.args.get('r1')
    # r_group_2_id = request.args.get('r2')
    chosen_mol[0] = request.args.get('r1')
    chosen_mol[1] = request.args.get('r2')
    return jsonify({'chosen_mol': chosen_mol}), chosen_mol


def return_chosen_molecules():
    """Returns the final choice of molecule for the end of the game as a json
    dict containing the tuple of the final molecule's R Group IDs.
    Call just /chosenmolecule.

    :returns: Tuple of the R Group IDs as a json dict. Access tuple with
    'chosen_mol' key.
    :rtype: json dict
    """
    if chosen_mol[0] is not None and chosen_mol[1] is not None:
        return jsonify({'chosen_mol': chosen_mol})
    else:
        return jsonify({})


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


def return_saved_molecules(molecule_info):
    """Returns the list of saved molecules as json dict of a list, containing
    tuples of the R Group IDs. Currently a global variable in place of database
    storage.
    Call just /savedmolecules.

    :return: List of tuples, containing the R Group IDs as a json dict.
    Access list with 'saved_mols' key.
    :rtype: json dict
    """
    list_saved_mols = [molecule_info[m]["keys"]
                       for m in molecule_info.keys()]
    return jsonify({'saved_mols': list_saved_mols})


def return_assayed_data():
    """Returns molecule_info in a restructured format, to facilatate easier
    plotting.
    Call just /getplotdata.

    :returns: A json dictionary, with a key 'assay_dict' and a list as a value.
    :rtype: json dict
    """
    data = {}
    for k, v in molecule_info.items():
        try:
            assay_dict = v['assays']
        except Exception:
            assay_dict = {}
        try:
            descriptor_dict = v['descriptors']
        except Exception:
            descriptor_dict = {}
        metric_dict = {**assay_dict, **descriptor_dict}
        data[k] = metric_dict
    for k, v in data.items():
        v['--'] = 0
    return jsonify({'assay_dict': [data]})


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
