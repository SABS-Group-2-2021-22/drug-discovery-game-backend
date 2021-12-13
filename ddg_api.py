from flask import Flask, jsonify, request
from flask_cors import CORS

from src.Molecule import R_group, Molecule, FinalMolecule

app = Flask(__name__)
cors = CORS(app, resources={r"/*":
                            {"origins": "http://localhost:3000"}})

# Temporary storage of data
global saved_mols
saved_mols = []

global assayed_molecules
assayed_molecules = {}

global chosen_mol
chosen_mol = []

global money
money = [100000.0]

global time
time = [30.0]


@app.route("/")
def hello_world():
    return "<p>Hello, World!</p>"


@app.route("/update_time_money")
def update_time_and_money():
    """Updates the time and money values

    :returns: A json tuple of the latest money and time value.
    :rtype: json tuple
    """
    return jsonify((money[-1], time[-1]))


@app.route("/assays")
def run_assays():
    """Runs assays selected for a specific molecule, tracking the reduction of
    time and money as a result.
    The longest time of the assays being run is taken from the total amount of
    time.
    The sum of the cost of the assays is taken from the total amount of money.
    Pass R group IDs as queries: /assays?r1=A01&r2=B10... and Yes or No for
    each assay,
    ...&pic50=Yes&clearance_mouse=No&clearance_human=Yes&logd=No&pampa=Yes.

    :returns: A json dictionary of which molecules have been assayed, indexed
    by the concatenated string of their R Group IDs, with the values being
    which assays have been run and their values for all assays
    :rtype: json dict
    """

    # The prices of each assay
    assay_prices = {
        "pic50": 70.0,
        "clearance_mouse": 7000.0,
        "clearance_human": 9000.0,
        "logd": 1000.0,
        "pampa": 700.0,
    }
    # How long each assay takes in weeks
    assay_times = {
        "pic50": 1.0,
        "clearance_mouse": 3.0,
        "clearance_human": 3.5,
        "logd": 1.5,
        "pampa": 1.0,
    }
    r_group_1_id = request.args.get('r1')
    r_group_2_id = request.args.get('r2')
    assay_list = []
    # Make a list of the assays being run
    for label in ['pic50',
                  'clearance_mouse',
                  'clearance_human',
                  'logd',
                  'pampa']:
        if request.args.get(label) == "Yes":
            assay_list.append(label)

    drug_mol = FinalMolecule(r_group_1_id, r_group_2_id)
    drug_properties = {label: drug_mol.drug_properties()[
        label] for label in assay_list}
    molecule_key = tuple2str((r_group_1_id, r_group_2_id))
    if molecule_key not in assayed_molecules.keys():
        assayed_molecules[molecule_key] = drug_properties
    else:
        for label in assay_list:
            assayed_molecules[molecule_key][label] = drug_properties[label]
    if money[-1] - sum([assay_prices[p] for p in assay_list]) < 0:
        pass
    else:
        money.append(money[-1] - sum([assay_prices[p] for p in assay_list]))
    if time[-1] - max([assay_times[p] for p in assay_list]) < 0:
        pass
    else:
        time.append(time[-1] - max([assay_times[p] for p in assay_list]))
    return jsonify({'assay_dict': assayed_molecules})


def tuple2str(tuple_in):
    """Converts a tuple into a string

    :param tuple_in: tuple to convery
    :type tuple_in: tuple
    :returns: concatenated string version of the tuple
    :rtype: str
    """
    string = ''
    for i in tuple_in:
        string += str(i)
    return string


@app.route("/choose", methods=['POST'])
def choose_molecule():
    """Saves the final choice of molecule for the end of the game as a tuple of
    the final molecule's R Group IDs.
    Pass R group IDs as queries: /choose?r1=A01&r2=B10

    :returns: Tuple of the R Group IDs as a json dict. Access tuple with
    'chosen_mol' key.
    :rtype: json dict
    """
    # Prevents overwrite if a molecule has already been chosen
    if len(chosen_mol) > 0:
        return jsonify({'chosen_mol': chosen_mol})
    else:
        r_group_1_id = request.args.get('r1')
        r_group_2_id = request.args.get('r2')
        chosen_mol.append((r_group_1_id, r_group_2_id))
        return jsonify({'chosen_mol': chosen_mol})


@app.route("/chosenmolecule")
def return_chosen_molecules():
    """Returns the final choice of molecule for the end of the game as a json
    dict containing the tuple of the final molecule's R Group IDs.

    :returns: Tuple of the R Group IDs as a json dict. Access tuple with
    'chosen_mol' key.
    :rtype: json dict
    """
    return jsonify({'chosen_mol': chosen_mol})


@app.route("/save", methods=['POST'])
def save_molecule():
    """Saves a molecule from the front-end, storing in the back-end,
    representing it as a pair of R Group IDs. Uses POST API call.
    Pass R group IDs as queries: /save?r1=A01&r2=B10

    :return: List of tuples, containing the R Group IDs as a json dict. Access
    list with 'saved_mols' key
    :rtype: json dict
    """
    r_group_1_id = request.args.get('r1')
    r_group_2_id = request.args.get('r2')
    if (r_group_1_id, r_group_2_id) not in saved_mols:
        saved_mols.append((r_group_1_id, r_group_2_id))
    print(saved_mols)
    return jsonify({'saved_mols': saved_mols})


@app.route("/savedmolecules")
def return_saved_molecules():
    """Returns the list of saved molecules as json dict of a list, containing
    tuples of the R Group IDs. Currently a global variable in place of database
    storage

    :return: List of tuples, containing the R Group IDs as a json dict.
    Access list with 'saved_mols' key.
    :rtype: json dict
    """
    return jsonify({'saved_mols': saved_mols})


# Pass R group IDs as queries: /molecule?r1=A01&r2=B10
@app.route("/molecule")
def molecule_img():
    """Returns bytestream image  and drug properties of compund molecule
    consisting of scaffold and up to two R Groups, if these are specified in
    the query parameters.
    Pass R group IDs as queries: /molecule?r1=A01&r2=B10
    If no R groups are specified this returns the scaffold.

    :return: JSON containing image of compound molecule as bytestream.
    :rtype: json dict
    """
    r_group_1_id = request.args.get('r1')
    r_group_2_id = request.args.get('r2')

    r_group_1 = None
    r_group_2 = None

    scaffold_smiles = 'O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]'

    base_molecule = Molecule(scaffold_smiles)

    # Add R group one
    if r_group_1_id is not None:        # Check R1 was included in query
        r_group_1 = R_group(r_group_1_id)
        if r_group_1 is not None:       # Check R1 id was valid
            base_molecule = r_group_1.add_r_group(base_molecule)

    # Add R group one
    if r_group_2_id is not None:        # Check R2 was included in query
        r_group_2 = R_group(r_group_2_id)
        if r_group_2 is not None:       # Check R2 id was valid
            base_molecule = r_group_2.add_r_group(base_molecule)

    bytestream = base_molecule.drawMoleculeAsByteStream(
        orient_with_scaffold=True, size=(800, 800)
    )
    if r_group_1_id is not None and r_group_2_id is not None:
        drug_mol = FinalMolecule(r_group_1_id, r_group_2_id)
        drug_property_dict = drug_mol.drug_properties()
    return jsonify({'img_html': f"data:;base64,{bytestream}",
                    'drug_props': drug_property_dict})
