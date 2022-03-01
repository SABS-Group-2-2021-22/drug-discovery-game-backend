from flask import Flask, jsonify, request
from flask_cors import CORS

from src.Molecule import R_group, Molecule, FinalMolecule

app = Flask(__name__)
cors = CORS(app, resources={r"/*":
                            {"origins": "http://localhost:3000"}})

# Temporary storage of data
global molecule_info
molecule_info = {}

global chosen_mol
chosen_mol = [None, None]

global money
money = [100000.0]

global time
time = [30.0]


@app.route("/")
def hello_world():
    return "<p>Hello, World!</p>"


# e.g. http://127.0.0.1:5000/get_all_mol_info
@app.route("/get_all_mol_info")
def get_all_mol_info():
    """Returns all assay, lipinski, filters and descriptor values for those run,
    for all molecules saved.
    Call just /get_all_mol_info

    :returns: Json dictionary of all saved information during a game
    :rtype: json dict
    """
    return jsonify(molecule_info)


# e.g. http://127.0.0.1:5000/update_time_money
@app.route("/update_time_money")
def update_time_and_money():
    """Updates the time and money values.
    Call just /update_time_money.

    :returns: A json tuple of the latest money and time value.
    :rtype: json tuple
    """
    return jsonify((money[-1], time[-1]))


def tuple2str(tuple_in):
    """Converts a tuple into a string.

    :param tuple_in: tuple to convert
    :type tuple_in: tuple
    :returns: concatenated string version of the tuple
    :rtype: str
    """
    string = ''
    for i in tuple_in:
        string += str(i)
    return string


# e.g. http://127.0.0.1:5000/lipinski?r1=A01&r2=B01
@app.route("/lipinski")
def run_lipinski():
    """Checks if molecule passes Lipinski's Rule of 5:
        MW < 500.0
        h_acc <= 10
        h_don <= 5
        logP < 5
    Saves the information to the global molecule_info dict and returns the
    information as its own dict.
    Pass R Group IDs as queries: /lipinski?r1=A01&r2=B01.

    :returns: A json dictionary of the molecule, indexed
    by the concatenated string of its R Group IDs, with the values being
    True or False for each descriptor relevant to Lipinski's Rule of 5.
    :rtype: json dict
    """
    lipinski_list = ['MW', 'logP', 'h_acc', 'h_don']
    r_group_1_id = request.args.get('r1')
    r_group_2_id = request.args.get('r2')
    molecule_key = tuple2str((r_group_1_id, r_group_2_id))
    drug_mol = FinalMolecule(r_group_1_id, r_group_2_id)
    drug_lipinski = drug_mol.lipinski(drug_mol.descriptors())
    lipinski_dict = {molecule_key: drug_lipinski}
    for label in lipinski_list:
        if "lipinski" in molecule_info[molecule_key].keys():
            pass
        else:
            molecule_info[molecule_key]["lipinski"] = {}
        molecule_info[molecule_key]["lipinski"][label] = drug_lipinski[label]
    return jsonify({"lipinski": lipinski_dict})


# e.g. http://127.0.0.1:5000/assays?r1=A01&r2=B01&pic50=Yes&clearance_mouse=No&
# clearance_human=Yes&logd=No&pampa=Yes
@app.route("/assays")
def run_assays():
    """Runs assays selected for a specific molecule, tracking the reduction of
    time and money as a result.
    The longest time of the assays being run is taken from the total amount of
    time.
    The sum of the cost of the assays is taken from the total amount of money.
    Pass R group IDs as queries: /assays?r1=A01&r2=B01... and Yes or No for
    each assay,
    ...&pic50=Yes&clearance_mouse=No&clearance_human=Yes&logd=No&pampa=Yes.

    :returns: A json dictionary of the molecule assayed, indexed
    by the concatenated string of its R Group IDs, with the values being
    which assays have been run and their values.
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
    print(assay_list)
    drug_mol = FinalMolecule(r_group_1_id, r_group_2_id)
    drug_properties = {label: drug_mol.drug_properties()[
        label] for label in assay_list}
    molecule_key = tuple2str((r_group_1_id, r_group_2_id))
    assay_dict = {}
    assay_dict[molecule_key] = {}

    for label in assay_list:
        if molecule_key in molecule_info.keys():
            if "assays" in molecule_info[molecule_key].keys():
                pass
            else:
                molecule_info[molecule_key]["assays"] = {}
            molecule_info[molecule_key]["assays"][label] = drug_properties[
                label
                ]
        assay_dict[molecule_key][label] = drug_properties[label]
    # Ensures that if run assay button is pressed, this code is not run and
    # so no crash occurs
    if molecule_key in molecule_info.keys() and len(assay_list) > 0:
        if money[-1] - sum([assay_prices[p] for p in assay_list]) < 0:
            pass
        else:
            money.append(money[-1] - sum(
                [assay_prices[p] for p in assay_list]))
        if time[-1] - max([assay_times[p] for p in assay_list]) < 0:
            pass
        else:
            time.append(time[-1] - max([assay_times[p] for p in assay_list]))
    return jsonify({"assays": assay_dict})


# e.g. http://127.0.0.1:5000/descriptors?r1=A01&r2=B01
@app.route("/descriptors")
def run_descriptors():
    """Runs descriptors ('MW', 'logP', 'TPSA', 'HA', 'h_acc', 'h_don', 'rings')
    for molecule selected.
    Saves the information to the global molecule_info dict and returns the
    information as its own dict.
    Pass R Group IDs as queries: /descriptors?r1=A01&r2=B01.

    :returns: A json dictionary of the molecule, indexed
    by the concatenated string of its R Group IDs, with the values for each
    descriptor, with each key being its respective descriptor label.
    :rtype: json dict
    """
    desc_list = ['MW', 'logP', 'TPSA', 'HA', 'h_acc', 'h_don', 'rings']
    r_group_1_id = request.args.get('r1')
    r_group_2_id = request.args.get('r2')
    drug_mol = FinalMolecule(r_group_1_id, r_group_2_id)
    drug_desc = {label: drug_mol.descriptors()[label] for label in desc_list}
    molecule_key = tuple2str((r_group_1_id, r_group_2_id))
    desc_dict = {}
    desc_dict[molecule_key] = {}
    for label in desc_list:
        if "descriptors" in molecule_info[molecule_key].keys():
            pass
        else:
            molecule_info[molecule_key]["descriptors"] = {}
        molecule_info[molecule_key]["descriptors"][label] = drug_desc[label]
        desc_dict[molecule_key][label] = drug_desc[label]
    return jsonify({"descriptors": desc_dict})


# e.g. http://127.0.0.1:5000/filters?r1=A01&r2=B01
@app.route("/filters")
def run_filters():
    """Runs filters ('PAINS', 'ZINC', 'BRENK', 'NIH')for molecule selected.
    Saves the information to the global molecule_info dict and returns the
    information as its own dict.
    Pass R Group IDs as queries: /filters?r1=A01&r2=B01

    :returns: A json dictionary of the molecule, indexed
    by the concatenated string of its R Group IDs, with the values for each
    descriptor, with each key being its respective descriptor label.
    :rtype: json dict
    """
    filter_names = ['PAINS', 'ZINC', 'BRENK', 'NIH']
    r_group_1_id = request.args.get('r1')
    r_group_2_id = request.args.get('r2')
    drug_mol = FinalMolecule(r_group_1_id, r_group_2_id)
    drug_filters = drug_mol.filter_properties()
    molecule_key = tuple2str((r_group_1_id, r_group_2_id))
    filt_dict = {}
    filt_dict[molecule_key] = {}
    for label in filter_names:
        if "filters" in molecule_info[molecule_key].keys():
            pass
        else:
            molecule_info[molecule_key]["filters"] = {}
        molecule_info[molecule_key]["filters"][label] = drug_filters[label]
        filt_dict[molecule_key][label] = drug_filters[label]
    return jsonify({"filter_dict": filt_dict})


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
    return jsonify({'chosen_mol': chosen_mol})


@app.route("/chosenmolecule")
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


# http://127.0.0.1:5000/save?r1=A01&r2=B01
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
    molecule_key = tuple2str((r_group_1_id, r_group_2_id))
    if molecule_key not in molecule_info:
        molecule_info[molecule_key] = {}
        molecule_info[molecule_key]["keys"] = (r_group_1_id, r_group_2_id)
    list_saved_mols = [molecule_info[m]["keys"] for m in molecule_info.keys()]
    return jsonify({'saved_mols': list_saved_mols})


# http://127.0.0.1:5000/savedmolecules
@app.route("/savedmolecules")
def return_saved_molecules():
    """Returns the list of saved molecules as json dict of a list, containing
    tuples of the R Group IDs. Currently a global variable in place of database
    storage.
    Call just /savedmolecules.

    :return: List of tuples, containing the R Group IDs as a json dict.
    Access list with 'saved_mols' key.
    :rtype: json dict
    """
    list_saved_mols = [molecule_info[m]["keys"] for m in molecule_info.keys()]
    return jsonify({'saved_mols': list_saved_mols})


@app.route("/r-group-<string:r_group_id>")
def rgroup_img(r_group_id):
    """Returns image and stats R group specified by ID as a bytestream and dict
    to be rendered in a browser.
    Pass R Group ID as query: /r-group-A01

    :param r_group_id: ID number of R Group, eg. 'B26'
    :type r_group_id: String
    :return: Image and stats of R Group in a json dict.
    Access image bytestream with `img_html` key and stats with 'stats'
    :rtype: json dict
    """
    mol = R_group(r_group_id)
    bytestream = mol.drawMoleculeAsByteStream()
    stats_dict = mol.descriptors()
    return jsonify({'img_html': f"data:;base64,{bytestream}",
                    'stats': stats_dict})


@app.route("/molecule")
def molecule_img():
    """Returns bytestream image  and drug properties of compund molecule
    consisting of scaffold and up to two R Groups, if these are specified in
    the query parameters.
    Pass R Group IDs and image size as queries:
    /molecule?r1=A01&r2=B10&size=800,800
    If no R Groups are specified this returns the scaffold.

    :return: JSON containing image of compound molecule as bytestream.
    :rtype: json dict
    """
    r_group_1_id = request.args.get('r1')
    r_group_2_id = request.args.get('r2')
    size = eval(request.args.get('size'))

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
        orient_with_scaffold=True, size=size
    )
    if r_group_1_id is not None and r_group_2_id is not None:
        drug_mol = FinalMolecule(r_group_1_id, r_group_2_id)
        drug_property_dict = drug_mol.drug_properties()
    return jsonify({'img_html': f"data:;base64,{bytestream}",
                    'drug_props': drug_property_dict,
                    'descriptors': {},
                    'filters': {},
                    'assays_run': {
                       # 'pic50': 'N',
                       # 'clearance_mouse': 'N',
                       # 'clearance_human': 'N',
                       # 'logd': 'N',
                       # 'pampa': 'N',
                       # 'filters': 'N',
                       # 'descriptors': 'N'
                    }})


@app.route("/getplotdata")
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


@app.route("/getspiderdata")
def return_spider_data():
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


def numerise_params(prop_dict):
    """ Returns drug properties with all qualitative values transformed into
    numeric values

    returns: numerically transformed property dictionaries
    rtype: dict
    """
    clearance_dict = {
                'low (< 5.6)': 1,
                'medium (5.6-30.5)': 4,
                'low (< 3.7)': 1,
                'good': 1,
                'high (> 30.5)': 7,
                'fair': 4,
                'poor': 7,
                'low (< 12)': 1,
                'medium (12-44)': 4,
                'medium (5.6-30.5)': 4
    }
    pampa_dict = {
                'neg': 0,
                'poor': 1,
                'low': 2.5,
                'fair': 5.5,
                'med2high': 5.5,
                'good': 6.5,
                'best': 8
    }
    drug_properties = prop_dict

    for k, v in clearance_dict.items():
        if k == drug_properties['clearance_mouse']:
            drug_properties['clearance_mouse'] = v
        if k == drug_properties['clearance_human']:
            drug_properties['clearance_human'] = v
    for k, v in pampa_dict.items():
        if k == drug_properties['pampa']:
            drug_properties['pampa'] = v
        if k == drug_properties['logd']:
            drug_properties['logd'] = v
    return (drug_properties)


@app.route("/comparisontxt")
def comparison_txt():
    """ Returns comparison text depending on pic50, logd, and
    clearance_human of chosen molecule

    returns: json dict with text in value depending on metrix
    rtype: json dict
    """
    assay_list = ['pic50', 'clearance_mouse', 'clearance_human',
                  'logd', 'pampa']

    if chosen_mol[0] is not None and chosen_mol[1] is not None:
        r_group_1_id, r_group_2_id = chosen_mol[0], chosen_mol[1]
    else:
        r_group_1_id, r_group_2_id = 'A01', 'B01'
    drug_mol = FinalMolecule(r_group_1_id, r_group_2_id)
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
    print(comp_dict)
    return jsonify({'comparison': comp_dict})


@app.route("/reset")
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
