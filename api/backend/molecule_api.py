

from flask import jsonify, request

from src.Molecule import R_group, Molecule, FinalMolecule

from .utils import tuple2str, numerise_params

# global molecule_info
# molecule_info = {}

global chosen_mol
chosen_mol = [None, None]

global money
money = [100000.0]

global time
time = [30.0]


# TODO: refactor query passing here, direct via function args
def run_lipinski(molecule_info):
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
            print()
        else:
            molecule_info[molecule_key]["lipinski"] = {}
        molecule_info[molecule_key]["lipinski"][label] = drug_lipinski[label]
    return jsonify({"lipinski": lipinski_dict}), molecule_info


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
            time.append(time[-1] - max([assay_times[p]
                        for p in assay_list]))
    return jsonify({"assays": assay_dict})


def run_descriptors(molecule_info):
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
    drug_desc = {label: drug_mol.descriptors()[
        label] for label in desc_list}
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
    return jsonify({"descriptors": desc_dict}), molecule_info


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
                    'assays_run': {}
                    })


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
