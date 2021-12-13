from flask import Flask, jsonify, request
from flask_cors import CORS

from Molecule import R_group, Molecule, FinalMolecule

app = Flask(__name__)
cors = CORS(app, resources={r"/*":
                            {"origins": "http://localhost:3000"}})

global saved_mols
# Can't serialize sets using json
saved_mols = []

global chosen_mol
chosen_mol = []

global money
money = [100000.0]

global time
time = [30.0]

global assay_prices
assay_prices = {
    "pIC50": 70.0,
    "c_mouse": 7000.0,
    "c_human": 9000.0,
    "LogD": 1000.0,
    "PAMPA": 700.0,
}

global assay_times
assay_times = {
    "pIC50": 1.0,
    "c_mouse": 3.0,
    "c_human": 3.5,
    "LogD": 1.5,
    "PAMPA": 1.0,
}

global assayed_molecules
assayed_molecules = {}


@app.route("/")
def hello_world():
    return "<p>Hello, World!</p>"


@app.route("/update_time_money")
def update_time_and_money():
    return jsonify((money[-1], time[-1]))


# Pass Yes/No for each assay as query, with the molecule as a key:
# /assays?r1=A01&r2=B10&pIC50=Yes&c_mouse=No&c_human=Yes&LogD=No&PAMPA=Yes
@app.route("/assays", methods=['POST'])
def run_assays():
    r_group_1_id = request.args.get('r1')
    r_group_2_id = request.args.get('r2')
    assay_list = []
    for label in ['pIC50', 'c_mouse', 'LogD', 'PAMPA']:
        if request.args.get(label) == "Yes":
            assay_list.append(label)
    if (r_group_1_id, r_group_2_id) in assayed_molecules:
        assayed_molecules[(r_group_1_id, r_group_2_id)].extend(assay_list)
    else:
        assayed_molecules[(r_group_1_id, r_group_2_id)] = assay_list
    if money[-1] - sum([assay_prices[p] for p in assay_list]) < 0:
        pass
    else:
        money.append(money[-1] - sum([assay_prices[p] for p in assay_list]))
    if time[-1] - max([assay_times[p] for p in assay_list]) < 0:
        pass
    else:
        time.append(time[-1] - max([assay_times[p] for p in assay_list]))
    return jsonify((money[-1], time[-1]))


# Pass molecule ID as query: /choose?r1=A01&r2=B10
@app.route("/choose", methods=['POST'])
def choose_molecule():
    if len(chosen_mol) > 0:
        return jsonify({'chosen_mol': chosen_mol})
    else:
        r_group_1_id = request.args.get('r1')
        r_group_2_id = request.args.get('r2')
        chosen_mol.append((r_group_1_id, r_group_2_id))
        return jsonify({'chosen_mol': chosen_mol})


@app.route("/chosenmolecule")
def return_chosen_molecules():
    return jsonify({'chosen_mol': chosen_mol})


# Pass molecule IDs as queries: /save?r1=A01&r2=B10
@app.route("/save", methods=['POST'])
def save_molecule():
    r_group_1_id = request.args.get('r1')
    r_group_2_id = request.args.get('r2')
    if (r_group_1_id, r_group_2_id) not in saved_mols:
        saved_mols.append((r_group_1_id, r_group_2_id))
    print(saved_mols)
    return jsonify({'saved_mols': saved_mols})


@app.route("/savedmolecules")
def return_saved_molecules():
    return jsonify({'saved_mols': saved_mols})


@app.route("/r-group-<string:r_group_id>")
def rgroup_img(r_group_id):
    """Returns image and stats R group specified by ID as a
    bytestream and dict to be rendered in a browser.

    :param r_group_id: ID number of R Group, eg. 'B26'
    :type r_group_id: String
    :return: Image and stats of R Group in a json dict.
    Access image bytestream with `img_html` key and stats with 'stats'
    :rtype: json dict
    """
    # smile_string = get_smile_string(r_group_id)
    # mol = Molecule(smile_string)
    mol = R_group(r_group_id)
    bytestream = mol.drawMoleculeAsByteStream()
    stats_dict = mol.descriptors()
    return jsonify({'img_html': f"data:;base64,{bytestream}",
                    'stats': stats_dict})


@app.route("/byte-img")
def byte_image():
    smile_string = 'O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]'
    mol = Molecule(smile_string)
    bytestream = mol.drawMoleculeAsByteStream()
    return f'<img src="data:;base64,{bytestream}"/>'


# Pass R group IDs as queries: /molecule?r1=A01&r2=B10
@app.route("/molecule")
def molecule_img():
    """Returns bytestream image  and drug properties
    of compund molecule consisting of
    scaffold and up to two R Groups, if these are specified in
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
    # molecule_smiles = scaffold_smiles

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


def add_r_group(base_molecule_as_smiles, r_group_as_smiles, r_group_nr=1):
    """Adds smile string of R group to scaffold smile string
    and returns combined compound as SMILE string.

    :param base_molecule_as_smiles: Scaffold molecule as SMILE string
    :type base_molecule_as_smiles: String
    :param r_group_as_smiles: R Group to be added to molecule as SMILE string
    :type r_group_as_smiles: String
    :param r_group_nr: R group number determines location where R-group is
    attached to scaffold, defaults to 1
    :type r_group_nr: int, optional
    :return: Compound molecule as SMILE string
    :rtype: String
    """
    new_mol = base_molecule_as_smiles + '.' + r_group_as_smiles
    connection_site = str(7+r_group_nr)
    new_mol = new_mol.replace(f'[*:{r_group_nr}]', connection_site)
    new_mol = new_mol.replace('('+connection_site+')', connection_site)
    return new_mol


def get_smile_string(r_group_id):
    """Retrieves smile string of rgroup by id

    :param r_group_id: R Group ID
    :type r_group_id: String
    :param group_number: R Group, defaults to 1
    :type group_number: int, optional
    :return: smile string of R group
    :rtype: String
    """
    group_number = 1 if r_group_id[0] == 'A' else 2
    try:
        mol = R_group(r_group_id, group_number)
        smile_string = mol.get_smile_string
    except ValueError('R group ID is invalid'):
        smile_string = None
    return smile_string
