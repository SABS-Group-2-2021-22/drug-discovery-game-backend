from flask import Flask, jsonify, request
from flask_cors import CORS

from Molecule import R_group, Molecule

app = Flask(__name__)
cors = CORS(app, resources={r"/*":
                            {"origins": "http://localhost:3000"}})


@app.route("/")
def hello_world():
    return "<p>Hello, World!</p>"


@app.route("/r-group-<string:r_group_id>")
def rgroup_img(r_group_id):
    """Returns image of R group specified by ID as a
    bytestream to be rendered in a browser.

    :param r_group_id: ID number of R Group, eg. 'B26'
    :type r_group_id: String
    :return: Image of R Group as bytstream in a json dict.
    Access image bytestream with `img_html` key
    :rtype: json dict
    """
    smile_string = get_smile_string(r_group_id)
    mol = Molecule(smile_string)
    bytestream = mol.drawMoleculeAsByteStream()
    return jsonify({'img_html': f"data:;base64,{bytestream}"})


@app.route("/byte-img")
def byte_image():
    smile_string = 'O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]'
    mol = Molecule(smile_string)
    bytestream = mol.drawMoleculeAsByteStream()
    return f'<img src="data:;base64,{bytestream}"/>'


# Pass R group IDs as queries: /molecule?r1=A01&r2=B10
@app.route("/molecule")
def molecule_img():
    """Returns bytestream image of compund molecule consisting of
    scaffold and up to two R Groups, if these are specified in
    the query parameters.
    Pass R group IDs as queries: /molecule?r1=A01&r2=B10
    If no R groups are specified this returns the scaffold.

    :return: JSON containing image of compound molecule as bytestream.
    :rtype: json dict
    """
    r_group_1_id = request.args.get('r1')
    r_group_2_id = request.args.get('r2')

    scaffold_smiles = 'O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]'
    molecule_smiles = scaffold_smiles

    # Add R group one
    if r_group_1_id is not None:        # Check R1 was included in query
        r_group_1 = get_smile_string(r_group_1_id)
        if r_group_1 is not None:       # Check R1 id was valid
            molecule_smiles = add_r_group(molecule_smiles,
                                          r_group_1,
                                          r_group_nr=1)

    # Add R group one
    if r_group_2_id is not None:        # Check R2 was included in query
        r_group_2 = get_smile_string(r_group_2_id)
        if r_group_2 is not None:       # Check R2 id was valid
            molecule_smiles = add_r_group(molecule_smiles,
                                          r_group_2,
                                          r_group_nr=2)

    molecule = Molecule(molecule_smiles)
    bytestream = molecule.drawMoleculeAsByteStream()
    return jsonify({'img_html': f"data:;base64,{bytestream}"})


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
