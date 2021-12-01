from flask import Flask, jsonify
from flask_cors import CORS

from src.Molecule import R_group, Molecule

app = Flask(__name__)
cors = CORS(app, resources={r"/r-group-*":
                            {"origins": "http://localhost:3000"}})


@app.route("/")
def hello_world():
    return "<p>Hello, World!</p>"


@app.route("/r-group-<string:r_group_id>")
def rgroup_img(r_group_id):
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


def get_smile_string(r_group_id, group_number=1):
    """Retrieves smile string of rgroup by id

    :param r_group_id: R Group ID
    :type r_group_id: String
    :param group_number: R Group, defaults to 1
    :type group_number: int, optional
    :return: smile string of R group
    :rtype: String
    """
    mol = R_group(r_group_id, group_number)
    smile_string = mol.get_smile_string
    return smile_string
