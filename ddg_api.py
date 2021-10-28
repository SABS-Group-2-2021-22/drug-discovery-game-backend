import flask
from flask import Flask
from flask_cors import CORS
import base64
import io

from flask import send_file

import Molecule
from Molecule import R_group

app = Flask(__name__)
cors = CORS(app, resources={r"/r-group-*": {"origins": "http://localhost:3000"}})
# CORS

@app.route("/")
def hello_world():
    return "<p>Hello, World!</p>"


@app.route("/r-group-<string:r_group_id>")
def rgroup_img(r_group_id):
    smile_string = get_smile_string(r_group_id)
    mol = Molecule.Molecule(smile_string)
    bytestream = mol.drawMoleculeAsByteStream()
    return send_file(
        io.BytesIO(bytestream),
        mimetype='image/png'
    )    
    # return f'<img src="data:;base64,{bytestream}"/>'
    # return f"data:;base64,{bytestream}"
    # img_filename = f'./images/{r_group_id}.png'
    # return flask.send_file(img_filename, mimetype='image/png')


@app.route("/byte-img")
def byte_image():
    smile_string = 'O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]'
    mol = Molecule.Molecule(smile_string)
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
    # smile_string = Molecule.R_group(r_group_id,
                                    # group_number).extract_smilefromcsv()
    return smile_string
