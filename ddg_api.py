import flask
from flask import Flask
import base64

import Molecule

app = Flask(__name__)

@app.route("/")
def hello_world():
    return "<p>Hello, World!</p>"


@app.route("/r-group-<string:r_group_id>")
def rgroup_img(r_group_id):
    img_filename = f'./images/{r_group_id}.png'
    return flask.send_file(img_filename, mimetype='image/png')


@app.route("/byte-img")
def byte_image():
    smile_string = 'O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]'
    mol = Molecule.Molecule(smile_string)
    bytestream = mol.drawMoleculeAsByteStream()
    return f'<img src="data:;base64,{bytestream}"/>'

@property
def get_smile_string(r_group_id, group_number=1):
    """Retrieves smile string of rgroup by id

    :param r_group_id: R Group ID
    :type r_group_id: String
    :param group_number: R Group, defaults to 1
    :type group_number: int, optional
    :return: [description]
    :rtype: [type]
    """
    smile_string = 'XXX'
    
    return smile_string





def something():
    """Retrieves smile string of rgroup by id

    :param r_group_id: R Group ID
    :type r_group_id: String
    :param r_group_id: R Group ID
    :type r_group_id: String

    """
