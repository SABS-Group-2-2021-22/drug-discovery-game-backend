from flask import Flask, jsonify, request
from flask_cors import CORS

from api.backend.backend_api import API 

app = Flask(__name__)
cors = CORS(app, resources={r"/*":
                            {"origins": "http://localhost:3000"}})


@app.route("/")
def hello_world():
    return API.hello_world()


@app.route("/get_all_mol_info")
def get_all_mol_info():
    return API.get_all_mol_info()


@app.route("/update_time_money")
def update_time_and_money():
    return API.update_time_and_money()


# TODO: refactor query passing here, direct via function args
@app.route("/lipinski")
def run_lipinski():
    return API.run_lipinski()

@app.route("/assays")
def run_assays():
    return API.run_assays()

@app.route("/descriptors")
def run_descriptors():
    return API.run_descriptors()


@app.route("/filters")
def run_filters():
    return API.run_filters()


@app.route("/choose", methods=['GET', 'POST'])
def choose_molecule():
    return API.choose_molecule()


@app.route("/chosenmolecule")
def return_chosen_molecules():
    return API.return_chosen_molecules()


@app.route("/save", methods=['GET', 'POST'])
def save_molecule():
    return API.save_molecule()


@app.route("/savedmolecules")
def return_saved_molecules():
    return API.return_saved_molecules()


@app.route("/r-group-<string:r_group_id>")
def rgroup_img(r_group_id):
    return API.rgroup_img(r_group_id)

@app.route("/molecule")
def molecule_img():
    return API.molecule_img()

@app.route("/getplotdata")
def return_assayed_data():
    return API.return_assayed_data()


@app.route("/getspiderdata")
def return_spider_data():
    return API.return_spider_data()


@app.route("/comparisontxt")
def comparison_txt():
    return API.comparison_txt()

@app.route("/reset")
def reset():
    return API.reset()


@app.route("/users/authenticate", methods=['POST'])
def authenticate_login():
    return API.authenticate_login()