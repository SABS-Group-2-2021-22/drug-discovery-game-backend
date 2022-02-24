from flask import Flask, jsonify, request, render_template_string
from flask_cors import CORS


# import api.backend.backend_api as api
import api.backend as api

from src.user import User

app = Flask(__name__)

cors = CORS(app, 
            resources={r"/*":{"origins": "http://localhost:3000"}},
            )

global sessions
sessions = {}

@app.route("/")
def hello_world():
    return api.hello_world()


@app.route("/get_all_mol_info")
def get_all_mol_info():
    return api.get_all_mol_info()


@app.route("/update_time_money")
def update_time_and_money():
    return api.update_time_and_money()


# TODO: refactor query passing here, direct via function args
@app.route("/lipinski")
def run_lipinski():
    return api.run_lipinski()


@app.route("/assays")
def run_assays():
    return api.run_assays()


@app.route("/descriptors")
def run_descriptors():
    return api.run_descriptors()


@app.route("/filters")
def run_filters():
    return api.run_filters()


@app.route("/choose", methods=['GET', 'POST'])
def choose_molecule():
    return api.choose_molecule()


@app.route("/chosenmolecule")
def return_chosen_molecules():
    return api.return_chosen_molecules()


@app.route("/save", methods=['GET', 'POST'])
def save_molecule():
    username = request.get_json()['username']
    session_molecule_info = sessions[username].get_molecule_info()
    response, updated_mol_dict = api.save_molecule(session_molecule_info)
    sessions[username].update_molecule_info(updated_mol_dict)
    return response


@app.route("/savedmolecules")
def return_saved_molecules():
    return api.return_saved_molecules()


@app.route("/r-group-<string:r_group_id>")
def rgroup_img(r_group_id):
    return api.rgroup_img(r_group_id)


@app.route("/molecule")
def molecule_img():
    return api.molecule_img()


@app.route("/getplotdata")
def return_assayed_data():
    return api.return_assayed_data()


@app.route("/getspiderdata")
def return_spider_data():
    return api.return_spider_data()


@app.route("/comparisontxt")
def comparison_txt():
    return api.comparison_txt()


@app.route("/reset")
def reset():
    return api.reset()


@app.route("/users/authenticate", methods=['POST'])
def authenticate_login():
    auth_response, user = api.authenticate_login()
    if user.username not in sessions:
        sessions[user.username] = user
    return auth_response
