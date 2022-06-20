from flask import Flask, request
from flask_cors import CORS
import base64
import json


# import api.backend.backend_api as api
import api.backend as api

app = Flask(__name__)

cors = CORS(app,
            resources={r"/*": {"origins": "http://localhost:3000"}},
            )

global sessions
sessions = {}


# TODO: refactor query passing here, direct via function args
@app.route("/lipinski")
def run_lipinski():
    username = json.loads(request.headers['username'])['username']
    session_molecule_info = sessions[username].get_molecule_info()
    response, updated_mol_dict = api.run_lipinski(session_molecule_info)
    sessions[username].update_molecule_info(updated_mol_dict)
    return response


@app.route("/descriptors")
def run_descriptors():
    username = json.loads(request.headers['username'])['username']
    session_molecule_info = sessions[username].get_molecule_info()
    response, updated_mol_dict = api.run_descriptors(session_molecule_info)
    sessions[username].update_molecule_info(updated_mol_dict)
    return response


@app.route("/choose", methods=['GET', 'POST'])
def choose_molecule():
    username = json.loads(request.headers['username'])['username']
    response, new_chosen_molecule = api.choose_molecule()
    sessions[username].set_chosen_molecule(new_chosen_molecule)
    return response


@app.route("/save", methods=['GET', 'POST'])
def save_molecule():
    username = json.loads(request.headers['username'])['username']
    session_molecule_info = sessions[username].get_molecule_info()
    response, updated_mol_dict = api.save_molecule(session_molecule_info)
    sessions[username].update_molecule_info(updated_mol_dict)
    return response


@app.route("/r-group-<string:r_group_id>")
def rgroup_img(r_group_id):
    return api.rgroup_img(r_group_id)


@app.route("/molecule")
def molecule_img():
    return api.molecule_img()


@app.route("/getspiderdata")
def return_spider_data():
    username = json.loads(request.headers['username'])['username']
    session_chosen_mol = sessions[username].get_chosen_molecule()
    return api.return_spider_data(session_chosen_mol)


@app.route("/comparisontxt")
def comparison_txt():
    username = json.loads(request.headers['username'])['username']
    session_chosen_mol = sessions[username].get_chosen_molecule()
    return api.comparison_txt(session_chosen_mol)


@app.route("/reset")
def reset():
    return api.reset()


@app.route("/users/authenticate", methods=['POST'])
def authenticate_login():
    auth_response, user = api.authenticate_login()
    if user.username not in sessions:
        sessions[user.username] = user
    return auth_response


@app.route("/save_game_data", methods=['GET'])
def save_game_data():
    username = json.loads(request.headers['username'])['username']
    sessions[username].save_game()
    return api.save_game_data()


@app.route("/sketcher_save_molecule")
def sketcher_save_molecule():
    username = json.loads(request.headers['username'])['username']
    session_molecule_info = sessions[username].get_molecule_info()
    mol_block = base64.b64decode(request.args.get('mol'))
    response, updated_mol_dict = api.sketcher_save_molecule(
        mol_block, session_molecule_info)
    sessions[username].update_molecule_info(updated_mol_dict)
    print(response)
    return response


@app.route("/sketcher_choose", methods=['POST'])
def sketcher_choose():
    username = json.loads(request.headers['username'])['username']
    response, new_chosen_molecule = api.sketcher_choose_molecule()
    sessions[username].set_chosen_molecule(new_chosen_molecule)
    return response


@app.route("/sketcher_comparisontxt")
def sketcher_comparisontxt():
    username = json.loads(request.headers['username'])['username']
    session_chosen_mol = sessions[username].get_chosen_molecule()
    return api.sketcher_comparison_txt(session_chosen_mol)


@app.route("/sketcher_getspiderdata")
def sketcher_getspiderdata():
    username = json.loads(request.headers['username'])['username']
    session_chosen_mol = sessions[username].get_chosen_molecule()
    return api.sketcher_return_spider_data(session_chosen_mol)
