import base64
import json
from src.user import User
from glob import glob
from flask import Flask, jsonify, request, send_from_directory, send_file, Response
from flask_cors import CORS
from flask_cors import cross_origin
from http.server import HTTPServer, SimpleHTTPRequestHandler, test


# import api.backend.backend_api as api
import api.backend as api
import sys

import os



app = Flask(__name__)

cors = CORS(
    app,
    resources={r"/*": {"origins":
                       "*" 
                    #    "http://localhost:3000"
                       }},
)
class CORSRequestHandler (SimpleHTTPRequestHandler):
    def end_headers (self):
        self.send_header('Access-Control-Allow-Origin', '*')
        SimpleHTTPRequestHandler.end_headers(self)

if __name__ == '__main__':
    test(CORSRequestHandler, HTTPServer, port=int(sys.argv[1]) if len(sys.argv) > 1 else 8000)



global sessions
sessions = {}

@app.route("/load_sessions")
def load_sessions():
    session_files = glob('src/saved_data/*.json')
    for file in session_files:
        username = file.split('/')[-1].split('_')[0]
        with open(file) as f:
            session_dict = json.load(f)
            sessions[username] = User(username, session_dict[username])
    return None

@app.route("/reset_and_clear_sessions")
def clear_sessions():
    sessions = {}
    return None


# It is likely that Nele will hate this function
def check_user(username):
    if username not in sessions:
        sessions[username] = User(username)
    return None


@app.route("/lipinski")
def run_lipinski():
    """API function for running run_lipinski() function.

    Pass R Group IDs as queries: /lipinski?r1=A01&r2=B01.

    :return: A json dictionary of the molecule, indexed
    by the concatenated string of its R Group IDs, with the values being
    True or False for each descriptor relevant to Lipinski's Rule of 5.
    :rtype: json dict
    """
    username = json.loads(request.headers['username'])['username']
    check_user(username)
    session_molecule_info = sessions[username].get_molecule_info()
    response, updated_mol_dict = api.run_lipinski(session_molecule_info)
    sessions[username].update_molecule_info(updated_mol_dict)
    return response


@app.route("/descriptors")
def run_descriptors():
    """API call for running run_descriptors() function.

    Pass R Group IDs as queries: /descriptors?r1=A01&r2=B01.

    :return: A json dictionary of the molecule, indexed
    by the concatenated string of its R Group IDs, with the values for each
    descriptor, with each key being its respective descriptor label.
    :rtype: json dict
    """
    username = json.loads(request.headers['username'])['username']
    check_user(username)
    session_molecule_info = sessions[username].get_molecule_info()
    response, updated_mol_dict = api.run_descriptors(session_molecule_info)
    sessions[username].update_molecule_info(updated_mol_dict)
    return response


@app.route("/choose", methods=["GET", "POST"])
def choose_molecule():
    """API call for running choose_molecule() function.

    Pass R group IDs as queries: /choose?r1=A01&r2=B10.

    :return: Json dict contained list of the R Group IDs and the list on its
    own.
    :rtype: json dict
    """
    username = json.loads(request.headers['username'])['username']
    check_user(username)
    response, new_chosen_molecule = api.choose_molecule()
    sessions[username].set_chosen_molecule(new_chosen_molecule)
    return response


@app.route("/save", methods=["GET", "POST"])
def save_molecule():
    """API call for running save_molecule() function.

    :return: List of tuples, containing the R Group IDs as a json dict
    :rtype: json dict
    """
    username = json.loads(request.headers['username'])['username']
    check_user(username)
    session_molecule_info = sessions[username].get_molecule_info()
    response, updated_mol_dict = api.save_molecule(session_molecule_info)
    sessions[username].update_molecule_info(updated_mol_dict)
    return response

@app.route("/r-group")
def rgroup_img():
    """API call for running rgroup_img() function.

    Pass R Group ID as query: /r-group-A01.

    :param r_group_id: ID number of R Group, eg. 'B26'
    :type r_group_id: String
    :return: Image and stats of R Group in a json dict.
    Access image bytestream with `img_html` key and stats with 'stats'
    :rtype: json dict
    """
    return api.rgroup_img()


@app.route("/molecule")
def molecule_img():
    """API call for running molecule_img() function.

    Pass R Group IDs and image size as queries:
    /molecule?r1=A01&r2=B10&size=800,800

    :return: JSON containing image of compound molecule as bytestream, drug
    properties and empty dictionaries for descriptors, filters and assays run.
    :rtype: json dict
    """
    return api.molecule_img()


@app.route("/getspiderdata")
def return_spider_data():
    """API call for running return_spider_data() function.

    Call /getspiderdata.

    :return: A json dictionary containing a list of 2 dictionaries, one
    containing chosen mol parameters and the other containing reference
    drug parameters.
    :rtype: json dict
    """
    username = json.loads(request.headers['username'])['username']
    check_user(username)
    session_chosen_mol = sessions[username].get_chosen_molecule()
    return api.return_spider_data(session_chosen_mol)


@app.route("/comparisontxt")
def comparison_txt():
    """API call for running comparison_text() function.

    Call /comparisontxt.

    :return: json dict with text in value depending on metric
    :rtype: json dict
    """
    username = json.loads(request.headers['username'])['username']
    check_user(username)
    session_chosen_mol = sessions[username].get_chosen_molecule()
    return api.comparison_txt(session_chosen_mol)

@app.route("/infotxt")
def info_text():
    """API call for running info_text() function.

    Call /infontxt.

    :return: json dict with text in value depending on help button
    :rtype: json dict
    """
    return api.info_text()


@app.route("/reset")
def reset():
    """API call for running reset() function.

    Call /reset.

    :return: The new empty molecule_info dictionary as a json dict
    :rtype: json dict
    """
    username = json.loads(request.headers['username'])['username']
    check_user(username)
    sessions[username].update_molecule_info({})
    session_molecule_info = sessions[username].get_molecule_info()
    return api.reset(session_molecule_info)


@app.route("/authenticate", methods=["POST"])
def authenticate_login():
    """API call for running authenticate_login() function.

    Call /users/authenticate.

    :return: The authentication response.
    :rtype: json dict
    """
    request_data = request.get_json()
    non_jsonified_auth_response, user = api.authenticate_login(request_data)
    if user.username not in sessions:
        sessions[user.username] = user
    else:
        non_jsonified_auth_response['user_status'] = 'Exists'
    # sessions[user.username] = user
    auth_response = jsonify(non_jsonified_auth_response)
    return auth_response

@app.route("/logout", methods=["POST"])
def logout():
    """API call for running authenticate_login() function.

    Call /users/authenticate.

    :return: The authentication response.
    :rtype: json dict
    """
    username = json.loads(request.headers['username'])['username']
    # check_user(username)
    # sessions[username].save_game()
    return jsonify({})

@app.route("/loadgamestate", methods=["GET"])
def loadgamestate():
    #todo: add comments 
    
    username = json.loads(request.headers['username'])['username']
    # check_user(username)
    # sessions[username].save_game()
    filename = 'src/saved_data/' + username + '.json'
    try:
        with open(filename, 'r') as f:
            game_data = json.load(f)
            sessions[username] = User(username, game_data[username])
        return jsonify(game_data)
    except:
        return jsonify({'error':'There is an error retrieving your data from the previous game'}) 

# TODO api.save_game_data() does not exist so currently returning nothing


@app.route("/save_game_data", methods=["GET", "POST"])
def save_game_data():
    """API call for running save_game_data() (currently does not exist).

    Call /save_game_data.

    :return: _description_
    :rtype: _type_
    """
    username = json.loads(request.headers['username'])['username']
    game_data = request.json
    check_user(username)
    sessions[username].save_game(game_data)
    return username


@app.route("/sketcher_save_molecule")
def sketcher_save_molecule():
    """API call for sketcher_save_molecule() function.

    Call /sketcher_save_molecule?mol=<BYTESTREAM>

    :return: A json dictionary of the bytestream, descriptors, filters and
    Lipinksi rules for the sketched molecule.
    :rtype: json dict
    """
    username = json.loads(request.headers['username'])['username']
    check_user(username)
    session_molecule_info = sessions[username].get_molecule_info()
    mol_block = base64.b64decode(request.args.get("mol"))
    response, updated_mol_dict = api.sketcher_save_molecule(
        mol_block, session_molecule_info
    )
    sessions[username].update_molecule_info(updated_mol_dict)
    return response


@app.route("/sketcher_choose", methods=["POST"])
def sketcher_choose():
    """API call for sketcher_choose_molecule() function.

    Call /sketcher_choose?id=<ID>&smiles=<SMILES>.

    :return: List of the ID and smiles of the chosen sketcher molecule as a
    json dict. Access list with 'chosen_mol' key.
    :rtype: json dict
    """
    username = json.loads(request.headers['username'])['username']
    check_user(username)
    response, new_chosen_molecule = api.sketcher_choose_molecule()
    sessions[username].set_chosen_molecule(new_chosen_molecule)
    return response


@app.route("/sketcher_comparisontxt")
def sketcher_comparisontxt():
    """API call for sketcher_comparison_txt() function.

    Call /sketcher_comparisontxt.

    :return: Json dict with comparison text with values for each metric.
    :rtype: json dict
    """
    username = json.loads(request.headers['username'])['username']
    check_user(username)
    session_chosen_mol = sessions[username].get_chosen_molecule()
    return api.sketcher_comparison_txt(session_chosen_mol)


@app.route("/sketcher_getspiderdata")
def sketcher_getspiderdata():
    """API call for sketcher_return_spider_data() function.

    Call /sketcher_getspiderdata.

    :return: A json dictionary containing a list of 2 dictionaries, one
     containing chosen mol.parameters and the other containing reference
     drug parameters
    :rtype: json dict
    """
    username = json.loads(request.headers['username'])['username']
    check_user(username)
    session_chosen_mol = sessions[username].get_chosen_molecule()
    return api.sketcher_return_spider_data(session_chosen_mol)

@app.route("/docking-<path:filename>")
def serve_pdb_file(filename):
    # def serve_pdb_file(pdb_filename):
        import os
        print('CP1', os.getcwd())
        filepath = './static/ligand_docks/' + filename
        with open(filepath, 'rb') as file:
                byte_stream = file.read()
        base64_encoded = base64.b64encode(byte_stream).decode('utf-8')
        print(byte_stream)
        return Response(byte_stream, mimetype='application/octet-stream', headers={"Content-Disposition": "attachment;filename=" + filename})

# def convert_pdb_files_to_byte_stream():
#     directory = 'path/to/your/directory'
#     encoded_files = {}

#     for filename in os.listdir(directory):
#         if filename.endswith('.pdb'):
#             pdb_file_path = os.path.join(directory, filename)
#             with open(pdb_file_path, 'rb') as file:
#                 byte_stream = file.read()
#             encoded_files[filename] = base64.b64encode(byte_stream).decode('utf-8')

#     return jsonify(encoded_files)
#    """API call for running serve_pdb_file() function.

#     Pass PDB filename as query: /docking-pdb_filename.

#     :param pdb_filename: PDB file of docked molecule, eg. '6lu7.pdb'
#     :type r_group_id: String
#     :return: Image and stats of R Group in a json dict.
#     Access image bytestream with `img_html` key and stats with 'stats'
#     :rtype: json dict
#     """