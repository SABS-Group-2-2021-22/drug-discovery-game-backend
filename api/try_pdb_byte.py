import os
import base64
from flask import jsonify

def convert_pdb_files_to_byte_stream(directory):
    encoded_files = {}  # Dictionary to store encoded file content

    for filename in os.listdir(directory):
        if filename.endswith('.pdb'):
            pdb_file_path = os.path.join(directory, filename)
            with open(pdb_file_path, 'rb') as file:
                byte_stream = file.read()
            base64_encoded = base64.b64encode(byte_stream).decode('utf-8')
            encoded_files[filename] = base64_encoded  # Add to dictionary
            print(f"Converted {filename} to byte stream and encoded to base64.")

    return jsonify(encoded_files)

# Example usage (assuming you're in a Flask route)
directory = '/Users/sanazkazeminia/Drug_disc_game/drug-discovery-game-backend/static/ligand_docks'
response = convert_pdb_files_to_byte_stream(directory)
