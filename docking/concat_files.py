import os

# Set the paths to the input PDB directory and the 7leu PDB file
pdb_dir = '/Users/apunt/repos/drug-discovery-game-backend/static/ligand_docks/'
pdb_file = '/Users/apunt/repos/drug-discovery-game-backend/static/3ehy_no_ligand.pdb'

# Loop through all PDB files in the input directory
for pdb_name in os.listdir(pdb_dir):
    if pdb_name.endswith('dock1.pdb'):
        # Read in the contents of the current PDB file
        with open(os.path.join(pdb_dir, pdb_name), 'r') as f:
            pdb_contents = f.read()
        # Read in the contents of the 7leu PDB file
        with open(pdb_file, 'r') as f:
            mmp12_contents = f.read()
        # Concatenate the two PDB files
        new_contents = pdb_contents + mmp12_contents

        # Create a new filename with "_concatenated" suffix
        new_name = pdb_name.replace('.pdb', '_concatenated.pdb')

        # Write out the concatenated PDB file
        with open(os.path.join(pdb_dir, new_name), 'w') as f:
            f.write(new_contents)