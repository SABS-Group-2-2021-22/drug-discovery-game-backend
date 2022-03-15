

class User(object):
    def __init__(self, username):
        self.username = username

        self.money = 100000.0
        self.time = 30.0
        self.chosen_mol = [None, None]
        self.molecule_info = {}

    def __repr__(self):
        return f'class instance User\n' + \
               f'User: {self.username} \n' + \
               f'Money: {self.money} \n' + \
               f'Time: {self.time} \n' + \
               f'Chosen Mol: {str(self.chosen_mol)} \n' + \
               f'Molecule Info: {str(self.molecule_info)} \n'
             
    def get_molecule_info(self):
        return self.molecule_info

    def get_chosen_molecule(self):
        return self.chosen_mol
    
    def set_chosen_molecule(self, chosen_mol):
        self.chosen_mol = chosen_mol
        return chosen_mol

    def update_molecule_info(self, update_dict):
        self.molecule_info = update_dict
        return None
