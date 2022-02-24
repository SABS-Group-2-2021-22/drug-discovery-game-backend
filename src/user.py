

class User(object):
    def __init__(self, username):
        self.username = username

        self.money = 100000.0
        self.time = 30.0
        self.chosen_mol = {}
        self.molecule_info = {}

    def __repr__(self):
        return f'User: {self.username} \n' + \
               f'Money: {self.money} \n' + \
               f'Time: {self.time} \n' + \
               f'Chosen Mol: {str(self.chosen_mol)} \n' + \
               f'Molecule Info: {str(self.molecule_info)} \n'
                

    def get_molecule_info(self):
        return self.molecule_info

    def update_molecule_info(self, update_dict):
        self.molecule_info = update_dict
        return None