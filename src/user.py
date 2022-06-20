
import json
import datetime


class User(object):
    def __init__(self, username):
        self.username = username

        self.money = 100000.0
        self.time = 30.0
        self.chosen_mol = [None, None]
        self.molecule_info = {}

    def __repr__(self):
        return 'class instance User\n' + \
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

    def save_game(self):
        user_as_dict = {self.username: {
           'money': self.money,
           'time': self.time,
           'chosen_mol': self.chosen_mol,
           'molecule_info': self.molecule_info
                                    }}
        filename = 'out/' + self.username + '_'\
            + datetime.datetime.now().strftime('%Y%m%d_%H%M') + '.json'
        with open(filename, 'w+') as fp:
            json.dump(user_as_dict, fp)
