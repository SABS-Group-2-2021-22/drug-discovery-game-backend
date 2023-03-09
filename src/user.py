
import json
import datetime


class User(object):
    """User class that defines the users of the game. Stores the game information.
    """
    def __init__(self, username):
        """Initiliases User object using new username.

        :param username: Username
        :type username: str
        """
        self.username = username

        self.money = 100000.0
        self.time = 30.0
        self.chosen_mol = [None, None]
        self.molecule_info = {}

    def __repr__(self):
        """Represents User object and all its properties as a string

        :return: String representation of User
        :rtype: str
        """
        return 'class instance User\n' + \
               f'User: {self.username} \n' + \
               f'Money: {self.money} \n' + \
               f'Time: {self.time} \n' + \
               f'Chosen Mol: {str(self.chosen_mol)} \n' + \
               f'Molecule Info: {str(self.molecule_info)} \n'

    def get_molecule_info(self):
        """Returns molecules saved and their info during a game

        :return: Molecule information
        :rtype: dict
        """
        return self.molecule_info

    def get_chosen_molecule(self):
        """Returns chosen molecule as final molecule

        :return: Chosen molecule as a list of keys
        :rtype: list
        """
        return self.chosen_mol

    def set_chosen_molecule(self, chosen_mol):
        """Set chosen molecule in object as chosen_mol

        :param chosen_mol: The chosen mol, sent from the front-end
        :type chosen_mol: list
        :return: The chosen molecule
        :rtype: list
        """
        self.chosen_mol = chosen_mol
        return chosen_mol

    def update_molecule_info(self, update_dict):
        """Updated molecule_info with a new dictionary, update_dict

        :param update_dict: The updated saved molecules and their information
        :type update_dict: dict
        :return: None
        :rtype: None
        """
        self.molecule_info = update_dict
        return None

    def save_game(self):
        """Saves user information as a dictionaruy and saves to file. File is
        named after user and the time of saving.
        """
        user_as_dict = {self.username: {
           'money': self.money,
           'time': self.time,
           'chosen_mol': self.chosen_mol,
           'molecule_info': self.molecule_info
                                    }}
        filename = 'src/saved_data/' + self.username + '_'\
            + datetime.datetime.now().strftime('%Y%m%d_%H%M') + '.json'
        with open(filename, 'w+') as fp:
            json.dump(user_as_dict, fp)
