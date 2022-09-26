from re import A
import unittest
import api.backend as api
from src.user import User


class TestUserAPI(unittest.TestCase):
    """Tests for api.backed.user_api functions
    """
    def test_authenticate_login(self):
        test_data = {'username': 'User'}
        test_request_data, test_user = api.authenticate_login(test_data)
        self.assertEqual(test_request_data, test_data)
        true_user = User('User')
        self.assertEqual(test_user.chosen_mol, true_user.chosen_mol)
        self.assertEqual(test_user.username, true_user.username)
        self.assertEqual(test_user.money, true_user.money)
        self.assertEqual(test_user.time, true_user.time)
        self.assertEqual(test_user.molecule_info, true_user.molecule_info)
