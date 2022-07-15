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
        print(true_user)
        print(test_user)
        self.assertEqual(test_user, true_user)
