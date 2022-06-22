import unittest
from src.user import User


class TestUser(unittest.TestCase):
    """Tests for User class"""
    def test___repr__(self):
        """Tests that __repr__ of User class outputs correct string
        """
        test_user = User('test')
        test_repr = repr(test_user)
        true_repr = 'class instance User\n' + \
               'User: test \n' + \
               'Money: 100000.0 \n' + \
               'Time: 30.0 \n' + \
               'Chosen Mol: [None, None] \n' + \
               'Molecule Info: {} \n'
        self.assertEqual(test_repr, true_repr)

    # def test_get_molecule_info(self):
