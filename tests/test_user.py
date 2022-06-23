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

    def test_get_and_update_molecule_info(self):
        test_user = User('test')
        true_molecule_info = test_user.get_molecule_info()
        self.assertEqual(true_molecule_info, {})
        test_user.update_molecule_info({'test': 42})
        new_molecule_info = test_user.get_molecule_info()
        self.assertEqual(new_molecule_info, {'test': 42})

