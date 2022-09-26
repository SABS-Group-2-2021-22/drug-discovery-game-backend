from api.backend.utils import tuple2str, numerise_params
import unittest


class TestUtils(unittest.TestCase):
    """Tests for api.backed.utils functions
    """
    def test_tuple2str(self):
        test_tuple = ("AC", "BD")
        test_str = tuple2str(test_tuple)
        true_str = "ACBD"
        self.assertEqual(test_str, true_str)

    def test_numerise_params(self):
        input_dict = {
            'clearance_human': 'fair',
            'clearance_mouse': 'poor',
            'pampa': 'low',
            'logd': 'best'
        }
        true_dict = {
            'clearance_human': 4,
            'clearance_mouse': 7,
            'pampa': 2.5,
            'logd': 8,
        }
        test_dict = numerise_params(input_dict)
        self.assertEqual(true_dict, test_dict)
