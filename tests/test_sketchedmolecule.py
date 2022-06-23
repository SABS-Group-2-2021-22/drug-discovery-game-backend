import unittest
from src.sketchedMolecule import sketchedMolecule

TEST_MOL_BLOCK = """
 OpenBabel06232214542D

 11 11  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  6 11  1  0  0  0  0
  6  7  2  0  0  0  0
  7  8  1  0  0  0  0
  8  9  2  0  0  0  0
  9 10  1  0  0  0  0
 10 11  2  0  0  0  0
M  END
"""


class TestSketcherMolecule(unittest.TestCase):
    """Tests for sketchedMolecule class
    """
    def test_smiles(self):
        """Tests that smiles() function returns correct smiles
        """
        test_molecule = sketchedMolecule(TEST_MOL_BLOCK)
        smiles = test_molecule.smiles
        self.assertEqual(smiles, "CCCCOc1ccccc1")

    def test_drawMoleculeAsByteStream(self):
        """Tests that molecule images are drawn as byte-streams, does not
        check image is correct though
        """
        test_molecule = sketchedMolecule(TEST_MOL_BLOCK)
        for size in None, 100:
            test_byte_stream =\
                 test_molecule.drawMoleculeAsByteStream(size=size)
            self.assertIsInstance(test_byte_stream, str)
            self.assertGreater(len(test_byte_stream), 1000)

    def test_filter_properties(self):
        """Checks that PAINS and other filters are correct in passing or
        failing a molecule
        """
        test_molecule = sketchedMolecule(TEST_MOL_BLOCK)
        test_filters = test_molecule.filter_properties()
        true_filters = {'PAINS': 'passing',
                        'ZINC': 'passing',
                        'BRENK': ['Aliphatic_long_chain'],
                        'NIH': 'passing'}
        self.assertEqual(test_filters, true_filters)

    def test_descriptors(self):
        """Tests that descriptors for a molecule are calculated correctly
        """
        test_molecule = sketchedMolecule(TEST_MOL_BLOCK)
        test_desc = test_molecule.descriptors()
        true_desc = {
            'mol': "CCCCOc1ccccc1",
            'MW': 150.104465068,
            'logP': 2.8655000000000017,
            'TPSA': 9.23,
            'HA': 11,
            'h_acc': 1,
            'h_don': 0,
            'rings': 1}
        self.assertEqual(test_desc, true_desc)

    def test_lipinski(self):
        """Ensures function is able to check if a molecule passes the Lipinski
        rules are not
        """
        test_molecule = sketchedMolecule(TEST_MOL_BLOCK)
        test_lipinski = test_molecule.lipinski()
        true_lipinski = {
            'MW': True,
            'h_acc': True,
            'h_don': True,
            'logP': True,
            }
        self.assertEqual(test_lipinski, true_lipinski)

    def test_drug_properties(self):
        """Checks that function is able to predict values for the assay values
        of a molecule. This should be updated if the models are retrained.
        """
        test_molecule = sketchedMolecule(TEST_MOL_BLOCK)
        test_drug_props = test_molecule.drug_properties()
        true_drug_props = {
            'pic50': '4.6772484184733205',
            'clearance_mouse': 'medium (5.6-30.5)',
            'clearance_human': 'low (< 12)',
            'logd': '0.11891286727456941',
            'pampa': 'med2high',
        }
        self.assertEqual(test_drug_props, true_drug_props)
