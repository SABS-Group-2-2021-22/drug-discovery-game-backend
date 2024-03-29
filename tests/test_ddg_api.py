# import unittest
# import api.service
# import json


# class TestAPI(unittest.TestCase):
#     # def test_hello_world(self):
#     #     with api.service.app.test_client() as c:
#     #         rv = c.get('/')
#     #         json_data = rv.get_data()
#     #         result = b"<p>Hello, World!</p>"
#     #         self.assertEqual(json_data, result)

#     def test_get_all_mol_info_filled(self):
#         with api.service.app.test_client() as c:
#             c.get('/reset')
#             c.post('/save?r1=A01&r2=B01')
#             c.get('/assays?r1=A01&r2=B01&pic50=Yes&clearance_mouse=No&'
#                   'clearance_human=Yes&logd=No&pampa=Yes')
#             rv = c.get('/get_all_mol_info')
#             json_data = rv.get_json()
#             result = {'A01B01': {
#                 'assays': {
#                     'clearance_human': 'low (< 12)',
#                     'pampa': 'low', 'pic50': '6.5'},
#                 'keys': ['A01', 'B01']
#                     }
#                     }
#             self.assertEqual(json_data, result)

#     def test_get_all_mol_info_empty(self):
#         with api.app.test_client() as c:
#             c.get('/reset')
#             c.post('/save?r1=A01&r2=B01')
#             c.get('/assays?r1=A01&r2=B01&pic50=Yes&clearance_mouse=No&'
#                   'clearance_human=Yes&logd=No&pampa=Yes')
#             c.get('/reset')
#             rv = c.get('/get_all_mol_info')
#             json_data = rv.get_json()
#             result = {}
#             self.assertEqual(json_data, result)

#     def test_get_all_mol_info_empty_again(self):
#         with api.app.test_client() as c:
#             c.get('/reset')
#             rv = c.get('/get_all_mol_info')
#             json_data = rv.get_json()
#             result = {}
#             self.assertEqual(json_data, result)

#     def test_update_time_money_changed(self):
#         with api.app.test_client() as c:
#             c.get('/reset')
#             c.post('/save?r1=A01&r2=B01')
#             c.get('/assays?r1=A01&r2=B01&pic50=Yes&clearance_mouse=No&'
#                   'clearance_human=Yes&logd=No&pampa=Yes')
#             rv = c.get('/update_time_money')
#             json_data = rv.get_json()
#             result = [90230.0, 26.5]
#             self.assertEqual(json_data, result)

#     def test_update_time_money_unchanged(self):
#         with api.app.test_client() as c:
#             c.get('/reset')
#             rv = c.get('/update_time_money')
#             json_data = rv.get_json()
#             result = [100000.0, 30.0]
#             self.assertEqual(json_data, result)

#     def test_tuple2str(self):
#         true_result = "AO1B01"
#         test_result = api.tuple2str(("AO1", "B01"))
#         self.assertEqual(true_result, test_result)

#     def test_run_lipinski(self):
#         with api.app.test_client() as c:
#             c.get('/reset')
#             c.post('/save?r1=A01&r2=B01')
#             rv = c.get('/lipinski?r1=A01&r2=B01')
#             json_data = rv.get_json()
#             result = {'lipinski': {'A01B01': {
#                 'MW': True,
#                 'h_acc': True,
#                 'h_don': True,
#                 'logP': True
#                 }}}
#             self.assertEqual(json_data, result)
#             rv = c.get('/lipinski?r1=A01&r2=B01')
#             json_data = rv.get_json()
#             self.assertEqual(json_data, result)

#     def test_run_assays(self):
#         with api.app.test_client() as c:
#             c.get('/reset')
#             c.post('/save?r1=A01&r2=B01')
#             rv = c.get('/assays?r1=A01&r2=B01&pic50=Yes&clearance_mouse=No&'
#                        'clearance_human=Yes&logd=No&pampa=Yes')
#             json_data = rv.get_json()
#             result = {'assays': {'A01B01': {
#                 'clearance_human': 'low (< 12)',
#                 'pampa': 'low',
#                 'pic50': '6.5'
#                 }}}
#             self.assertEqual(json_data, result)
#             rv = c.get('/assays?r1=A01&r2=B01&pic50=Yes&clearance_mouse=Yes&'
#                        'clearance_human=Yes&logd=No&pampa=Yes')
#             json_data = rv.get_json()
#             result = {'assays': {'A01B01': {
#                 'clearance_human': 'low (< 12)',
#                 'clearance_mouse': 'medium (5.6-30.5)',
#                 'pampa': 'low',
#                 'pic50': '6.5'}}}
#             self.assertEqual(json_data, result)

#     def test_run_descriptors(self):
#         with api.app.test_client() as c:
#             c.get('/reset')
#             c.post('/save?r1=A01&r2=B01')
#             rv = c.get('/descriptors?r1=A01&r2=B01')
#             json_data = rv.get_json()
#             result = {'descriptors': {'A01B01': {
#                 'HA': 28,
#                 'MW': 397.09839370800006,
#                 'TPSA': 103.7,
#                 'h_acc': 4,
#                 'h_don': 3,
#                 'logP': 3.033400000000001,
#                 'rings': 3
#                 }}}
#             self.assertEqual(json_data, result)
#             rv = c.get('/descriptors?r1=A01&r2=B01')
#             json_data = rv.get_json()
#             self.assertEqual(json_data, result)

#     def test_run_filters(self):
#         with api.app.test_client() as c:
#             c.get('/reset')
#             c.post('/save?r1=A01&r2=B01')
#             rv = c.get('/filters?r1=A01&r2=B01')
#             json_data = rv.get_json()
#             result = {'filter_dict': {'A01B01': {
#                 'BRENK': 'passing',
#                 'NIH': 'passing',
#                 'PAINS': 'passing',
#                 'ZINC': 'passing'
#                 }}}
#             self.assertEqual(json_data, result)
#             rv = c.get('/filters?r1=A01&r2=B01')
#             json_data = rv.get_json()
#             self.assertEqual(json_data, result)

#     def test_choose_molecule(self):
#         with api.app.test_client() as c:
#             c.get('/reset')
#             rv = c.post('/choose?r1=A01&r2=B10')
#             json_data = rv.get_json()
#             result = {'chosen_mol': ['A01', 'B10']}
#             self.assertEqual(json_data, result)

#     def test_return_chosen_molecules(self):
#         with api.app.test_client() as c:
#             c.get('/reset')
#             c.post('/choose?r1=A01&r2=B10')
#             rv = c.get('chosenmolecule')
#             json_data = rv.get_json()
#             result = {'chosen_mol': ['A01', 'B10']}
#             self.assertEqual(json_data, result)
#             c.get('/reset')
#             rv = c.get('chosenmolecule')
#             json_data = rv.get_json()
#             result = {}
#             self.assertEqual(json_data, result)

#     def test_save_molecule(self):
#         with api.app.test_client() as c:
#             c.get('/reset')
#             rv = c.post('/save?r1=A04&r2=B05')
#             json_data = rv.get_json()
#             result = {'saved_mols': [['A04', 'B05']]}
#             self.assertEqual(json_data, result)

#     def test_return_saved_molecules(self):
#         with api.app.test_client() as c:
#             c.get('/reset')
#             c.post('/save?r1=A04&r2=B05')
#             rv = c.get('/savedmolecules')
#             json_data = rv.get_json()
#             result = {'saved_mols': [['A04', 'B05']]}
#             self.assertEqual(json_data, result)

#     def test_rgroup_img(self):
#         with api.app.test_client() as c:
#             rv = c.get('/r-group-B10')
#             json_data = rv.get_json()
#             print(json_data)
#             with open('tests/test_rgroup_img.txt') as true_file:
#                 true_data = json.load(true_file)
#             self.assertEqual(json_data, true_data)

#     def test_molecule_img(self):
#         with api.app.test_client() as c:
#             c.get('/reset')
#             rv = c.get('/molecule?r1=A04&r2=B04&size=800,800')
#             json_data = rv.get_json()
#             with open('tests/test_molecule_full.txt') as true_file:
#                 true_data = json.load(true_file)
#             self.assertEqual(json_data, true_data)

#     def test_reset(self):
#         with api.app.test_client() as c:
#             rv = c.get('/reset')
#             json_data = rv.get_json()
#             result = {"new_info": {}}
#             self.assertEqual(json_data, result)

#     def test_return_assayed_data(self):
#         with api.app.test_client() as c:
#             c.get('/reset')
#             c.post('/save?r1=A01&r2=B01')
#             c.get('/descriptors?r1=A01&r2=B01')
#             c.get('/assays?r1=A01&r2=B01&pic50=Yes&clearance_mouse=Yes&'
#                   'clearance_human=Yes&logd=Yes&pampa=Yes')
#             rv = c.get('/getplotdata')
#             json_data = rv.get_json()
#             result = {"assay_dict": [{'A01B01': {
#                         '--': 0,
#                         'HA': 28,
#                         'MW': 397.09839370800006,
#                         'TPSA': 103.7,
#                         'clearance_human': "low (< 12)",
#                         'clearance_mouse': "medium (5.6-30.5)",
#                         'h_acc': 4,
#                         'h_don': 3,
#                         'logP': 3.033400000000001,
#                         'logd': 0.3,
#                         'pampa': "low",
#                         'pic50': "6.5",
#                         'rings': 3
#             }}]}
#             self.assertDictEqual(json_data, result)

#     def test_numerise_params(self):
#         with api.app.test_client() as c:
#             c.get('/reset')
#             c.post('/choose?r1=A02&r2=B02')
#             rv = c.get('/getspiderdata')
#             json_data = rv.get_json()
#             result = {'param_dict': [{
#                 'clearance_human': 1,
#                 'clearance_mouse': 1,
#                 'logd': 0.54,
#                 'pampa': 5.5,
#                 'pic50': "Not Made"
#             },
#                 {
#                 'clearance_human': 1,
#                 'clearance_mouse': 1,
#                 'logd': 1.08,
#                 'pampa': 5.5,
#                 'pic50': "7.7"
#             }]}
#             self.assertDictEqual(json_data, result)

#     def test_comparison_txt(self):
#         with api.app.test_client() as c:
#             c.get('/reset')
#             c.post('/choose?r1=A03&r2=B03')
#             rv = c.get('/comparisontxt')
#             json_data = rv.get_json()
#             self.assertIsNotNone(json_data)
#             comp_dict = json_data['comparison']
#             with open('./tests/comp_txt_test.txt', 'r') as f:
#                 lines = f.readlines()
#                 self.assertEqual(str(lines[0]), comp_dict['pic50'])
#                 self.assertEqual(str(lines[1]), comp_dict['logd'])
#                 self.assertEqual(str(lines[2]), comp_dict['clearance_human'])
