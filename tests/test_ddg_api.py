import unittest
import ddg_api


class TestAPI(unittest.TestCase):
    def test_hello_world(self):
        with ddg_api.app.test_client() as c:
            rv = c.get('/')
            json_data = rv.get_data()
            result = b"<p>Hello, World!</p>"
            self.assertEqual(json_data, result)

    def test_get_all_mol_info_filled(self):
        with ddg_api.app.test_client() as c:
            c.post('/save?r1=A01&r2=B01')
            c.get('/assays?r1=A01&r2=B01&pic50=Yes&clearance_mouse=No&'
                  'clearance_human=Yes&logd=No&pampa=Yes')
            rv = c.get('/get_all_mol_info')
            json_data = rv.get_json()
            result = {'A01B01': {
                'assays': {
                    'clearance_human': 'low (< 12)',
                    'pampa': 'low', 'pic50': '6.5'},
                'keys': ['A01', 'B01']
                    }
                    }
            self.assertEqual(json_data, result)

    def test_get_all_mol_info_empty(self):
        with ddg_api.app.test_client() as c:
            c.post('/save?r1=A01&r2=B01')
            c.get('/assays?r1=A01&r2=B01&pic50=Yes&clearance_mouse=No&'
                  'clearance_human=Yes&logd=No&pampa=Yes')
            c.get('/reset')
            rv = c.get('/get_all_mol_info')
            json_data = rv.get_json()
            result = {}
            self.assertEqual(json_data, result)

    def test_get_all_mol_info_empty_again(self):
        with ddg_api.app.test_client() as c:
            rv = c.get('/get_all_mol_info')
            json_data = rv.get_json()
            result = {}
            self.assertEqual(json_data, result)

    def test_update_time_money(self):
        with ddg_api.app.test_client() as c:
            c.post('/save?r1=A01&r2=B01')
            c.get('/assays?r1=A01&r2=B01&pic50=Yes&clearance_mouse=No&'
                  'clearance_human=Yes&logd=No&pampa=Yes')
            rv = c.get('/update_time_money')
            json_data = rv.get_json()
            result = [90230.0, 26.5]
            self.assertEqual(json_data, result)


    def test_reset(self):
        with ddg_api.app.test_client() as c:
            rv = c.get('/reset')
            json_data = rv.get_json()
            result = {"new_info": {}}
            self.assertEqual(json_data, result)
