from flask import Flask
from flask_cors import CORS

app = Flask(__name__)
cors = CORS(app, resources={r"/api/*": {"origins": "*"}})   

@app.route('/test', methods=['GET'])
def test():
    return 'Test route is working', 200


