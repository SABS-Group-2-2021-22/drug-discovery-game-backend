from flask import jsonify, request
from src.user import User


def authenticate_login():
    request_data = request.get_json()
    user = User(request_data['username'])
    return jsonify(request_data), user