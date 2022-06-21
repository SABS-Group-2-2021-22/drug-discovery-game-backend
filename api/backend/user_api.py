from flask import jsonify, request
from src.user import User


def authenticate_login():
    """ Creates new instance of User object using user sent from the frontend
    as API call.

    :return: The request data and the new User object
    :rtype: json dict, User object
    """
    request_data = request.get_json()
    user = User(request_data['username'])
    return jsonify(request_data), user
