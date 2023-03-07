from src.user import User


def authenticate_login(request_data):
    """ Creates new instance of User object using username sent from the
    frontend as API call.

    :return: The request data and the new User object
    :rtype: json dict, User object
    """
    user = User(request_data['username'])
    return request_data, user
