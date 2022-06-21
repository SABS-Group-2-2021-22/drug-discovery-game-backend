# TODO: move this into backend function
def tuple2str(tuple_in):
    """Converts a tuple into a string.

    :param tuple_in: tuple to convert
    :type tuple_in: tuple
    :returns: concatenated string version of the tuple
    :rtype: str
    """
    string = ''
    for i in tuple_in:
        string += str(i)
    return string


def numerise_params(prop_dict):
    """ Returns drug properties with all qualitative values transformed into
    numeric values

    returns: numerically transformed property dictionaries
    rtype: dict
    """
    clearance_dict = {
        'low (< 5.6)': 1,
        'medium (5.6-30.5)': 4,
        'low (< 3.7)': 1,
        'good': 1,
        'high (> 30.5)': 7,
        'fair': 4,
        'poor': 7,
        'low (< 12)': 1,
        'medium (12-44)': 4,
        'medium (5.6-30.5)': 4
    }
    pampa_dict = {
        'neg': 0,
        'poor': 1,
        'low': 2.5,
        'fair': 5.5,
        'med2high': 5.5,
        'good': 6.5,
        'best': 8

    }
    drug_properties = prop_dict

    for k, v in clearance_dict.items():
        if k == drug_properties['clearance_mouse']:
            drug_properties['clearance_mouse'] = v
        if k == drug_properties['clearance_human']:
            drug_properties['clearance_human'] = v
    for k, v in pampa_dict.items():
        if k == drug_properties['pampa']:
            drug_properties['pampa'] = v
        if k == drug_properties['logd']:
            drug_properties['logd'] = v
    return (drug_properties)
