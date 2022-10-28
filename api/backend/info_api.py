from flask import jsonify


def info_text():
    """Returns info help text as vales in a dictionary, with keys as the
    info button labels

    :return: JSON containing help info text
    :rtype: json dict
    """
    info_dict = {}
    with open("./src/info_text/info_box_build.txt", "r", encoding='UTF-8') as f:
        q = f.readlines()
        lines = [i.rstrip("\n") for i in q]
        info_dict["build"] = lines
    with open("./src/info_text/info_box_assay.txt", "r", encoding='UTF-8') as f:
        q = f.readlines()
        lines = [i.rstrip("\n") for i in q]
        info_dict["assay"] = lines
    with open("./src/info_text/info_box_analysis.txt", "r", encoding='UTF-8') as f:
        q = f.readlines()
        lines = [i.rstrip("\n") for i in q]
        info_dict["analysis"] = lines

    return jsonify({"info_dict": info_dict})
