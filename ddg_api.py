import flask
from flask import Flask
import base64


app = Flask(__name__)

@app.route("/")
def hello_world():
    return "<p>Hello, World!</p>"


@app.route("/r-group-<string:r_group_id>")
def rgroup_img(r_group_id):
    img_filename = f'./images/{r_group_id}.png'
    return flask.send_file(img_filename, mimetype='image/png')