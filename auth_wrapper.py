import os
from flask import Flask, request, Response
from CNPs_MFSA_Web import app 

USERNAME = os.getenv("DASH_AUTH_USERNAME", "cnps")
PASSWORD = os.getenv("DASH_AUTH_PASSWORD", "123456")

def check_auth(user, pwd):
    return user == USERNAME and pwd == PASSWORD

def authenticate():
    return Response(
        'Login required',
        401,
        {'WWW-Authenticate': 'Basic realm="Login Required"'}
    )

def requires_auth(f):
    def decorated(*args, **kwargs):
        auth = request.authorization
        if not auth or not check_auth(auth.username, auth.password):
            return authenticate()
        return f(*args, **kwargs)
    return decorated

server = app.server
server.before_request_funcs = {
    None: [requires_auth]
}
