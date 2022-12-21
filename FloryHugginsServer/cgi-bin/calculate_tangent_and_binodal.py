#!/usr/bin/python3

import sys
import json
from http.server import HTTPServer, SimpleHTTPRequestHandler
from sympy.solvers import nsolve
from phase_solve import find_guesses_for_common_tangents, find_binodal_points, find_spinodal_points
from utils import F_mix_s

hostName = "localhost"
serverPort = 8000

k = 1.38 * 10 ** -23

class MyServer(SimpleHTTPRequestHandler):
    def do_POST(self):
        self.send_response(200)
        self.send_header("Content-type", "application/json")
        self.end_headers()
        length = int(self.headers['Content-Length'])
        params = json.loads(self.rfile.read(length))
        result = solve_for_tangent_and_spinodal(params["na"], params["nb"], params["chi"], params["T"] * k)
        f = open("tangenet and binodal", 'w+')
        json.dump(result, f)
        f.close()
        f= open("tangenet and binodal", 'r')
        self.wfile.write(f.read().encode('utf-8'))
        return
    def do_GET(self):
        return SimpleHTTPRequestHandler.do_GET(self)

def solve_for_tangent_and_spinodal(NA, NB, chi, kT):
    guesses = find_guesses_for_common_tangents(NA, NB, chi, kT)
    input_guesses_y = [F_mix_s(i, NA, NB, chi, kT) for i in guesses]
    input_guesses = tuple(guesses + input_guesses_y)
    binodal = find_binodal_points(NA, NB, chi, kT, input_guesses)
    spinodal_xvalues = find_spinodal_points(NA, NB, chi, kT)
    spinodal_yvalues = [float(F_mix_s(i, NA, NB, chi, kT)) for i in spinodal_xvalues]
    spinodal = spinodal_xvalues + spinodal_yvalues
    return [binodal, spinodal]

if __name__ == "__main__":
    webServer = HTTPServer((hostName, serverPort), MyServer)
    print("Server started http://%s:%s" % (hostName, serverPort))
    try:
        webServer.serve_forever()
    except KeyboardInterrupt:
        pass
