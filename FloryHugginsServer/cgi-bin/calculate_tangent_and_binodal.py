#!/usr/bin/python3

import sys
import json
from http.server import HTTPServer, SimpleHTTPRequestHandler
from sympy.solvers import nsolve
from sympy import Symbol, log

hostName = "localhost"
serverPort = 8000

class MyServer(SimpleHTTPRequestHandler):
    def do_POST(self):
        self.send_response(200)
        self.send_header("Content-type", "text/html")
        self.end_headers()
        print("get here!")
    def do_GET(self):
        return SimpleHTTPRequestHandler.do_GET(self)
def solve_for_tangent_and_binodal(NA, NB, chi, kT):
    x1 = Symbol('x1')
    x2 = Symbol('x2')
    y1 = Symbol('y1')
    y2 = Symbol('y2')

    def d_d_F_mix_s(phi,NA,NB,chi,kT):
        return -2*chi + 1./(NA*phi) + 1./(NB - NB*phi)

    def d_F_mix_s(phi,NA,NB,chi,kT):
        return kT*(chi*(1.-2*phi) + 1/NA*log(phi)+ 1/NA -1/NB -1/NB*log(1-phi))

    def F_mix_s(phi,NA,NB,chi,kT):
        return kT*(chi*phi*(1.-phi) + phi/NA*log(phi) + (1.-phi)/NB*log(1-phi))

    try:
        sol = nsolve((
        y1-F_mix_s(x1,NA,NB,chi,kT),
        y2-F_mix_s(x2,NA,NB,chi,kT),
        d_F_mix_s(x1,NA,NB,chi,kT)-d_F_mix_s(x2,NA,NB,chi,kT),
        d_F_mix_s(x1,NA,NB,chi,kT)-(y2-y1)/(x2-x1)),
        (x1,x2,y1,y2),(0.01,0.99,-0.01,-0.01))
        tangent=[float(sol[0]),float(sol[1]),float(sol[2]),float(sol[3])]
    except:
        tangent=["nan", "nan", "nan", "nan"]

    try:
        sol = nsolve((
        d_d_F_mix_s(x1,NA,NB,chi,kT),
        d_d_F_mix_s(x2,NA,NB,chi,kT),
        y1-F_mix_s(x1,NA,NB,chi,kT),
        y2-F_mix_s(x2,NA,NB,chi,kT)),
        (x1,x2,y1,y2),(1./4.,3./4.,-0.01,-0.01))
        binodal=[float(sol[0]),float(sol[1]),float(sol[2]),float(sol[3])]
    except:
        binodal=["nan","nan", "nan", "nan"]

    return [tangent, binodal]

if __name__ == "__main__":
    webServer = HTTPServer((hostName, serverPort), MyServer)
    print("Server started http://%s:%s" % (hostName, serverPort))
    try:
        webServer.serve_forever()
    except KeyboardInterrupt:
        pass
    '''
    parameters = json.load(sys.stdin)
    result = solve_for_tangent_and_binodal(parameters["na"], parameters["nb"], parameters["chi"], 1.0)
    #with open("tangent_and_binodals", 'w') as f:
    #    print(result, file=f)
    print("Content-Type: application/json")
    print("Access-Control-Allow-Origin: *")
    print()
    print(json.dumps(result))
    '''
