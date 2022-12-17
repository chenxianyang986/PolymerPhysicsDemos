import numpy as np
from utils import F_mix_s
import matplotlib.pyplot as plt
import argparse
from phase_solve import find_guesses_for_common_tangents, find_binodal_points
k = 1.3807 * 10 ** -23

def run():
    #one test data NA = 200 NB = 100 chi = 5/T T = (250, 350)
    #other test data NA = 6000 NB = 5000 chi = -7.74/T + 0.018, T = (438, 450)
    #other test data NA=2000 NB=1500 chi=-9.2*10^-4+0.722/T, T = (240, 390)
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', help="please provide the test number (1, 2, 3)")
    parser.add_argument('--temp', help="please provide temperature")
    args = parser.parse_args()
    match int(args.test):
        case 2:
            test_2(float(args.temp))
        case 3:
            test_3(float(args.temp))
        case 4:
            test_4(float(args.temp))
        case 5:
            test_5(float(args.temp))
        case _:
            test_1(float(args.temp))

def test_1(T):
    if T < 200 or T > 1000:
        print("provide a temperature between 200K and 1000K")
        return
    x = np.linspace(0.01, 0.99, 1001)
    y = []
    NA = 200
    NB = 100
    chi = 5/T
    for i in x:
        y.append(F_mix_s(i, NA, NB, chi, k*T))
    plt.plot(x, y, label="Free energy curve")
    guesses = find_guesses_for_common_tangents(NA, NB, chi, k*T)
    input_guesses_y = [F_mix_s(i, NA, NB, chi, k*T) for i in guesses]
    input_guesses = tuple(guesses + input_guesses_y)
    binodal_points = find_binodal_points(NA, NB, chi, k*T, input_guesses)
    ''' 
    uncomment this for plot the guesses
    guesses_yvalue = []
    for i in guesses:
        guesses_yvalue.append(F_mix(i, NA, NB, chi, k*T))
    plt.plot(guesses, guesses_yvalue, 'o', label="transitional points")
    '''
    plt.plot([binodal_points[0], binodal_points[1]], [F_mix_s(binodal_points[0], NA, NB, chi, k*T), F_mix_s(binodal_points[1], NA, NB, chi, k*T)], 'o-', label="common tangent line")
    plt.legend()
    plt.show()

def test_2(T):
    if T < 400 or T > 460:
        print("provide a temperature between 400K and 460K")
        return
    x = np.linspace(0.01, 0.99, 1001)
    y = []
    NA = 6000
    NB = 5000
    chi = -7.74/T + 0.018
    for i in x:
        y.append(F_mix_s(i, NA, NB, chi, k*T))
    plt.plot(x, y, label="Free energy curve")
    guesses = find_guesses_for_common_tangents(NA, NB, chi, k*T)
    input_guesses_y = [F_mix_s(i, NA, NB, chi, k*T) for i in guesses]
    input_guesses = tuple(guesses + input_guesses_y)
    binodal_points = find_binodal_points(NA, NB, chi, k*T, input_guesses)
    plt.plot([binodal_points[0], binodal_points[1]], [F_mix_s(binodal_points[0], NA, NB, chi, k*T), F_mix_s(binodal_points[1], NA, NB, chi, k*T)], 'o-', label="common tangent line")
    plt.legend()
    plt.show()

def test_3(T):
    if T < 240 or T > 700:
        print("provide a temperature between 240K and 700K")
        return
    x = np.linspace(0.01, 0.99, 1001)
    y = []
    NA = 2000
    NB = 1500
    chi = -9.2*10**-4+0.722/T
    for i in x:
        y.append(F_mix_s(i, NA, NB, chi, k*T))
    plt.plot(x, y, label="Free energy curve")
    guesses = find_guesses_for_common_tangents(NA, NB, chi, k*T)
    input_guesses_y = [F_mix_s(i, NA, NB, chi, k*T) for i in guesses]
    input_guesses = tuple(guesses + input_guesses_y)
    binodal_points = find_binodal_points(NA, NB, chi, k*T, input_guesses)
    plt.plot([binodal_points[0], binodal_points[1]], [F_mix_s(binodal_points[0], NA, NB, chi, k*T), F_mix_s(binodal_points[1], NA, NB, chi, k*T)], 'o-', label="common tangent line")
    plt.legend()
    plt.show()

def test_4(T):
    if T < 200 or T > 1000:
        print("provide a temperature between 250K and 1000K")
        return
    x = np.linspace(0.01, 0.99, 1001)
    y = []
    NA = 100
    NB = 200
    chi = 5/T
    for i in x:
        y.append(F_mix_s(i, NA, NB, chi, k*T))
    plt.plot(x, y, label="Free energy curve")
    guesses = find_guesses_for_common_tangents(NA, NB, chi, k*T)
    input_guesses_y = [F_mix_s(i, NA, NB, chi, k*T) for i in guesses]
    input_guesses = tuple(guesses + input_guesses_y)
    binodal_points = find_binodal_points(NA, NB, chi, k*T, input_guesses)
    plt.plot([binodal_points[0], binodal_points[1]], [F_mix_s(binodal_points[0], NA, NB, chi, k*T), F_mix_s(binodal_points[1], NA, NB, chi, k*T)], 'o-', label="common tangent line")
    plt.legend()
    plt.show()

def test_5(T):
    
    if T < 200 or T > 600:
        print("provide a temperature between 250K and 600K")
        return
    
    x = np.linspace(0.01, 0.99, 1001)
    y = []
    NA = 200
    NB = 200
    chi = 5/T
    for i in x:
        y.append(F_mix_s(i, NA, NB, chi, k*T))
    plt.plot(x, y, label="Free energy curve")
    
    guesses = find_guesses_for_common_tangents(NA, NB, chi, k*T)
    input_guesses_y = [F_mix_s(i, NA, NB, chi, k*T) for i in guesses]
    input_guesses = tuple(guesses + input_guesses_y)
    binodal_points = find_binodal_points(NA, NB, chi, k*T, input_guesses)
    plt.plot([binodal_points[0], binodal_points[1]], [F_mix_s(binodal_points[0], NA, NB, chi, k*T), F_mix_s(binodal_points[1], NA, NB, chi, k*T)], 'o-', label="common tangent line")
    
    plt.legend()
    plt.show()

if __name__ == "__main__":
    run()
