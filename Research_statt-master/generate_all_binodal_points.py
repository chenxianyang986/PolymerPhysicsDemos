import numpy as np
from phase_solve import find_guesses_for_common_tangents, find_binodal_points, find_spinodal_points
from utils import F_mix_s
import matplotlib.pyplot as plt
k = 1.3807 * 10 ** -23
eps = np.finfo(np.float64).eps

def main():
    
    T = np.append(np.linspace(200, 340, 20), np.linspace(340, 345, 30))
    x = np.linspace(0.01, 0.99, 1001)
    NA = 100
    NB = 200
    chi_values = 5/T
    '''
    T = np.append(np.linspace(438.9, 439.4, 40), np.linspace(439.4, 450, 30))
    x = np.linspace(0.01, 0.99, 1001)
    NA = 6000
    NB = 5000
    chi_values = -7.74/T + 0.018
    '''
    bpts = []
    spts = []
    for i in range(len(chi_values)):
        print(i)
        guesses = find_guesses_for_common_tangents(NA, NB, chi_values[i], k*T[i])
        input_guesses_y = [F_mix_s(j, NA, NB, chi_values[i], k*T[i]) for j in guesses]
        input_guesses = tuple(guesses + input_guesses_y)
        binodal_points = find_binodal_points(NA, NB, chi_values[i], k*T[i], input_guesses)
        spinodal_points = find_spinodal_points(NA, NB, chi_values[i], k*T[i])
        spts += [[k, T[i]] for k in spinodal_points[:2]]
        bpts += [[k, T[i]] for k in binodal_points[:2]]
    
    plt_bpts = np.array(bpts, dtype=object)
    plt_spts = np.array(spts, dtype=object)
    
    a = np.array(plt_bpts[:, 0], dtype=np.float64)
    b = np.array(plt_bpts[:, 1], dtype=np.float64)
    
    z = np.polyfit(a, b, 10, rcond = len(x)*eps)
    p = np.poly1d(z)
    
    plt.plot(a, b, '.', x, p(x), '--', label="binodal points")
    plt.plot(plt_spts[:, 0], plt_spts[:, 1], '.', label = "spinodal points" )
    plt.xlabel("concentrations")
    plt.ylabel("temperatures K")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()