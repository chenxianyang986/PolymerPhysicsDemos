import numpy as np
from phase_solve import find_guesses_for_common_tangents, find_binodal_points
from utils import F_mix_s
import matplotlib.pyplot as plt
k = 1.3807 * 10 ** -23
eps = np.finfo(np.float64).eps

def main():
    T = np.linspace(200, 400, 50)
    x = np.linspace(0.01, 0.99, 1001)
    NA = 200
    NB = 100
    chi_values = 5/T
    bpts = []
    for i in range(len(chi_values)):
        print(i)
        guesses = find_guesses_for_common_tangents(NA, NB, chi_values[i], k*T[i])
        input_guesses_y = [F_mix_s(j, NA, NB, chi_values[i], k*T[i]) for j in guesses]
        input_guesses = tuple(guesses + input_guesses_y)
        binodal_points = find_binodal_points(NA, NB, chi_values[i], k*T[i], input_guesses)
        bpts += [[k, T[i]] for k in binodal_points[:2]]
    plt_bpts = np.array(bpts)
    z = np.polyfit(np.float64(plt_bpts[:, 0]), np.float64(plt_bpts[:, 1]), 2, rcond = len(x)*eps)
    p = np.poly1d(z)
    plt.plot(plt_bpts[:, 0], plt_bpts[:, 1], '.', x, p(x), '--')
    plt.xlabel("concentrations")
    plt.ylabel("temperatures K")
    plt.show()

if __name__ == "__main__":
    main()