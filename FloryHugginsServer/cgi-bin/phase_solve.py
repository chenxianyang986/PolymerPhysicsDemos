import numpy as np
from sympy import Symbol, nsolve
from utils import findLocalMaximaMinima, F_mix_s, d_F_mix_s, d_d_F_mix_s
import matplotlib.pyplot as plt

def find_guesses_for_common_tangents(NA, NB, chi, kT, range=(0.01, 0.99), number_of_points=1001):
  xval = np.linspace(range[0], range[1], number_of_points)
  yval = []
  sol = []
  for i in xval:
    yval.append(F_mix_s(i, NA, NB, chi, kT))
  zero_deriv_solminIdx, zero_deriv_max= findLocalMaximaMinima(len(yval), yval)
  min_of_zero_min = max([F_mix_s(xval[i], NA, NB, chi, kT) for i in zero_deriv_solminIdx])
  #print(F_mix_s(xval[zero_deriv_max[0]], NA, NB, chi, kT),F_mix_s(xval[zero_deriv_solminIdx[0]], NA, NB, chi, kT))
  yval_max_zero_deriv = F_mix_s(xval[zero_deriv_max[0]], NA, NB, chi, kT) if len(zero_deriv_max) == 1 else min_of_zero_min
  #if len(zero_deriv_max) == 1:
    #print(abs((min_of_zero_min - yval_max_zero_deriv)/min_of_zero_min))
  if len(zero_deriv_solminIdx) != 2 or abs((min_of_zero_min - yval_max_zero_deriv)/min_of_zero_min) < 0.8:
    yval.clear()
    for i in xval:
      yval.append(d_F_mix_s(i, NA, NB, chi, kT))
    solminIdx, solmaxIdx = findLocalMaximaMinima(len(yval), yval)
    yval.clear()
    for i in xval:
      yval.append(d_d_F_mix_s(i, NA, NB, chi, kT))
    sec_deriv_min, _ = findLocalMaximaMinima(len(yval), yval)
    for i in solminIdx + solmaxIdx:
      if xval[i] < xval[sec_deriv_min]:
        value = xval[zero_deriv_solminIdx[0]] if len(zero_deriv_solminIdx) == 1 and zero_deriv_solminIdx < sec_deriv_min else 0.01
        #xval[abs(sec_deriv_min[0] * 2 - zero_deriv_solminIdx[0])]
        sol.append((xval[i]*3 + value)/4)
      else:
        value = xval[zero_deriv_solminIdx[0]] if len(zero_deriv_solminIdx) == 1 and zero_deriv_solminIdx > sec_deriv_min else 0.99
        #xval[len(xval) - abs((sec_deriv_min[0] * 2 - zero_deriv_solminIdx[0]) - len(xval))]
        sol.append((xval[i]*3 + value)/4)
  else:
    for i in zero_deriv_solminIdx:
      sol.append(xval[i])
  return sol

def find_spinodal_points(NA, NB, chi, kT, range=(0.01, 0.99), number_of_points=1001):
  xval = np.linspace(range[0], range[1], number_of_points)
  yval = []
  sol = []
  for i in xval:
    yval.append(d_F_mix_s(i, NA, NB, chi, kT))
  sec_deriv_min, sec_deriv_max = findLocalMaximaMinima(len(yval), yval)
  for i in sec_deriv_min+sec_deriv_max:
    sol.append(xval[i])
  if len(sol) != 2:
    return []
  return [float(sol[0]), float(sol[1])]
  
def find_binodal_points(NA, NB, chi, kT, guesses):
  # guesses is in form of (x1, x2, y1, y2)
  x1 = Symbol('x1')
  x2 = Symbol('x2')
  y1 = Symbol('y1')
  y2 = Symbol('y2')
  ans = np.zeros((1, 4))
  
  if len(guesses) == 4:
    try:
      ans = nsolve((
          y1-F_mix_s(x1,NA,NB,chi,kT),
          y2-F_mix_s(x2,NA,NB,chi,kT), 
          d_F_mix_s(x1,NA,NB,chi,kT)-d_F_mix_s(x2,NA,NB,chi,kT),
          d_F_mix_s((x1 + x2)/2,NA,NB,chi,kT)-(y2-y1)/(x2-x1)),
          (x1,x2,y1,y2), guesses)
    except:
      ans = guesses
  else:
    print("Warning! no binodal points at this temperature")
    return []
    
  return [float(ans[0]), float(ans[1]), float(ans[2]), float(ans[3])]