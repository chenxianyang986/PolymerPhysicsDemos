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
  zero_deriv_solminIdx, _= findLocalMaximaMinima(len(yval), yval)
  if len(zero_deriv_solminIdx) != 2:
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
        '''
        if value == 0.01:
          sol.append((xval[i]*3 + value)/4)
        else:
          sol.append((xval[i]*3 + value)/4)
        '''
        sol.append((xval[i] + value)/2)
      else:
        value = xval[zero_deriv_solminIdx[0]] if len(zero_deriv_solminIdx) == 1 and zero_deriv_solminIdx > sec_deriv_min else 0.99
        '''
        if value == 0.01:
          sol.append((xval[i]*3 + value)/4)
        else:
          sol.append((xval[i]*3 + value)/4)
        '''
        sol.append((xval[i] + value)/2)
    '''
    if len(sol) == 0:
      value = int((sec_deriv_min[0] -zero_deriv_solminIdx[0]) / 2)
      sol.append(xval[sec_deriv_min[0]+value])
      sol.append(xval[sec_deriv_min[0]-value])
    '''
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
  return sol
  
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
    
  return ans

