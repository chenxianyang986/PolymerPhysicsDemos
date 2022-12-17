import numpy as np
from sympy import log
def findLocalMaximaMinima(n, arr):
 
    # Empty lists to store points of
    # local maxima and minima
    mx = []
    mn = []
 
    # Iterating over all points to check
    # local maxima and local minima
    for i in range(1, n-1):
 
        # Condition for local minima
        if(arr[i-1] > arr[i] < arr[i + 1]):
            mn.append(i)
 
        # Condition for local maxima
        elif(arr[i-1] < arr[i] > arr[i + 1]):
            mx.append(i)
    return mn, mx

def F_mix(phi,NA,NB,chi,kT):
    return kT*(chi*phi*(1.-phi) + phi/NA*np.log(phi) + (1.-phi)/NB*np.log(1-phi))

def d_F_mix(phi, NA, NB, chi, kT):
    return kT*(chi*(1.-2*phi) + 1/NA*np.log(phi)+ 1/NA -1/NB -1/NB*np.log(1-phi))

def F_mix_s(phi,NA,NB,chi,kT):
    return kT*(chi*phi*(1.-phi) + phi/NA*log(phi) + (1.-phi)/NB*log(1-phi))

def d_F_mix_s(phi,NA,NB,chi,kT):
    return kT*(chi*(1.-2*phi) + 1/NA*log(phi)+ 1/NA -1/NB -1/NB*log(1-phi))

def d_d_F_mix_s(phi, NA, NB, chi, kT):
    return kT*(1/(NA*phi)+1/(NB*(1-phi))-2*chi)