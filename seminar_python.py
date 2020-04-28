# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 14:23:13 2020

@author: Marsel Rus
"""


import numpy as np
import math

# =============================================================================
# Podatki
# =============================================================================

J0 = 0.7
k0 = 1300
M0 = 1150
N = 6
T = 0.17
t0 = 1.28
tk = 4.22

eul=np.exp(1)

# Definicija praznih matrik

Ji = []
kti = []

# Dodelitev vrednosti posameznih matrik

for i in range(N):
    J = J0*eul**(1 - 0.3*i)
    Ji.append(J)
    k = k0*eul**(1 - 0.3*i)
    kti.append(k)
    
MasMx = np.zeros((N,N),float)

for i in range(N):
    MasMx[i][i] = MasMx[i][i] + J[i][i]

print(MasMx[1][1])

# =============================================================================
# Definicija funkcij
# =============================================================================

def M(t):
    
    M1 = M0*(1 + np.tanh(-2*t/t0))
    return M1

def func(t, M0=1, T=1):
    if t < (T / 2):
        M2 = (2*t / T) * M0
        return M2
    else:
        M2 = (2*t / T) * M0 - 2*M0
        return M2

# =============================================================================
# Izračuni
# =============================================================================

