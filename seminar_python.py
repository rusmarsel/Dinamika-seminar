# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 14:23:13 2020

@author: Marsel Rus
"""


import numpy as np
import math
import matplotlib.pyplot as plt

from sympy import *
from IPython.display import display
from IPython.display import Image, display
from IPython.lib.latextools import latex_to_png

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

# =============================================================================
# Definicija praznih matrik
# =============================================================================
# Matrika vztrajnostnih momentov:
Ji = []
# Vektor togosti sistema
kti = []
# Masna matrika
MasMx = np.zeros((N,N),float)
# Togostna matrika
Kmx = np.zeros((N,N),float)
# Vrednosti funkcije momenta M1
fun_mom1 = []
# Vrednosti funkcije momenta M2
fun_mom2 = []
#Lastni vektorji
fi_last = np.zeros((N,N),float)

# Dodelitev vrednosti posameznih matrik

for i in range(N):
    J = J0*eul**(1 - 0.3*(i+1))
    Ji.append(J)
    k = k0*eul**(1 - 0.3*(i+1))
    kti.append(k)
    
# Zapis masne matrike:

for i in range(N):
    MasMx[i][i] = MasMx[i][i] + Ji[i]
    
# Zapis togostne matrike:
    
for i in range(N):
    if i == 0:
        Kmx[i][i] = kti[i]
        Kmx[i][i+1] = -kti[i]
        
    elif i < (N-1):
        Kmx[i][i-1] = -kti[i-1]
        Kmx[i][i] = kti[i-1]+kti[i]
        Kmx[i][i+1] = -kti[i]
        
    elif i < (N):
        Kmx[i][i-1] = -kti[i-1]
        Kmx[i][i] = kti[i-1]+kti[i]

# x, y, z = symbols('x y z')
# eq = r'F(k) = \int_{-\infty}^{\infty} f(x) e^{2\pi i k} dx'
# data = latex_to_png(eq, wrap=True)
# display(Image(data=data))

# =============================================================================
# Definicija funkcij
# =============================================================================


def moment1(t,M0):
    M1 = M0*(1 + np.tanh(-2*t/t0))
    return M1



def moment2(t,T,M0):
    if t < (T/2):
        M2 = (t / T) * M0
        return M2
    
    else:
        M2 = (t / T) * M0 - M0
        return M2

# print(moment2(0.17))

t_fun1 = np.linspace(0,tk,num=1000)    
t_fun2 = np.linspace(0,T,num=1000)

for i in range(len(t_fun1)):
    fun_mom1.append(moment1(t_fun1[i],M0))

for i in range(len(t_fun2)):
    fun_mom2.append(moment2(t_fun2[i],T,M0))

# plt.subplot(121)
# plt.plot(t_fun1,fun_mom1)
# plt.subplot(122)
# plt.plot(t_fun2,fun_mom2)
# plt.show()

Adin = np.matmul(np.linalg.inv(MasMx), Kmx)
# print(np.linalg.inv(MasMx))
# print((np.around(Adin)))

EIG = np.linalg.eigvals(Adin)
# print(EIG[3])

# print(np.around(Kmx,decimals=2)) 
# print(np.around(Adin,decimals=2))

# =============================================================================
# Izračuni
# =============================================================================

# for j in range(N):
#     print(EIG[j])

for j in range(N):
    fi_last[0][j] = 1
    for i in range(N-1):
        if i == 0:
            fi_last[i+1][j] = (fi_last[i][j]*(Kmx[i][i]-EIG[j]*MasMx[i][i]))/Kmx[i][i+1]
        elif i < (N-1):
            fi_last[i+1][j] = (fi_last[i][j]*(Kmx[i][i]-EIG[j]*MasMx[i][i])-Kmx[i][i-1])/Kmx[i][i+1]
        elif i == (N-1):
            fi_last[i+1][j] = (fi_last[i][j]*k[i+1][i])/(Kmx[i+1][i+1]-EIG[j]*Masmx[i+1][i+1])
        # print(fi_last)
        
print(np.around(fi_last,decimals=2))