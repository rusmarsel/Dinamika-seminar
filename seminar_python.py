# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 14:23:13 2020

@author: Marsel Rus
"""


import numpy as np
# import math
import matplotlib.pyplot as plt

# from sympy import *
# from IPython.display import Image, display
# from IPython.lib.latextools import latex_to_png

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

# Lastni vektorji
fi_last = np.zeros((N,N),float)

# Lastne vrednosti
las_vr = []

# Lastne frekvence
las_frek = []


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
# Definicija obremenitvenih funkcij
# =============================================================================


def moment1(t, M0, t0):
    '''Popis momenta M1 s funkcijo''' 
    M1 = M0*(1 + np.tanh(-2*t/t0))
    return M1


def moment2(t,T,M0):
    '''Popis momenta M2 s funkcijo'''
    if t < (T/2):
        M2 = (t / T) * 2*M0
        return M2
    
    else:
        M2 = (t / T) * 2*M0 - 2*M0
        return M2

t_fun1 = np.linspace(0,t0,num=1000)    
t_fun2 = np.linspace(0,T,num=1000)


# =============================================================================
# Grafičen izris funkcij
# =============================================================================


for i in range(len(t_fun1)):
    fun_mom1.append(moment1(t_fun1[i],M0,t0))


for i in range(len(t_fun2)):
    fun_mom2.append(moment2(t_fun2[i],T,M0))

# plt.subplot(121)
# plt.plot(t_fun1,fun_mom1)
# plt.subplot(122)
# plt.plot(t_fun2,fun_mom2)
# plt.show()

# =============================================================================
# =============================================================================
# Izračuni
# =============================================================================
# =============================================================================

Adin = np.matmul(np.linalg.inv(MasMx), Kmx)
EIG = np.linalg.eigvals(Adin)
las_vr = np.sort(EIG)

for i in range(N):
    # Izračun lastnih frekvenc sistema
    root = np.sqrt(las_vr[i])
    las_frek.append(root)


# print(las_frek)


for j in range(N):
    fi_last[0][j] = 1
    for i in range(N-1):
        if i == 0:
            fi_last[i+1][j] = (-fi_last[i][j]*(Kmx[i][i] - MasMx[i][i]*las_vr[j])) / Kmx[i][i+1]
        elif i < (N-1):
            fi_last[i+1][j] = (-Kmx[i][i-1]*fi_last[i-1][j] - (fi_last[i][j] * (Kmx[i][i] - MasMx[i][i]*las_vr[j]))) / Kmx[i][i+1]
        # print(fi_last)

i_iter = np.linspace(1,N,N)

# =============================================================================
# Izris razmerij amplitud pomikov glede na prvo komponento lastnih vektorjev
# =============================================================================

# plt.subplot(161)
# plt.plot(i_iter,fi_last[:,0])
# plt.subplot(162)
# plt.plot(i_iter,fi_last[:,1])
# plt.subplot(163)
# plt.plot(i_iter,fi_last[:,2])
# plt.subplot(164)
# plt.plot(i_iter,fi_last[:,3])
# plt.subplot(165)
# plt.plot(i_iter,fi_last[:,4])
# plt.subplot(166)
# plt.plot(i_iter,fi_last[:,5])
# plt.show()

# plt.figure()
# plt.plot(i_iter,fi_last[:,0])
# plt.grid()
# plt.title('Razmerja pomikov prvega lastnega vektorja')
# plt.xlim(1,6)
# plt.ylim(0,1.2)
# plt.xlabel('Lastna vrednost [/]')
# plt.ylabel('Razmerje [/]')
# plt.show()

# =============================================================================
# Modalne vrednosti
# =============================================================================

M_modalna = np.transpose(fi_last).dot(MasMx).dot(fi_last)
K_modalna = np.transpose(fi_last).dot(Kmx).dot(fi_last)

np.set_printoptions(suppress=True)
# print(np.around(M_modalna,decimals=2))
# print(np.around(K_modalna,decimals=2))

# =============================================================================
# Izračun integralske enačbe (odziv sistema)
# =============================================================================

# Korak trapezne interpolacije
h = .1
t_h = np.linspace(0,t0,100)

# print(moment1(t_h,M0,t0))

def trapezno(f, h):
    '''Funkcija za trapezno interpolacijo integrala'''
    trapez = (f[0] + f[-1])*h/2
    return trapez

def odziv(t, M0, t0):
    '''Funkcija za izračun odziva sistema na impulzno motnjo'''
    x_odziv = trapezno(moment1(t_h, M0, t0), h)
    return x_odziv

print(odziv(0, M0, t0))

