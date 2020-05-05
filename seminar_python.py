# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 14:23:13 2020

@author: Marsel Rus
"""

# import seminar_knjiznica as sk
import numpy as np
import matplotlib.pyplot as plt

# from sympy import *
from scipy.integrate import solve_ivp

# =============================================================================
# Podatki
# =============================================================================

J0 = 0.7    # [kgm^2]
k0 = 1300   # [Nm/rad]
M0 = 1150   # [Nm]
N = 6       # [/]
T = 0.17    # [s]
t0 = 1.28   # [s]
tk = 4.22   # [s]

# Fourirjeva aproksimacija funkcije momenta:
n_max = 200  # [/]
# t = 10;     # [s]

# Omega

w = 2*np.pi/T 

# Časovni intervali

t_min = 0       # [s]
t_max = 2*T     # [s]
dt = 0.001     # [s]

t_array = np.arange(t_min, t_max, dt)
t_impulz = np.arange(t_min, t0, dt)

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


# =============================================================================
# Definicija obremenitvenih funkcij
# =============================================================================


def moment1(t, t0, M0, show = False):
    '''Popis momenta M1 s funkcijo''' 
    M1 = np.ones_like(t_impulz)
    
    for (t, i) in zip(t_impulz, range(len(M1))):
        M1[i] = M0*(1 + np.tanh(-2*t/t0))
    
    if show:        
        plt.figure()
        plt.plot(t_impulz, M1)
        plt.grid()
        plt.xlabel('Čas [s]')
        plt.ylabel('Moment [Nm]')
        plt.show()
    
    return M1


def moment2(t,T,M0, show = False):
    '''Popis momenta M2 s funkcijo'''

    M = np.ones_like(t_array)

    c1 = 0
    c2 = 0 
    for (t, i) in zip(t_array, range(len(M))):
        if (t % T < T/2):
            M[i] = ((t-c1) / T) * 2*M0
            c2 = t
        
        
        elif (t % T > T/2 and t % T < T):
            M[i] = 2*M0*((t-c2) / T - 1/2) 
            c1 = t
        
    
    if show:        
        plt.figure()
        plt.plot(t_array, M)
        plt.grid()
        plt.xlabel('Čas [s]')
        plt.ylabel('Moment [Nm]')
        plt.show()
    
    return M
    
    
# =============================================================================
# Grafičen izris funkcij
# =============================================================================

moment1(t, t0, M0, show = False)
moment2(t, T, M0, show = False)

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


# =============================================================================
# Izris razmerij amplitud pomikov glede na prvo komponento lastnih vektorjev
# =============================================================================


def las_vek_izr(fi_last, N,show=False):
    '''Izriše grafe razmerij lastnih vektorjev'''
    
    i_iter = np.linspace(1,N,N)
    if show:
        for i in range(len(i_iter)):
            
            plt.figure()
            plt.plot(i_iter,fi_last[:,i])
            plt.grid()
            plt.title('Razmerja pomikov ' + str(i+1) + '. lastnega vektorja')
            plt.xlim(1,N)
            plt.ylim(min(fi_last[:,i]),max(fi_last[:,i]))
            plt.xlabel('Lastna vrednost [/]')
            plt.ylabel('Razmerje [/]')
            plt.show()
        
las_vek_izr(fi_last,N,show=False)

# =============================================================================
# Modalne vrednosti
# =============================================================================

M_modalna = np.transpose(fi_last).dot(MasMx).dot(fi_last)
K_modalna = np.transpose(fi_last).dot(Kmx).dot(fi_last)

np.set_printoptions(suppress=True)
# print(np.around(M_modalna,decimals=2))
# print(np.around(K_modalna,decimals=2))


def M_vzb_modalna(N):
   '''Izračuna matriko modalnih vrednosti momenta'''
   M_vzb = []
   for i in range(N):
       if i < 2:
           M_vzb.append([1])
       else:
           M_vzb.append([0])
   return np.matmul(np.transpose(fi_last), M_vzb)
    
# print(M_vzb_en(N))

# M_vzb_mod = (np.matmul(np.transpose(fi_last), M_vzb_en(N)))
# print(M_vzb_modalna(N))

# =============================================================================
# Izračun integralske enačbe (odziv sistema)
# =============================================================================

def koef_bj(M0, n_max):
    
    a_j=np.zeros(n_max)
    b_j=np.zeros(n_max)
    
    for j in range(n_max):
        
        if j==0: #a_0
        # Srednja vrednost funkcije momenta je 0 Nm
            a_j[j] = 0
        
        else:
            # Vrednosti a_j in b_j sta bili pridobljeni iz rešitve integrala v simbolični obliki v programu ___.py. 
            # Ker je vrednost a_j zelo majhna (skoraj nična), sem ji pripisal vrednost nič (razlog je lihost funkcije M2).
            # a_j[j] = 11.7647058823529*Piecewise((-5.6843418860808e-14, Eq(j, 0)), (31.1147913744655*sin(3.14159265358979*j)/j + 31.1147913744655*sin(9.42477796076938*j)/j - 9.90414570103852*cos(3.14159265358979*j)/j**2 + 9.90414570103852*cos(9.42477796076938*j)/j**2, True))
            a_j[j] = 0
            b_j[j] = 11.7647058823529*((-31.1147913744655*np.cos(3.14159265358979*j)/j + 9.90414570103852*np.sin(3.14159265358979*j)/j**2)+(-31.1147913744655*np.cos(3.14159265358979*j)/j - 9.90414570103852*np.sin(3.14159265358979*j)/j**2 + 9.90414570103852*np.sin(6.28318530717959*j)))
            
    return a_j, b_j

    
def primerjava(t, moment2, a_j, b_j, w, show = False):

    M_approx = np.ones_like(moment2) 
    
    for j in range(len(a_j+1)):
        if j==0:
            M_approx *= a_j[j]/2 
        else:
            M_approx += a_j[j]*np.cos(j*w*t) + b_j[j]*np.sin(j*w*t)

    if show:   
        plt.figure(dpi=70)
        plt.plot(t, M_approx, label="Aproksimacija momenta")
        plt.plot(t, moment2, label="Funkcija momenta")
        plt.xlabel('Čas [s]')
        plt.ylabel('Moment [Nm]')
        plt.grid()
        plt.legend(loc="upper right")
        plt.show()
    
    
def mom_four(t, b_j):
    M_four = np.ones_like(t_array) 
    
    for j in range(len(b_j+1)):
        if j==0:
            M_four *= 0
        else:
            M_four += b_j[j]*np.sin(j*w*t)
    
    return M_four
    
# Eksaktna funkcija vzbujevalnega momenta
M_exact = moment2(t, T, M0, show = False)  # Za prikaz -> show = True

a_j, b_j = koef_bj(M0, n_max)

primerjava(t_array, M_exact, a_j, b_j, w, show = False)  # Za prikaz -> show = True


# plt.figure()
# plt.plot(t_array, mom_four(t_array, b_j))

def beta(N, n_max):
        
    Beta = np.ones(N)
    
    for i in range(N):
        for j in range(n_max):
            
            Beta[i] += 1/(1-((j*w)/las_frek[i])**2)
        
    return Beta


def X_3(N, b_j, n_max):
    
    X = np.ones(N) 
    
    for i in range(N):
        for j in range(n_max):            
            X[i] += b_j[j] / (MasMx[i,i]*(las_frek[i])**2)
    
    return X


def odziv(t, n_max, N, X_3, beta):
    
    Beta = beta(N,n_max)
    X3 = X_3(N, b_j,n_max)
    
    odz = np.ones_like(t_array) 
    
    for i in range(N):
        for j in range(n_max):
            odz += np.sin(j*w*t)*Beta[i]*X3[i]*np.sin(j*w*t)
    
        plt.figure()
        plt.plot(t, odz)
        plt.show()
    return odz
        
odziv(t_array, n_max, N, X_3, beta)

# Beta = beta(N,n_max)
# X3 = X_3(N, b_j,n_max)
# print(Beta[0]*X3[0]*np.sin(1*w*t_array[1]))

# X = X_3(N, b_j, n_max)
# beta_1 = beta(N, n_max)

# print(X[0]*beta_1[0])

# def Odziv(N, n_max, moment2, beta):
    
#     Od = np.ones((len(t_array), N))
    
#     for i in range(N):
#         for j in range(n_max):
            
#             Od[i] += beta 
        
#     return Od



# def g_t(t, M_vzb_modalna, mom_four):
#     M_x = mom_four(t, b_j)
#     mod_m = M_vzb_modalna(N)
#     return mod_m * M_x[n_max]

# print(g_t(0.1,M_vzb_modalna, mom_four))


# def eta(t, N):
    
#     phi_i = np.ones(N)
    
#     for i in range(N):
#         G = g_t(t, M_vzb_modalna, mom_four)[i,0]
#         B = Beta[i]
#         OM = las_vr[i]
#         masa = M_modalna[i,i]
#         # print(B)
        
#         phi_i[i] = (180/np.pi)*np.sin(w*t)*(((G/masa) * B)/OM)

#     return phi_i
    
# Eta = eta(t_array, N)
# print(t)
# print(fi_last, "\n")
# print(Eta, "\n")
# print(np.matmul(fi_last, Eta))
       
# plt.figure()
# plt.plot(t_array, phi_i)
            



# b_j[j]*Beta[i]


# plt.figure()
# plt.plot(t_array, mom_four(t_array, b_j, w))


    
# def g_t(fi_last, mom_four):
#     np.transpose(fi_last)