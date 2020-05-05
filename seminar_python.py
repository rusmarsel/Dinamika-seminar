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
N = 32      # [/]
T = 0.17    # [s]
t0 = 1.28   # [s]
tk = 4.22   # [s]

# Fourirjeva aproksimacija funkcije momenta:
n_max = 300  # [/]


# Omega

w = 2*np.pi/T 

# Časovni intervali

t_min = 0       # [s]
t_max = 3*T     # [s]
dt = 0.001      # [s]
dt_d = .001     # [s]
t_d = 3*T       # [s]

t_vec = np.arange(t_min,t_d,dt_d)
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

moment1(t_array, t0, M0, show = False)
moment2(t_array, T, M0, show = False)

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
    las_frek.append(root/(2*np.pi))


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
        
    return M_approx
    

# def mom_four(t, b_j):
#     M_four = np.ones_like(t_array) 
    
#     for j in range(len(b_j+1)):
#         if j==0:
#             M_four *= 0
#         else:
#             M_four += b_j[j]*np.sin(j*w*t)
    
#     return M_four[n_max]

def mom_four(t, a_j, b_j, show=False):
    M_four = np.ones_like(moment2) 
    
    for j in range(len(a_j+1)):
        if j==0:
            M_four *= a_j[j]/2 
        else:
            M_four += a_j[j]*np.cos(j*w*t) + b_j[j]*np.sin(j*w*t)
            
    if show:
        plt.figure()
        plt.plot(t_array, M_four)
        plt.show()
    
    return M_four
    
# Eksaktna funkcija vzbujevalnega momenta
M_exact = moment2(t_array, T, M0, show = False)  # Za prikaz -> show = True

a_j, b_j = koef_bj(M0, n_max)

primerjava(t_array, M_exact, a_j, b_j, w, show = False)  # Za prikaz -> show = True


def vectorfield(t, y, N, a_j, b_j, M_modalna, K_modalna):
    eta = []
    
    for i in range(N):
        # Pripišemo pomike sistema:
        eta.append(y[2*i])
    
    np.array(eta)
    
    # Zapis vektorja momentov (zunanje obremenitve)
    moment = np.zeros((N,1))
    moment[0,0] = mom_four(t, a_j, b_j) # Vrednost momenta, aproksimiranega s fourirjevimi vrstami, pri času t
    moment[1,0] = mom_four(t, a_j, b_j)
    
    DD = np.transpose(moment) - np.transpose(np.dot(K_modalna, np.transpose(eta)))
    M_mod_inv = np.linalg.inv(M_modalna)
    
    dy = []
    
    z = 0
    
    for i in range(2*N):
        if i % 2 == 0: # Pripis prvega ali drugega odvoda, glede na vrednost i
            dy.append(y[i+1])
        else:
            dy.append(np.dot(M_mod_inv, np.transpose(DD))[z])
            z += 1
            
    return dy 


params = (N, a_j, b_j, K_modalna, M_modalna) # vhodni parametri

init_vals = []
for i in range(2*N):
    init_vals.append(0)

rez = solve_ivp(vectorfield, [t_vec[0], t_vec[-1]], init_vals, args=params, t_eval=t_vec)


ta = rez.t
odz = 0
odz = rez.y[0,:]
# for i in range(N):
#     odz += rez.y[(i*2),:]
      
# print(odz)
plt.figure()
plt.plot(ta, (180/np.pi)*(odz), label="Prikaz odziva")
plt.xlabel('Čas [s]')
plt.ylabel('Zasuk [°]')
plt.grid()
plt.legend(loc="upper right")
plt.show()
        