# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 14:23:13 2020

@author: Marsel Rus
"""

import numpy as np
import matplotlib.pyplot as plt

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
delta = 0.1 # [/] - Dušenje; namen: razširitev naloge
# k0 = (78e+9*np.pi*0.04**4)/(32*0.3) # Razširitev

# Omega

w = 2*np.pi/T 

# Časovni interval
t_min = 0      # [s]
t_max = 1.75*t0      # [s]
dt = 0.001     # [s]

t_impulz = np.arange(t_min, t_max, dt)

eul=np.exp(1)
np.set_printoptions(suppress=True)
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
    '''Popis obremenitvenega momenta M1 s funkcijo''' 
    M1 = np.zeros_like(t)    
    for i in range(len(M1)):
        if t[i]<t0:
            M1[i] = M0*np.tanh(-2*t[i]/t0)
    
    if show:        
        plt.figure()
        plt.plot(t_impulz, M1)
        plt.grid()
        plt.xlabel('Čas [s]')
        plt.ylabel('Moment [Nm]')
        plt.show()    
    return M1

# =============================================================================
# Grafičen izris funkcij
# =============================================================================
moment1(t_impulz, t0, M0, show = False)

# =============================================================================
# =============================================================================
# Izračuni
# =============================================================================
# =============================================================================

# Lastne vrednosti
Adin = np.matmul(np.linalg.inv(MasMx), Kmx)
# Kvadrati lastnih frekvenc
EIG = np.linalg.eigvals(Adin)
# Uredim po velikosti
las_vr = np.sort(EIG)

# Izračun lastnih frekvenc sistema
for i in range(N):
    root = np.sqrt(las_vr[i])
    las_frek.append(root)

# Izračun matrike lasntih vektorjev sistema 
# Normirano glede na prvo lastno vrednost vektorjev
for j in range(N):
    fi_last[0][j] = 1
    for i in range(N-1):
        if i == 0:
            fi_last[i+1][j] = (-fi_last[i][j]*(Kmx[i][i] - MasMx[i][i]*las_vr[j])) / Kmx[i][i+1]
        elif i < (N-1):
            fi_last[i+1][j] = (-Kmx[i][i-1]*fi_last[i-1][j] - (fi_last[i][j] * (Kmx[i][i] - MasMx[i][i]*las_vr[j]))) / Kmx[i][i+1]


# =============================================================================
# Modalne vrednosti
# =============================================================================
# Izračun lastne masne (M_modalna) in lastne togostne (K_modalna) matrike
M_modalna = np.transpose(fi_last).dot(MasMx).dot(fi_last)
K_modalna = np.transpose(fi_last).dot(Kmx).dot(fi_last)

# Ta vektor uporabim kot zapis faktorjev obremenitve pri prehodu v modalne koordinate
sile = np.zeros((N,1))
sile[0,0] = 1
sile[1,0] = 1
# Izračun vpliva obremenitve na posamezno prostostno stopnjo v modalnih koordinatah    
fact = np.matmul(np.transpose(fi_last), sile)
fact = np.matmul(np.linalg.inv(M_modalna),fact)

# =============================================================================
# Izračun integralske enačbe
# =============================================================================
def tangens(t):
    '''Zapis obremenitvene funkcije''' # Funkcija naredi isto kot moment1
    tan = np.zeros(len(t))
    for i in range(len(t)):
        if t[i]<t0:
            tan[i] = np.tanh(-2*(t[i]/t0))
        else:
            tan[i] = 0
    return tan


def g_t(t, PS_t):
    '''Prenosna funkcija'''
    gt = np.zeros(len(t))
    wn = las_frek[PS_t]
    for i in range(len(t)):
        gt[i] = np.sin(wn*t[i])
    return gt


def konv(t, PS):
    '''Izračun konvolucijskega integrala obremenitvene in prenosne funkcije'''
    eta = np.zeros((PS, len(t)))
    tan = tangens(t)
    for p in range(PS):
        gt = g_t(t,p)
        konv = np.convolve(tan,gt,mode='Full')
        eta[p,:] = fact[p]*M0*(1/las_frek[p])*konv[:len(t)]*dt
    return eta


def g_t_duseno(t, PS_t):
    '''Razširitev - prenosna funkcija z dušenjem'''
    gtd = np.zeros(len(t))
    wd = las_frek[PS_t]*np.sqrt(1-delta**2)
    for i in range(len(t)):
        gtd[i] = (eul**(-delta*wd*t[i]))*np.sin(wd*t[i])
    return gtd
  
    
def konv_dus(t, PS):
    '''Razširitev - konvolucija z dušenjem'''
    eta = np.zeros((PS, len(t)))
    tan = tangens(t)
    for p in range(PS):
        wd = las_frek[p]*np.sqrt(1-delta**2)
        gtd = g_t_duseno(t,p)
        konv = np.convolve(tan,gtd,mode='Full')
        eta[p,:] = fact[p]*M0*(1/wd)*konv[:len(t)]*dt
    return eta


# Zapis vektorjev in prikaz rezultatov
eta = konv(t_impulz,N)
etax = np.ones((N,len(t_impulz)))
# Vektor pomikov
x_t1 = np.zeros((N,len(t_impulz)))

# Preračun iz modalnih v fizikalne koordinate
for p in range(N):
    for i in range(len(t_impulz)):
            etax[:,i] = eta[:,i]
            x_t1[:,i] = np.matmul(fi_last, etax[:,i]) 
   
            
# Preračun iz modalnih v fizikalne koordinate - dušen sistem
etad = konv_dus(t_impulz,N)
etaxd = np.ones((N,len(t_impulz)))

# Vektor pomikov
x_t1d = np.zeros((N,len(t_impulz)))
for p in range(N):
    for i in range(len(t_impulz)):
            etaxd[:,i] = etad[:,i]
            x_t1d[:,i] = np.matmul(fi_last, etaxd[:,i]) 
            
# print(x_t1[0,int(1/dt)])


def izris(show=False):
    '''Izriše grafe pomikov sistema v časovnem intervalu'''
    if show:
        plt.figure()        
        for i in range(N):
            plt.plot(t_impulz, (180/np.pi)*x_t1d[i,:])
            # plt.plot(t_impulz, (180/np.pi)*x_t1d[i,:])  # Izris enega ali drugega načina (dušeno/nedušeno)
            plt.xlabel('Čas [s]')
            plt.ylabel('Zasuk [°]')
        plt.grid()
        plt.show()
    
    
izris(show=True)


# plt.figure()
# plt.plot(t_impulz, (180/np.pi)*x_t1[0,:], Label='Pomik $\\varphi_{1}$') #, label='Odziv $\\varphi_{'+str(i+1)+'}$')
# plt.plot(t_impulz, (180/np.pi)*x_t1[1,:], Label='Pomik $\\varphi_{2}$') 
# plt.plot(t_impulz, (180/np.pi)*x_t1[2,:], Label='Pomik $\\varphi_{3}$') 
# plt.plot(t_impulz, (180/np.pi)*x_t1[3,:], Label='Pomik $\\varphi_{4}$') 
# plt.plot(t_impulz, (180/np.pi)*x_t1[4,:], Label='Pomik $\\varphi_{5}$')
# plt.plot(t_impulz, (180/np.pi)*x_t1[5,:], Label='Pomik $\\varphi_{6}$') 
# plt.xlabel('Čas [s]')
# plt.ylabel('Zasuk [°]')
# plt.legend(loc="upper left")
# plt.grid()
# plt.show()


# print(max(x_t1[0,:])*(180/np.pi))

