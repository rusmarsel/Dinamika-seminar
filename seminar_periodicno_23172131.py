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

# Fourirjeva aproksimacija funkcije momenta:
n_max = 50  # [/]

# Omega

w = 2*np.pi/T 

# Časovni intervali

t_min = 0.       # [s]
t_max = T     # [s]
dt = 0.0005     # [s]

t_array = np.arange(t_min, t_max, dt)
# t_array = np.linspace(0,T,1000)
t_impulz = np.arange(t_min, t0, dt)

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
    '''Popis momenta M1 s funkcijo''' 
    M1 = np.ones_like(t_impulz)    
    for (t, i) in zip(t_impulz, range(len(M1))):
        M1[i] = M0*np.tanh(-2*t/t0)
        
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
moment1(t_array, t0, M0, show = False) # Kot je na sliki, ne uporabljam
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
def las_vek_izr(fi_last, N,show=True):
    '''Izriše grafe razmerij lastnih vektorjev'''    
    i_iter = np.linspace(1,N,N)
    if show:  
        plt.figure()
        for i in range(len(i_iter)): 
            plt.plot(i_iter,fi_last[:,i], label=str(i+1)+'. lastni vektor')
            plt.title('Razmerja pomikov ' + str(i+1) + '. lastnega vektorja')
            # plt.title('Razmerja amplitud pomikov glede na prvo lastno frekvenco')
            plt.xlim(1,N)
            plt.ylim(min(fi_last[:,i]),max(fi_last[:,i]))
        plt.xlabel('Lastna vrednost [/]')
        plt.ylabel('Razmerje [/]')
        plt.legend()
        plt.grid()
        plt.show()   
        return ()
las_vek_izr(fi_last,N,show=False)

# =============================================================================
# Modalne vrednosti
# =============================================================================
# Izračun lastne masne (M_modalna) in lastne togostne (K_modalna) matrike
M_modalna = np.transpose(fi_last).dot(MasMx).dot(fi_last)
K_modalna = np.transpose(fi_last).dot(Kmx).dot(fi_last)

# =============================================================================
# Izračun integralske enačbe
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
            # Ker je vrednost a_j zelo majhna (numerična ničla), sem ji pripisal vrednost nič (razlog je lihost funkcije M2).
            # a_j[j] = 11.7647058823529*Piecewise((-5.6843418860808e-14, Eq(j, 0)), (31.1147913744655*sin(3.14159265358979*j)/j + 31.1147913744655*sin(9.42477796076938*j)/j - 9.90414570103852*cos(3.14159265358979*j)/j**2 + 9.90414570103852*cos(9.42477796076938*j)/j**2, True))
            a_j[j] = 0
            b_j[j] = 2*M0*(-np.cos(j*np.pi)/(j*np.pi))            
    return a_j, b_j

    
def primerjava(t, moment2, a_j, b_j, w, show = False):
    '''Primerja eksaktno vredost z aproksimirano'''
    M_approx = np.ones_like(t_array)     
    for j in range(len(a_j+1)):
        if j==0:
            M_approx *= a_j[j]/2 
        else:
            M_approx += a_j[j]*np.cos(j*w*t) + b_j[j]*np.sin(j*w*t)
            
    if show:   
        plt.figure()
        plt.plot(t, M_approx, label="Aproksimacija momenta")
        plt.plot(t, moment2, label="Funkcija momenta")
        plt.xlabel('Čas [s]')
        plt.ylabel('Moment [Nm]')
        plt.grid()
        plt.legend(loc="upper right")
        plt.show()


# Eksaktna funkcija vzbujevalnega momenta (za preverjanje, koda iz vaj)
M_exact = moment2(t_array, T, M0, show = False)  # Za prikaz -> show = True
a_j, b_j = koef_bj(M0, n_max)
primerjava(t_array, M_exact, a_j, b_j, w, show = False)  # Za prikaz -> show = True


# =============================================================================
# Modalna obremenitev
# =============================================================================
# Ta vektor uporabim kot zapis faktorjev pri prehodu v modalne koordinate
sile = np.zeros((N,1))
sile[0,0] = 1
sile[1,0] = 1    
fact1 = np.matmul(np.transpose(fi_last), sile)
fact = np.matmul(np.linalg.inv(M_modalna),fact1)


def Eta_vec(t, PS):
    '''Izračun vektorja eta'''
    eta = np.ones(len(t))
    for i in range(len(t)):
        for j in range(n_max):
            if j==0:
                eta[i] = 0
            if ((j*w) < las_frek[PS]):       # primerja, dokler je j*w < w0i
                eta[i] += (fact[PS]/(las_vr[PS]))*b_j[j]*np.sin(j*w*t[i])/np.sqrt((1-((j*w)/las_frek[PS])**2)**2)         
            else:
                eta[i] += (fact[PS]/(las_vr[PS]))*b_j[j]*np.sin(j*w*t[i]-np.pi)/np.sqrt((1-((j*w)/las_frek[PS])**2)**2)
    return eta


Vec_Eta = np.ones((N,len(t_array)))
etax = np.ones((N,len(t_array)))
# Vektor pomikov
x_t = np.ones((N,len(t_array)))

# Preračun iz modalnih v fizikalne koordinate
for i in range(N):
    Vec_Eta[i,:] = Eta_vec(t_array, i)
    
for j in range(len(t_array)):
    etax[:,j] = Vec_Eta[:,j]
    x_t[:,j] = np.matmul(fi_last, etax[:,j])    


def izris(show=False):
    '''Izriše grafe pomikov sistema v časovnem intervalu'''
    if show:
        plt.figure()
        for i in range(N):        
            plt.plot(t_array, (180/np.pi)*x_t[i,:]) #, label='Odziv $\\varphi_{'+str(i+1)+'}$')
            plt.xlabel('Čas [s]')
            plt.ylabel('Zasuk [°]')
        plt.legend(loc="upper right")
        plt.grid()
        plt.show()

# print((180/np.pi)*min(x_t[0,:]))
# print((180/np.pi)*max(x_t[0,:]))
izris(show=True)

# plt.figure()
# plt.plot(t_array, (180/np.pi)*x_t[0,:], Label='Pomik $\\varphi_{1}$') #, label='Odziv $\\varphi_{'+str(i+1)+'}$')
# plt.plot(t_array, (180/np.pi)*x_t[1,:], Label='Pomik $\\varphi_{2}$') 
# plt.plot(t_array, (180/np.pi)*x_t[2,:], Label='Pomik $\\varphi_{3}$') 
# plt.plot(t_array, (180/np.pi)*x_t[3,:], Label='Pomik $\\varphi_{4}$') 
# plt.plot(t_array, (180/np.pi)*x_t[4,:], Label='Pomik $\\varphi_{5}$')
# plt.plot(t_array, (180/np.pi)*x_t[5,:], Label='Pomik $\\varphi_{6}$') 
# plt.xlabel('Čas [s]')
# plt.ylabel('Zasuk [°]')
# plt.legend(loc="upper right")
# plt.grid()
# plt.show()

# print((180/np.pi)*max(x_t[0,:]), (180/np.pi)*min(x_t[0,:]))

