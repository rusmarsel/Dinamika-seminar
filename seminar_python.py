# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 14:23:13 2020

@author: Marsel Rus
"""

# import seminar_knjiznica as sk
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
n_max = 80  # [/]

# Omega

w = 2*np.pi/T 

# Časovni intervali

t_min = 0.       # [s]
t_max = T     # [s]
dt = 0.0001     # [s]

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
            b_j[j] = 2*M0*(-np.cos(j*np.pi)/(j*np.pi))
            
    return a_j, b_j

    
def primerjava(t, moment2, a_j, b_j, w, show = False):

    M_approx = np.ones_like(t_array) 
    
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

        
# def mom_four(t, b_j):
#     M_four = np.ones_like(t_array) 
        
#     return M_approx


def mom_four(t_interval, b_j, show=False):
    '''Funkcija izračuna vrednosti momenta za dani časovni interval'''
    M_four = np.ones_like(t_interval) 
    for j in range(len(b_j+1)):
        if j==0:
            M_four *= 0
        else:
            M_four += b_j[j]*np.sin(j*w*t_interval)           
    if show:
        plt.figure()
        plt.plot(t_interval, M_four)
        plt.show()    
    return M_four


# Eksaktna funkcija vzbujevalnega momenta (za preverjanje, koda iz vaj)
M_exact = moment2(t_array, T, M0, show = False)  # Za prikaz -> show = True
a_j, b_j = koef_bj(M0, n_max)
primerjava(t_array, M_exact, a_j, b_j, w, show = False)  # Za prikaz -> show = True

# =============================================================================
# Modalna obremenitev
# =============================================================================
# To matriko uporabim kot zapis faktorjev pri prehodu v modalne koordinate
def M_vzb_modalna(N):
   '''Izračuna matriko modalnih koeficientov momenta'''
   M_vzb = np.zeros((N,1))
   for i in range(N):
       if i < 2:
           M_vzb[i,0] = 1
   M_vzb = np.matmul(np.transpose(fi_last), M_vzb)    # fi_last - matrika lastnih vektorjev
   M_vzb = np.matmul(np.linalg.inv(M_modalna), M_vzb) # M_modalna - matrika modalnih mas
   return M_vzb


m_breme = M_vzb_modalna(N)


# def b_j_modalni(P_S, b_j):
#     h_j=np.zeros(n_max)
#     for j in range(n_max):
#         if j==0:
#             h_j = h_j
#         else:
#             h_j[j] = b_j[j]*m_breme[P_S,0]/las_vr[P_S]   # bj/woi^2
#     return h_j


def b_js(t):
    moment = mom_four(t, b_j)    
    four = np.ones((N,len(t)))
    
    # sinus = np.ones_like(t)    
    # if (j*w) < las_frek[P_S]: # primerja, dokler je j*w < w0i
    #         sinus += np.sin(j*w*t_int)
    #     else:
    #         sinus += np.sin(j*w*t_int - np.pi)  
    
    for i in range(len(t)):
        for p in range(N):
            for j in range(n_max):
                    if (j*w) < las_frek[p]: # primerja, dokler je j*w < w0i
                        four[p,i] = (np.sin(j*w*t[i])*(moment[i]*m_breme[p])*(1/(1-((j*w)/las_frek[p])**2)))/las_vr[p]
                    else:
                        four[p,i] = (np.sin(j*w*t[i]-np.pi)*(moment[i]*m_breme[p])*(1/(1-((j*w)/las_frek[p])**2)))/las_vr[p]         
        
    # hij = b_j_modalni(P_S, b_j)
    # sinus = np.ones_like(t)
    
    # for i in range(len(t)):
    #     for j in range(n_max):
    #         if j == 0:
    #             hij=hij
    #         else:
    #             sinus[i] += hij[j]*np.sin(j*w*t[i])
    return four

# hjs = b_js(t_array,0)
# plt.plot(t_array, hjs)

# bj = np.ones((N,len(t_array)))    
# for j in range(N):
#     bj[j,:] = b_js(t_array, j)
# def breme (t_int, P_S):
#     '''Izračuna obremenitev v časovnem intervalu za vsako prostostno stopnjo posebej'''
#     mom_fourier = mom_four(t_int, b_j)
#     hi_modalna = np.ones((P_S,len(t_int)))    
#     for i in range(len(t_int)):        
#         br = m_breme * mom_fourier[i]    
#         hi_t_mod = np.dot(np.transpose(fi_last), br)            
#         for j in range(P_S):            
#             hi_modalna[j,i] = hi_t_mod[j]*(1/M_modalna[j,j])       #     np.dot(np.linalg.inv(MasMx), hi_t_mod)[j]  
#     return hi_modalna


# def h_i_t(t_interval, b_j, P_S, show = False):
#     '''Izračun vektorja momentov v modalnih koordinatah - IZRAČUNANO LE NA PROSTOSTNI STOPNJI P_S (med 0 in N-1)'''    
#     br = breme(t_interval, N)    
#     hi_mi = np.ones(len(t_interval))    
#     for i in range(len(t_interval)):
#         hi_mi[i] = br[P_S,i]        
#     if show:
#         plt.figure()
#         plt.plot(t_interval, hi_mi)
#         plt.show()        
#     return hi_mi
# # h_i_t(t_array, b_j, 0, show=True)


# def beta_j(P_S):
#     '''Izračun faktorja beta na prostostni stopnji'''
#     Beta=0    
#     for j in range(n_max):
#         Beta += 1/(1-((j*w)/las_frek[P_S])**2)     # las_frek nosi vrednosti w0i   
#     return Beta


# def sin_fi_t(t_int, P_S, show=False):
#     '''Izračun sinusa sin(jwt-fi) v časovnem intervalu za eno prostostno stopnjo'''
#     sinus = np.ones_like(t_int)    
#     for j in range(len(b_j+1)):  
#         if (j*w) < las_frek[P_S]: # primerja, dokler je j*w < w0i
#             sinus += np.sin(j*w*t_int)
#         else:
#             sinus += np.sin(j*w*t_int - np.pi)            
#     if show:
#         plt.figure()
#         plt.plot(t_int, sinus)
#         plt.show()    
#     return sinus


# def b_jt(t, P_S):
#     bj = np.ones_like(t)
#     for i in range(n_max):
#         bj += b_j[i]*np.sin(j*w*t)*m_breme[P_S]       
#     return bj


# def eta_t(t, P_S):
#     '''Izračun vektorja eta za vse prostostne stopnje v časovnem intervalu t'''    
#     # eta = np.ones((N,len(t)))
#     eta = b_js(t)            
#     for j in range(P_S):        
#         # Be = beta_j(j)
#         # sinus = sin_fi_t(t, j)
#         # hi = h_i_t(t, b_j, j)
        
#         # for i in range(len(t)):
#         #     eta[j, i] = Be*sinus[i]*bj[j,i] # las_vrvzame vrednost w0i^2
        
#         # plt.figure()
#         # plt.plot(t,eta[j,:])
#         # plt.show()
#       return eta

x_t = np.ones((N,len(t_array)))
Vec_Eta = b_js(t_array)

for i in range(N):
    for j in range(len(t_array)):
        eta = Vec_Eta[:,j]
        
        x_t[:,j] = np.dot(fi_last, eta)        

plt.figure()
plt.plot(t_array, (180/np.pi)*x_t[0,:])
plt.xlabel('Čas [s]')
plt.ylabel('Zasuk [°]')
plt.show()
