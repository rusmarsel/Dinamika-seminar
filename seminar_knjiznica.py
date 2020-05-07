# -*- coding: utf-8 -*-
"""
Created on Sat May  2 23:18:32 2020

@author: Mars
"""

import numpy as np
import matplotlib.pyplot as plt

from sympy import *

M0 = 1150   # [Nm]
T = 0.17    # [s]
w = 2*np.pi/T
cleni = 300

a, b, t, tau, M, pi, j = symbols('a, b, t, tau, M, pi, j')

f_mom1 = 2*M0*(t/tau)
f_mom2 = 2*M0*(t/tau) - 2*M0

# I_1 = integrate(f_mom1, (t, a, b))
# I_2 = integrate(f_mom2, (t, a, b))
# podatki = {a: 0, b: T/2, tau: T, M: M0}
# podatki2 = {a: T/2, b: T, tau: T, M: M0}
# I1_analiticno = I_1.subs(podatki).evalf() 
# I2_analiticno = I_2.subs(podatki2).evalf()   
# print(I1_analiticno)
# print(I2_analiticno)
# a_0 = (2/T)*(I1_analiticno + I2_analiticno)


# print(a_0, "Srednja vrednost funkcije je 0, člen zanemarim","\n")

# I_1 = integrate(f_mom1*cos(j*2*np.pi/T*t), (t, a, b))
# I_2 = integrate(f_mom2*cos(j*2*np.pi/T*t), (t, a, b))
# podatki = {a: 0, b: T/2, tau: T, M: M0}
# podatki2 = {a: T/2, b: T, tau: T, M: M0}
# I1_analiticno = I_1.subs(podatki).evalf() 
# I2_analiticno = I_2.subs(podatki2).evalf()   
# # print(I1_analiticno)
# # print(I2_analiticno)
# a_j = (2/T)*(I1_analiticno + I2_analiticno)


# print(a_j, "Člena integrala se med seboj odštejeta. Člen je tako enak 0, lahko ga zanemarim.", "\n")

   
I_1 = integrate(f_mom1*sin(j*t*(2*pi)/tau), (t, a, b))
I_2 = integrate(f_mom2*sin(j*t*(2*pi)/tau), (t, a, b))
podatki = {a: 0, b: T/2, tau: T, M: M0}
podatki2 = {a: T/2, b: T, tau: T, M: M0}
I1_analiticno = I_1.subs(podatki).evalf() 
I2_analiticno = I_2.subs(podatki2).evalf()   
# print(I1_analiticno)
# print(I2_analiticno)
b_j = (2/T)*(I1_analiticno + I2_analiticno)


print(b_j)

# Vzamem samo uporabne j-je, ki so med 0 in N. Ta člen nesem tudi v glavni program.
# b_j = 11.7647058823529*((-31.1147913744655*cos(3.14159265358979*j)/j + 9.90414570103852*sin(3.14159265358979*j)/j**2)+(-31.1147913744655*cos(3.14159265358979*j)/j - 9.90414570103852*sin(3.14159265358979*j)/j**2 + 9.90414570103852*sin(6.28318530717959*j)))



# print(moment1(t_h,M0,t0))

# def trapezno(f, h):
#     '''Funkcija za trapezno interpolacijo integrala'''
#     trapez = (f[0] + f[-1])*h/2
#     return trapez

# def odziv(t, M0, t0):
#     '''Funkcija za izračun odziva sistema na impulzno motnjo'''
#     x_odziv = trapezno(moment1(t_h, M0, t0), h)
#     return x_odziv

# print(odziv(0, M0, t0))


# Korak trapezne interpolacije
# h = .1
# t_h = np.linspace(0,t0,100)
bee = np.zeros(cleni)
for j in range(cleni):
    if j==0:
        bee[j]=bee[j]
    else:
        bee[j] = 1150*(-np.cos(j*np.pi)/(j*np.pi) + np.sin(j*np.pi)/(j**2*np.pi**2)+(-np.cos(j*np.pi)/(j*np.pi) - np.sin(j*np.pi)/(j**2*np.pi**2) + np.sin(2*j*np.pi)/(j**2*np.pi**2)))
    
t0 = 0
tmax = T
dt = .0001
t_array = np.arange(t0,tmax,dt)

M_approx = np.ones_like(t_array)
for j in range(cleni):
    M_approx += bee[j]*np.sin(j*w*t_array)
    
plt.figure()
plt.plot(t_array, M_approx)
plt.show()