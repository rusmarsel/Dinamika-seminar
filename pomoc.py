

# def beta(N, n_max):
        
#     Beta = np.ones(N)
    
#     for i in range(N):
#         for j in range(n_max):
            
#             Beta[i] += 1/(1-((j*w)/las_frek[i])**2)
        
#     return Beta


# def X_3(N, b_j, n_max):
    
#     X = np.ones(N) 
    
#     for i in range(N):
#         for j in range(n_max):            
#             X[i] += b_j[j] / (MasMx[i,i]*(las_frek[i])**2)
    
#     return X


# def odziv(t, n_max, N, X_3, beta):
    
#     Beta = beta(N,n_max)
#     X3 = X_3(N, b_j,n_max)
    
#     odz = np.ones_like(t_array) 
    
#     for i in range(N):
#         for j in range(n_max):
#             odz += np.sin(j*w*t)*Beta[i]*X3[i]*np.sin(j*w*t)
    
#         plt.figure()
#         plt.plot(t, odz)
#         plt.show()
#     return odz
        
# odziv(t_array, n_max, N, X_3, beta)

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

# def vectorfield(t, y, N, a_j, b_j, M_modalna, K_modalna):
    
# Vectorfield function takes a state space vector as input and returns its derivative
# at given time instance. Arbitrary arguments can be passed for the derivative calculation
#
# inputs
# y:        state-space vector; in the current case of prism and cylinder it consists of
#           y = [y[0], y[1], y[2], y[3]] = [x, dx/dt, s, ds/dt] 
# t:        time instance; this is automatically provided by odeint function
# params:   a tuple of additional arguments
#
# outputs
# dy:       a derivative of a state-space vector; in the current case it consists of
#           dy = [dx/dt, d^2x/dt^2, ds/dt, d^2s/dt^2]

    
    #define matrix A according to devised equilibrium eqs