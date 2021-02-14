import numpy as np
import matplotlib.pyplot as plt

Ti = 150e6
rho_H2O = 1000
rho_bone = 1.85 
I_H20 = 75 #eV
I_bone = 91.9 #eV
c = 299792458
#me = 0.51099895000e6 #eV
me = 9.1094E-31*6.242e18 #eV
mp = 1.6727E-27*6.242e18 #eV 
re = 2.8179403227e-15 #m

#ne_H20 = 3.385161894888596e23
#ne_bone = 5.9065706308435045e23

#ne_H20 = 3.342916491E29 #é/m^3

def f(T, ne, I, rho):
    gamma = T/(mp*c**2) + 1
    beta = np.sqrt(1- 1/gamma**2)
    Te_max = (2 * me * c**2 * (gamma**2 - 1)) / (1 + 2 * gamma * (me / mp) + (me / mp)**2)
    Scol = 2*np.pi*(re**2)*me*(c**2)*ne*(1/(beta**2))*(np.log(2*me*(c**2)*(beta**2)*(gamma**2)*Te_max/(I**2))-2*(beta**2))
        
    return 1/(Scol/rho)
#===================================================================
# #Énumérer les constantes importantes
# r_e = 2.8179E-15 #m
# m_e = 9.1094E-31 #kg
# c = 3.00E8 #m^2
# n_e = 3.342916491E29 #é/m^3
# I = 75 #eV
# m_p = 1.6726E-27 #kg
# conv = 6.242E18 #eV
# Ti=15010**6 #eV
# rho = 1000 #kg/m^3



# def S_col(T):
#     g = (T / (conv * m_p * c  2)) + 1
#     b2 = ((((T / (conv * m_p * c  2)) + 1)  2) - 1) / (((T / (conv * m_p * c  2)) + 1)  2)
#     Te_max = ((2 * conv * m_e * c  2) * ((((T / (conv * m_p * c  2)) + 1)  2) - 1)) / (
#             1 + (2 * (T / (conv * m_p * c  2) + 1) * (m_e / m_p)) + ((m_e / m_p)  2))
#     return rho/(2 * m.pi * (r_e  2) * m_e * conv * (c  2) * n_e * (1 / (b2)) * (np.log(((conv * 2 * m_e * (c  2) * (b2) * (g  2) * Te_max) / (I  2)) - 2 * (b2))))
#====================================================================
    
def trapeze(N, a, b, ne, I, rho): 
    h = (b-a)/N
    #méthode du trapèze
    s = 0.5*f(a, ne, I, rho) + 0.5*f(b, ne, I, rho)
    
    for k in range(1,N):
        s += f(a+k*h, ne, I, rho)
        
        #print("Valeur calculée = ", h*s, "( N=",N,")")
    return h*s
    
def simpson(N, a, b, ne, I, rho):
    h = (b-a)/N
    #méthode de Simpson
    s1 = 0.0
    for k in range(1,N,2):
        s1 += f(a+k*h, ne, I, rho)
    s2 = 0.0
    for k in range(2,N,2):
        s2 += f(a+k*h, ne, I, rho)
        
    return (f(a, ne, I, rho)+f(b, ne, I, rho)+4.0*s1+2.0*s2)*h/3.0

previous_value = 0
N = 0
percent_diff = 100
liste_value = []
liste_N = []


#while percent_diff > 0.01:
# for i in range(100):    
#     N += 1
#     value = trapeze(N, 3e6, Ti, ne_H20, I_H20, rho_H2O)
#     liste_value += [value]
#     liste_N += [N]
#     diff = value - previous_value
#     percent_diff = 100*diff/value
#     #print(value)
#     previous_value = value
    
# print(N)
# plt.plot(liste_N, liste_value)
# plt.show()

# previous_value = 0
# N = 0
# percent_diff = 100
# liste_value = []
# liste_N = []

# #while percent_diff > 0.01
# for i in range(100):
#     N += 1
#     value = simpson(N, 3e6, Ti, ne_H20, I_H20, rho_H2O)
#     liste_value += [value]
#     liste_N += [N]
#     #diff = value - previous_value
#     #percent_diff = 100*diff/value
#     #print(value)
#     previous_value = value
    
# print(N)
# plt.plot(liste_N, liste_value)

# plt.show()

print(trapeze(1000, 3e6, Ti, ne_H20, I_H20, rho_H2O))
# print(simpson(100, 1e-16, Ti, ne_H20, I_H20, rho_H2O))
        


