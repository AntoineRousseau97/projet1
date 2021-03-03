import numpy as np
from scipy import integrate


def g(E,x):
    v = 4*((1/x)**12-(1/x)**6)
    
    return (E - v)**(1/2)

def f(E, gamma, n):
    x_plus = 2**(1/6) * (np.sqrt(E+1)/E-1/E)**(1/6)
    x_minus = 2**(1/6) * (-((np.sqrt(E+1)+1)/E))**(1/6)
    integral, error = integrate.quad(lambda x :g(E,x), x_plus, x_minus)
    return gamma * integral - (n+1/2)*3.141592

def epsilon(gamma, n, error):
    x1 = -1
    x2 = -0.0000000001
    fx1 = f(x1, gamma, n)
    fx2 = f(x2, gamma, n)
    N = 0
    
    while np.abs(x2-x1) > error:
        x_prime = 0.5*(x1 + x2)
        new_f = f(x_prime, gamma, n)
        
        if new_f == 0: 
            print("wtf")
            break
        elif new_f*fx1 > 0:
             x1 = x_prime
             fx1 = new_f
        elif new_f*fx2 > 0:
            x2 = x_prime
            fx2 = new_f

        N += 1
            
    return 0.5*(x1+x2), N
            
    
answer = epsilon(150, 1, 1e-8)
print(answer[0])
print(answer[1])
print(f(answer[0], 150, 1))
    
            