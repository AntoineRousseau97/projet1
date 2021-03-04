import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import timeit


def V(x):
    return 4*((1/x)**12-(1/x)**6)

def g(E,x):
    v = V(x)
    return (E - v)**(1/2), v

def x_plusminus(E):
    x_minus = 2**(1/6) * (np.sqrt(E+1)/E-1/E)**(1/6)
    x_plus = 2**(1/6) * (-((np.sqrt(E+1)+1)/E))**(1/6)
    return x_plus, x_minus

def f(E, gamma, n):
    x_plus, x_minus = x_plusminus(E) 
    # x_plus = 2**(1/6) * (np.sqrt(E+1)/E-1/E)**(1/6)
    # x_minus = 2**(1/6) * (-((np.sqrt(E+1)+1)/E))**(1/6)
    integral, error = integrate.quad(lambda x :g(E,x)[0], x_minus, x_plus)
    return gamma * integral - (n+1/2)*3.141592

def bissection(gamma, n, error):
    x1 = -1
    x2 = -0.0000000001
    fx1 = f(x1, gamma, n)
    fx2 = f(x2, gamma, n)
    N = 0
    
    while np.abs(x2-x1) > error:
        x_prime = 0.5*(x1 + x2)
        new_f = f(x_prime, gamma, n)
        
        if new_f == 0: 
            break
        elif new_f*fx1 > 0:
             x1 = x_prime
             fx1 = new_f
        elif new_f*fx2 > 0:
            x2 = x_prime
            fx2 = new_f

        N += 1
            
    return 0.5*(x1+x2), N

def secante(gamma, n, error):
    x0 = -1 
    x1 = -0.0000000001
    fx0 = f(x0, gamma, n)
    fx1 = f(x1, gamma, n)
    N = 0
    
    while np.abs(x1 - x0) > error:
        N += 1
        x_next = x1 - fx1*((x1-x0)/(fx1-fx0))
        x0 = x1
        fx0 = fx1
        x1 = x_next
        fx1 = f(x1, gamma, n)
        
    return x1, N

def time_algorithm(algo_type, target):
    E_list = []
    start_time = timeit.default_timer()
    
    for i in range(20):
        if algo_type == "bissection":
            E = [bissection(150, i, target)[0]]
                   
        elif algo_type == "secante":
            E = [secante(150, i, target)[0]]
        
        E_list += [E]
        
        end_time = timeit.default_timer()
            
    return (end_time-start_time)

print(time_algorithm("bissection", 1e-8))
print(time_algorithm("secante", 1e-8))
