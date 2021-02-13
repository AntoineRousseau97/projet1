# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 18:10:19 2021

@author: Gabriel
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

rho =
N = 
Z = 1
r_e = 2.8179403227 * 10**(-15)
m_e = 511000
c = 299792458
m_p = 938272000
beta = 
gamma = 
I = 
T_max = (2 * m_e * c**2 * (gamma**2 - 1)) / (1 + 2 * gamma * (m_e / m_p) + (m_e / m_p)**2)
T = 

def int_fif(T):
    return 
    
def n_e(rho):
    return (rho * N * Z * Int_fif) / (2 * np.pi * r_e**2 * m_e * c**2 * 1/beta**2 * (np.ln((2 * m_e**2 * beta**2 * gamma**2 * T_max) / I**2) - 2 * beta**2))