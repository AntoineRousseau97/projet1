import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

avogadro = 6.022e23 #particle/mol
c = 299792458
me = 0.51099895000e6 #eV
mp = 938.54e6 #eV
mp_gram = mp*1.78266192162790e-36
re = 2.8179403227e-15 #m

list_atom_H20 = [1, 8]
list_fraction_H20 = [0.11894, 0.888194]
list_molar_mass_H20 = [1.0079, 15.999]
rho_water = 1 #g/cm^3
I_H20 = 75 #eV

list_atom_bone = [1, 6, 7, 8, 12, 15, 16, 20]   
list_fraction_bone = [0.063984, 0.278, 0.027, 0.410016, 0.002, 0.07, 0.002, 0.147]
list_molar_mass_bone = [1.007, 12.01, 14.007, 15.999, 24.305, 30.974, 32.06, 40.078]
rho_bone = 1.85 #g/cm^3
I_bone = 91.9 #eV


def calc_ne(list_atom_no, list_fraction, list_molar_mass, rho):
    i = 0 
    ne = 0
    for atom_no in list_atom_no:
        nb_electron = atom_no*(list_fraction[i]/list_molar_mass[i])*avogadro
        #ne += (rho/list_molar_mass[i])*nb_electron
        ne += rho*nb_electron
        i += 1
    return ne

def plot_Scol(ne, I, figname):
    v = np.linspace(2e6, 0.99999*c, 10000)  
    beta = v/c
    gamma = 1/np.sqrt(1-beta**2)
    E = 1e6*(gamma-1)*mp_gram*c**2  
    #print(E.index(1e-3))
    #print(np.where())
    #E = (gamma-1)*mp*c**2
    Te_max = (2 * me * c**2 * (gamma**2 - 1)) / (1 + 2 * gamma * (me / mp) + (me / mp)**2)
    Scol = 2*np.pi*(re**2)*me*(c**2)*ne*(1/(beta**2))*(np.log(2*me*(c**2)*(beta**2)*(gamma**2)*Te_max/(I**2))-2*(beta**2))
    f = plt.figure()
    plt.plot(E, Scol, "k")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Énergie (eV)")
    plt.ylabel(r"Pouvoir d’arrêt collisionnel ($MeV \cdot cm^2 / g$)")
    plt.show()
    f.savefig(figname, dpi = 600, bbox_inches = "tight")
    
ne_H20 = calc_ne(list_atom_H20, list_fraction_H20, list_molar_mass_H20, rho_water)
ne_bone = calc_ne(list_atom_bone, list_fraction_bone, list_molar_mass_bone, rho_bone)

print(ne_H20)
print(ne_bone)  
    
plot_Scol(ne_H20, I_H20, "Scol_H20.png")
plot_Scol(ne_bone, I_bone, "Scol_bone.png")


#======================================================================================
# fraction_H = 0.11894 #g_H / gH20
# fraction_O = 0.888106 #g_O / gH207
# molar_mass_H = 1.0079
# molar_mass_O = 15.999
# molar_mass_H20 = 2*molar_mass_H + molar_mass_O

# nb_electron_H = (fraction_H/molar_mass_H)*avogadro*1
# nb_electron_O = (fraction_O/molar_mass_O)*avogadro*8
# nb_electron_H20 = nb_electron_H + nb_electron_O

# #ne_water = (rho_water/molar_mass_H20)*nb_electron_H20 
#======================================================================================


    
# def n_e(rho):
#     return (rho * N * Z * Int_fif) / (2 * np.pi * r_e**2 * m_e * c**2 * 1/beta**2 * (np.ln((2 * m_e**2 * beta**2 * gamma**2 * T_max) / I**2) - 2 * beta**2))




