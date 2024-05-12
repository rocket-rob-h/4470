# # Find R(i) constant for each molecule species
# R_O2 = R_gas(Mw_O2)
# print(f'The R Gas constant for Oxygen is {R_O2:.2f} J/kg*K')
# R_N2 = R_gas(Mw_N2)
# print(f'The R Gas constant for Nitrogen is {R_N2:.2f} J/kg*K')
#
# CV_diatomic_O2 = CV_diatomic(R_O2, thetaV_O2,thetaT_O2)
# print(f'CV_diatomic Oxygen is {CV_diatomic_O2:.2f}')


# def CV_diatomic(R,theta_v, T):
#     outside = ((3/2)*R)+R
#     num = (theta_v**2/T) * np.exp((theta_v/T))
#     denom = np.exp((theta_v/T)-1)**2
#     return outside + ((num/denom) * R)





#Preamble
import numpy as np
#Question 1

#Given
V_inf = 5200 # m/s
T_inf = 1000 # K
P_inf = 500 # Pa
Diameter = 0.075 # m
mol_N2 = 0.79
mol_O2 = 0.21
Mw_O2 = 32
Mw_N2 = 28
Ru = 8314
U2 = 0 # 0m/s

#Kruger
thetaV_O2 = 2270 # K
thetaT_O2 = 11390 #K
thetaV_N2 = 3390 # K

''' Post Shock Temperature'''
''' Post Shock Density'''
'''Sensible Energy'''
'''Post Shock Velocity'''
'''Composition by Mole Fraction'''

def R_gas(Mw):
    return 8314/Mw


def Ideal_gas_law_rho(Pressure,R, Temperature):
    return Pressure / (R * Temperature)

rho1 = Ideal_gas_law_rho(P_inf)
print(f'rho1{rho1:.2f}')

def P2(P1, rho1, U1, rho2,U2):
    return P1 + rho1 * U1**2 - rho2 * U2**2

def h2(h1, U1, U2):
    return h1 + (0.5 * U1**2) - (0.5 * U2**2)

def U2(rho1, U1, rho2):
    return (rho1 * U1)/(rho2)

def rho2(P2, Z2, R, T2):
    return (P2) / (Z2 * R * T2)


'''Free steam can be modelled using perfect gas - surely we have to find the real CP initial?'''
# 500 Pa = 0.00493462 atm

#Cengel Appendix says air
Cp_1000K = 1142
Cv_1000K = 855
gamma_1000 = 1.336

h1 = Cp_1000K *




















