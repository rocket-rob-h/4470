#Preamble
import numpy as np
#Question 1

#Given
U_inf = 5200 # m/s
T_inf = 1000 # K
P_inf = 500 # Pa
Diameter = 0.075 # m
mol_N2 = 0.79
mol_O2 = 0.21
Mw_O2 = 32
Mw_N2 = 28
Ru = 8314
U2 = 0 # 0m/s
Mw_air = 28.97 #g/mol


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

R_air = R_gas(Mw_air)
print(f'R for Air is {R_air:.2f}')

def Ideal_gas_law_rho(Pressure,R, Temperature):
    return Pressure / (R * Temperature)

rho1 = Ideal_gas_law_rho(P_inf, 287 ,T_inf)
print(f'rho1 is {rho1:.8f} kg/m^3')

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

h1 = Cp_1000K * T_inf
print(f'h1 Enthalpy 1 is {h1:.2f} J/kg.K')

def Pmom(P1,rho1,U1):
    return P1 + (rho1 * U1)

PMOM = Pmom(P_inf, rho1, U_inf**2)
print(f'Pmom is {PMOM:.2f} Pa')

PMOM_atm = PMOM/101325
print(f'PMOM in atm is {PMOM_atm} atm')

def Htot(h1,U1):
    return h1 + (0.5 * U1**2)

HTOT = Htot(h1,U_inf)
print(f'Enthalpy HTOT is {HTOT/1000000:.2f} MJ/kg.K')

'''Iteration 1'''

# PMOM = P2_1
# HTOT = h2_1

# Now to Fig 1.176c to find T2 and z1
T2_1 = 5750 # K
z1_1 = 1.313

# Rho 2 now with these variables

rho2_1 = rho2(PMOM,z1_1,R_air, T2_1)
print(f'Density 2 Iteration 1 rho2 is {rho2_1:8f} kg/m^3')

# Find U2
U2_1 = U2(rho1,U_inf,rho2_1)
print(f'First iteration of U2 is {U2_1:.2f} m/s')

'''Iteration 2 at this U2'''

P2_2 = PMOM - rho2_1 * U2_1**2
print(f'P2 Iteration 2 is {P2_2:.2f} Pa')
P2_2_atm = P2_2/101325
print(f'P2_2_atm is {P2_2_atm:.3f}')

# Now to Fig 1.176c to find T2_2 and z2_2

h2_2 = HTOT - (0.5 * U2_1**2)
print(f'h2_2 Enthalpy Iteration 2 is {h2_2/1000000:.2f} MJ/kg.K')

T2_2 = 5687 # K
z2 = 1.314

rho2_2 = rho2(P2_2,z2, R_air,T2_2)
print(f'rho2_2 Density 2 Iteration 2 is {rho2_2:6f} kg/m^3')

U2_2 = U2(rho2_1,U2_1,rho2_2)
print(f'Converged U2 is {U2_2:.2f} m/s')


'''Iteration 3 for U2_2'''
P2_3 = PMOM - (rho2_2 * U2_2**2)
print(f'P2_3 is {P2_3:.2f} Pa')

PMOM_3 = P2_3/101325
print(f'PMOM_3 in atm is {PMOM_3}')

h2_3 = HTOT - (0.5 * U2_2**2)
print(f'h2_3 Enthalpy Iteration 3 is {h2_3/1000000:.3f} MJ/kg.K')

print()
print('The Pressure and Enthalpy have converged to the full accuracy of these 70 year old graphs after 2 iterations.')

print(f'The post shock density is {rho2_2:.7f} kg/m^3')
print(f'The post shock pressure is {PMOM_3:.6f} Pa')
print(f'The post shock velocity is {U2_2:2f} m/s')
print(f'The post shock temperature is {T2_2:.2f} K')

print('Now find sensible enthalpy - is that what I found?')
print(f'The sensible post shock enthalpy is {h2_3} J/kg.K')
print(f'Now to find the dissociated composition at this Temperature and Pressure')

'''CEA Values'''


# Nitrogen

N_HoT_0K = 111543.3
# N_neg_FT_Ho_1000K  = 37658.6
# N_HoT_1000K = 116511.5
N_neg_FT_Ho_5687K = 263171.6
N_HoT_5678K = 140700

# N2 Nitrogen

N2_HoT_0K = -2072.3
# N2_neg_FT_Ho_1000K  = 41307.7
# N2_HoT_1000K = 5129.8
N2_neg_FT_Ho_5687K = 346000.8
N2_HoT_5678K = 46200


# O2 Oxygen

O2_HoT_0K = -2074.7
# O2_neg_FT_Ho_1000K  = 50691.7
# O2_HoT_1000K = 5426.9
O2_neg_FT_Ho_5687K = 371000.1
O2_HoT_5678K = 50200

# O Oxygen

O_HoT_0K = 57949.1
# O_neg_FT_Ho_1000K  = 39459.6
# O_HoT_1000K = 63108.4
O_neg_FT_Ho_5687K = 274002.6
O_HoT_5678K = 86820.2

# Nitrogen Oxygen

NO_HoT_0K = 19403.0
# NO_neg_FT_Ho_1000K  = 51865.5
# NO_HoT_1000K = 31096.3
NO_neg_FT_Ho_5687K = 374003.2
NO_HoT_5678K = 68700.9

joules = 4.184

print(f'P2_2_atm is {P2_2_atm:.6f}')

print()
'''Gibbs Boys'''
N2FT = N2_HoT_0K - N2_neg_FT_Ho_5687K
print(f'N2FT is {N2FT} cal and in Joules is {N2FT*joules:.2f} J')

NFT = N_HoT_0K - N_neg_FT_Ho_5687K
print(f'NFT is {NFT} cal and in Joules is {NFT*joules:.2f} J')

OFT = O_HoT_0K - O_neg_FT_Ho_5687K
print(f'OFT is {OFT:.2f} cal and in Joules is {OFT*joules:.2f} J')

O2FT = O2_HoT_0K - O2_neg_FT_Ho_5687K
print(f'O2FT is {O2FT} cal and in Joules is {O2FT*joules:.2f} J')

NOFT = NO_HoT_0K - NO_neg_FT_Ho_5687K
print(f'NOFT is {NOFT} cal and in Joules is {NOFT*joules:.2f} J')

N2FT_J = N2FT*joules
NFT_J = NFT*joules
OFT_J = OFT*joules
O2FT_J = O2FT*joules
NOFT_J = NOFT*joules

print()

GibbsNitrogen = (2 * NFT_J) + (-1 * N2FT_J)
print(f'Gibbs for N2 to N is {GibbsNitrogen:.2f} J/mol')

GibbsOxygen = (2 * OFT_J) + (-1 * O2FT_J)
print(f'Gibbs for O2 to O is {GibbsOxygen:.2f} J/mol')

GibbsNO = (1 * NOFT_J) + (-1 * NFT_J) + (-1 * OFT_J)
print(f'Gibbs for (-1)N + -1(O) to (1) NO is {GibbsNO:.2f} J/mol')

print()

def Gibbs(DeltaGibbs, T):
    exponent = 8.314 * T
    return np.exp(-DeltaGibbs / exponent)

KN = Gibbs(GibbsNitrogen,T2_2)
print(f'KN is {KN:.7f} atm')

KO = Gibbs(GibbsOxygen,T2_2)
print(f'KO is {KO:.7f} atm')

KNO = Gibbs(GibbsNO,T2_2)
print(f'KNO is {KNO:.7f} atm')

#DALTONS LAW



