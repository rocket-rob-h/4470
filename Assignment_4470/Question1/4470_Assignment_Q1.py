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





# Change in enthalpy between the pre shock and post shock state between T-1000K and T 5600 K
# The sensible enthalpy of a mixture can be obtained from the following:
# a) the sensible enthalpy for each species as given by the formulas of
# statistical mechanics for example, Eqs. (11.62), (11.63), and (11.92); and
# b) knowledge of the equilibrium composition described in terms of pi , Xi ,
# hi , or ci .
# 2. The zero-point energy can be treated as an effective value by using the
# heats of formation at absolute zero in its place. Therefore, Eq. (11.105) or
# (11.106) can be construed as the enthalpy of a gas mixture.

''' 
H'f,
col/mole
T0 -2072.3

-225159.0

'''
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

# Products - Reactants

#Correct around the different way

# GibbsNitrogen = (-1 * N2FT_J) + (2 * NFT_J)
# print(f'Gibbs for N2 to N is {GibbsNitrogen:.2f} J/mol')
#
# GibbsOxygen = (-1 * O2FT_J) + (2 * OFT_J)
# print(f'Gibbs for O2 to O is {GibbsOxygen:.2f} J/mol')
#
# GibbsNO = (-1 * NFT_J) + (-1 * OFT_J) + (1 * NOFT_J)
# print(f'Gibbs for (-1)N + -1(O) to (1) NO is {GibbsNO:.2f} J/mol')

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
print(f'KN is {KN:.5f} atm')

KO = Gibbs(GibbsOxygen,T2_2)
print(f'KO is {KO:.5f} atm')

KNO = Gibbs(GibbsNO,T2_2)
print(f'KNO is {KNO:.5f} atm')

#DALTONS LAW

yes = KN * KO
print(yes)


#
# import numpy as np
# from scipy.optimize import least_squares
#
# # Constants
# KN = 0.01895        # atm
# KO = 185.64018      # atm
# KNO = 0.31424       # atm
# P2_2_atm = 0.432992 # total pressure atm
#
# # Define the function for residuals
# def residuals(vars):
#     P_N2, P_N, P_O2, P_O, P_NO = vars
#     eq1 = (P_N**2 / P_N2) - KN   # Equilibrium for N
#     eq2 = (P_O**2 / P_O2) - KO   # Equilibrium for O
#     eq3 = (KNO * P_N * P_O) - P_NO  # Formation of NO
#     eq4 = (P_N2 + P_N + P_O2 + P_O + P_NO) - P2_2_atm  # Dalton's law
#     eq5 = ((2 * P_N2 + P_N + P_NO) / (2 * P_O2 + P_O + P_NO)) - (0.79 / 0.21)  # Conservation of nuclei ratio
#     return [eq1, eq2, eq3, eq4, eq5]
#
# # Bounds to ensure all pressures are non-negative
# bounds = (0.00000001, 1)
#
# # Initial guess based on known values
# initial_guess = [0.51083, 0.16241, 0.00017813, 0.31302, 0.0061796]
#
# # Perform the least squares optimization
# result = least_squares(residuals, initial_guess, bounds=bounds)
#
# # Print the results
# if result.success:
#     print("Solution found:")
#     print(f"P_N2: {result.x[0]}, P_N: {result.x[1]}, P_O2: {result.x[2]}, P_O: {result.x[3]}, P_NO: {result.x[4]}")
# else:
#     print("Solution not found:", result.message)


# import numpy as np
# from scipy.optimize import fsolve
#
# # Constants
# KN = 0.01895        # atm
# KO = 185.64018      # atm
# KNO = 0.31424       # atm
# P2_2_atm = 0.432992 # total pressure atm
#
# # Define the equations for fsolve
# def equations(vars):
#     P_N2, P_N, P_O2, P_O, P_NO = vars
#     eq1 = (P_N**2 / P_N2) - KN
#     eq2 = (P_O**2 / P_O2) - KO
#     eq3 = KNO * P_N * P_O - P_NO
#     eq4 = P_N2 + P_N + P_O2 + P_O + P_NO - P2_2_atm
#     eq5 = (P_N2 + P_N + P_NO) / (P_O2 + P_O + P_NO) - (0.79 / 0.21)
#     return [eq1, eq2, eq3, eq4, eq5]
#
# # Provided exact answers as initial guesses
# initial_guesses = [0.51083, 0.16241, 0.00017813, 0.31302, 0.0061796]
#
# # Solve the equations using fsolve
# solution = fsolve(equations, initial_guesses)
#
# # Print the solutions
# print(f"Solution found: P_N2: {solution[0]}, P_N: {solution[1]}, P_O2: {solution[2]}, P_O: {solution[3]}, P_NO: {solution[4]}")


#
# import sympy as sym
# P_O2, P_O, P_N2, P_N, P_NO = sym.symbols('P_O2, P_O, P_N2, P_N, P_NO')
#
# eq1 = sym.Eq(185.6*P_O2 - P_O**2,0)
# eq2 = sym.Eq(0.01895 *P_N2 - P_N**2,0)
# eq3 = sym.Eq(0.314 * P_N*P_O - P_NO,0)
# eq4 = sym.Eq(P_N2 + P_N + P_O2 + P_O + P_NO, P2_2_atm)
# eq5 = sym.Eq(2 * P_N2 + P_N + P_NO) / (2 * P_O2 + P_O + P_NO), (0.79 / 0.21)
#
# result = sym.solve([eq1,eq2,eq3,eq4,eq5],(P_O2, P_O, P_N2, P_N, P_NO))
#
#
# import sympy as sym
#
# # Define symbols
# P_O2, P_O, P_N2, P_N, P_NO = sym.symbols('P_O2 P_O P_N2 P_N P_NO')
# P2_2_atm = 0.432992  # Assuming this is a constant, set to the value previously used
#
# # Define equations
# eq1 = sym.Eq(185.6 * P_O2 - P_O**2, 0)  # Quadratic equation for oxygen
# eq2 = sym.Eq(0.01895 * P_N2 - P_N**2, 0)  # Quadratic equation for nitrogen
# eq3 = sym.Eq(0.314 * P_N * P_O - P_NO, 0)  # Equation for nitric oxide
# eq4 = sym.Eq((P_N2 + P_N + P_O2 + P_O + P_NO), 0.432992)  # Dalton's Law of partial pressures
# eq5 = sym.Eq((2 * P_O2 + P_O + P_NO) / (2 * P_N2 + P_N + P_NO), 0.2658)  # Ratio of nitrogen to oxygen
#
# # eq5 = sym.Eq((2 * P_N2 + P_N + P_NO) / (2 * P_O2 + P_O + P_NO), (0.79 / 0.21))  # Ratio of nitrogen to oxygen
#
# # Attempt to solve the equations
# result = sym.solve([eq1, eq2, eq3, eq4, eq5], (P_O2, P_O, P_N2, P_N, P_NO), dict=True)
# # Output results
# if result:
#     print("Solutions found:")
#     for sol in result:
#         print(sol)
# else:
#     print("No solutions found or solve failed to converge.")


#
# eq1 = sym.Eq(P_O2 + P_O -  P_N2 - P_N + P_NO,0)
# eq2 = sym.Eq(2 * P_O2 + P_O + 2 * P_N2 +P_NO,0.43)
# eq3 = sym.Eq(185.6*P_O2 - P_O**2,0)
# eq4 = sym.Eq(0.01895 *P_N2 - P_N**2,0)
# eq5 = sym.Eq(0.314 * P_N*P_O - P_NO,0)

#
# import sympy as sym
#
# # Define symbols
# P_O2, P_O, P_N2, P_N, P_NO = sym.symbols('P_O2 P_O P_N2 P_N P_NO')
#
# # Equations based on the provided relationships
# eq1 = sym.Eq(7.52 * P_O2 + 3.76 * P_O - 2 * P_N2 - P_N + 2.76 * P_NO, 0)
# eq2 = sym.Eq( P_O2 + P_O + P_N2 + P_NO, 0.43)
# eq3 = sym.Eq(244.64 * P_O2 - P_O**2, 0)
# eq4 = sym.Eq(0.0234 * P_N2 - P_N**2, 0)
# eq5 = sym.Eq(0.227 * P_N * P_O - P_NO, 0)
# #
# # Attempt to solve the equations
# result = sym.solve([eq1, eq2, eq3, eq4, eq5], (P_O2, P_O, P_N2, P_N, P_NO))

# # Check and display the results
# if result:
#     print("Solutions found:")
#     for sol in result:
#         print(sol)
# else:
#     print("No solutions found or solve failed to converge.")

#
# import sympy as sym
#
# # Define symbols
# P_O2, P_O, P_N2, P_N, P_NO = sym.symbols('P_O2 P_O P_N2 P_N P_NO')
#
# # Equations based on the provided relationships
# eq1 = sym.Eq(7.52 * P_O2 + 3.76 * P_O - 2 * P_N2 - P_N + 2.76 * P_NO, 0)
# eq2 = sym.Eq(2 * P_O2 + P_O + 2 * P_N2 + P_NO, 0.43299)
# eq3 = sym.Eq(185.6 * P_O2 - P_O**2, 0)
# eq4 = sym.Eq(0.01895 * P_N2 - P_N**2, 0)
# eq5 = sym.Eq(0.314 * P_N * P_O - P_NO, 0)
#

# P_O2, P_O, P_N2, P_N, P_NO = sym.symbols('P_O2 P_O P_N2 P_N P_NO')
# P2_2_atm = 0.432992  # Assuming this is a constant, set to the value previously used
#
# # Define equations
# eq1 = sym.Eq(185.6 * P_O2 - P_O**2, 0)  # Quadratic equation for oxygen
# eq2 = sym.Eq(0.01895 * P_N2 - P_N**2, 0)  # Quadratic equation for nitrogen
# eq3 = sym.Eq(0.314 * P_N * P_O - P_NO, 0)  # Equation for nitric oxide
# eq4 = sym.Eq((P_N2 + P_N + P_O2 + P_O + P_NO), 0.432992)  # Dalton's Law of partial pressures
# eq5 = sym.Eq((2 * P_O2 + P_O + P_NO) / (2 * P_N2 + P_N + P_NO), 0.2658)  # Ratio of nitrogen to oxygen
#


import sympy as sym

# Define symbols
P_O2, P_O, P_N2, P_N, P_NO = sym.symbols('P_O2 P_O P_N2 P_N P_NO')

# Equations based on the provided relationships
eq1 = sym.Eq(7.52 * P_O2 + 3.76 * P_O - 2 * P_N2 - P_N + 2.76 * P_NO, 0)
eq2 = sym.Eq( P_O2 + P_O + P_N2 + P_NO, 0.43)
eq3 = sym.Eq(244.64 * P_O2 - P_O**2, 0)
eq4 = sym.Eq(0.0234 * P_N2 - P_N**2, 0)
eq5 = sym.Eq(0.227 * P_N * P_O - P_NO, 0)
#
# Attempt to solve the equations
result = sym.solve([eq1, eq2, eq3, eq4, eq5], (P_O2, P_O, P_N2, P_N, P_NO))
