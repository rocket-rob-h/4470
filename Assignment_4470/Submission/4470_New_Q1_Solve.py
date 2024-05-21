# from scipy.optimize import fsolve

# # Define the system of equations
# def system(vars):
#     P_O2, P_O, P_N2, P_N, P_NO = vars
#     eq1 = 185.6 * P_O2 - P_O**2  # Rearranged from P_O**2 = 185.6 * P_O2
#     eq2 = 0.01895 * P_N2 - P_N**2  # Rearranged from P_N**2 = 0.01895 * P_N2
#     eq3 = 0.314 * P_N * P_O - P_NO  # Rearranged from P_NO = 0.314 * P_N * P_O
#     eq4 = P_N2 + P_N + P_O2 + P_O + P_NO - 0.432992  # Dalton's Law of partial pressures
#     eq5 = (P_O2 + P_O + P_NO) / (P_N2 + P_N + P_NO) - 0.2658  # Ratio of oxygen to nitrogen
#     return [
#         eq1,  # eq1 = 0
#         eq2,  # eq2 = 0
#         eq3,  # eq3 = 0
#         eq4,  # eq4 = 0
#         eq5   # eq5 = 0
#     ]
#
# # Initial guess for the variables
# initial_guess = [1, 1, 1, 1, 1]
#
# # Solve the system of equations using fsolve
# solution = fsolve(system, initial_guess)
#
# # Check if the solutions are within the desired range and output results
# if all(0 < sol <= 1 for sol in solution):
#     P_O2, P_O, P_N2, P_N, P_NO = solution
#     print("Solutions found:")
#     print(f"P_O2 (eq1): {P_O2}")
#     print(f"P_O (eq1): {P_O}")
#     print(f"P_N2 (eq2): {P_N2}")
#     print(f"P_N (eq2): {P_N}")
#     print(f"P_NO (eq3): {P_NO}")
# else:
#     print("No valid solution found within the desired range.")
#
# # Display the solution
# solution

#
# Solver
# import numpy as np
# from scipy.optimize import fsolve
#
# # KN is 0.0189519 atm
# # KO is 185.6401803 atm
# # KNO is 0.3142388 atm
#
# def Partial_Pressures(P):
#
#     return [7.52 * P[0] + 3.76 * P[1] - 2*P[2] - P[3] + 2.76*P[4],
#             P[0] + P[1] + P[2] + P[4] - 0.432992,
#             185.6401803 * P[0] - P[1]**2,
#             0.0189519 * P[2] - P[3]**2,
#             0.3142388 * P[3] * P[1] - P[4]]
#
# Pressures = fsolve(Partial_Pressures,[1,1,1,1,1])
#
# P_O2 = Pressures [0]
# P_O = Pressures [1]
# P_N2 = Pressures [2]
# P_N = Pressures [3]
# P_NO = Pressures [4]
#
# print(f'P_O2 = {P_O2}')
# print(f'P_O = {P_O}')
# print(f'P_N2 = {P_N2}')
# print(f'P_N = {P_N}')
# print(f'P_NO = {P_NO}')

# Solver
import numpy as np
from scipy.optimize import fsolve

# Equilibrium constants (in atm)
KN = 0.0189519  # for N2 dissociation
KO = 220.6401803  # for O2 dissociation
KNO = 0.3142388  # for NO formation

# Define the system of equations
def Partial_Pressures(P):
    P_O2, P_O, P_N2, P_N, P_NO = P
    return [
        7.52 * P_O2 + 3.76 * P_O - 2 * P_N2 - P_N + 2.76 * P_NO,
        P_O2 + P_O + P_N2 + P_NO - 0.432992,
        KO * P_O2 - P_O**2,
        KN * P_N2 - P_N**2,
        KNO * P_N * P_O - P_NO
    ]

# Initial guess for the pressures
initial_guess = [1, 1, 1, 1, 1]

# Solve the system of equations
Pressures = fsolve(Partial_Pressures, initial_guess)

# Extract the partial pressures
P_O2, P_O, P_N2, P_N, P_NO = Pressures

# Print the results
print(f'P_O2 = {P_O2}')
print(f'P_O = {P_O}')
print(f'P_N2 = {P_N2}')
print(f'P_N = {P_N}')
print(f'P_NO = {P_NO}')

change to mole fractions and compare mole fractions with the ones below

# Target mole fractions
mole_fractions = {
    "N": 0.16241,
    "NO": 0.0061796,
    "N2": 0.51083,
    "O": 0.31302,
    "O2": 0.00017813
}



# old




# # Given partial pressures
# P_O2 = 0.000137
# P_O = 0.159587
# P_N2 = 0.269702
# P_N = 0.071490
# P_NO = 0.003585
#
# # Total pressure
# P_total = P_O2 + P_O + P_N2 + P_NO
#
# # Calculate mole fractions
# X_O2 = P_O2 / P_total
# X_O = P_O / P_total
# X_N2 = P_N2 / P_total
# X_N = P_N / P_total
# X_NO = P_NO / P_total
#
# # Print the mole fractions
# print(f'X_O2 = {X_O2}')
# print(f'X_O = {X_O}')
# print(f'X_N2 = {X_N2}')
# print(f'X_N = {X_N}')
# print(f'X_NO = {X_NO}')
