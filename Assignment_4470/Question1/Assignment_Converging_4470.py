# import sympy as sym
#
# # Define symbols
# P_O2, P_O, P_N2, P_N, P_NO = sym.symbols('P_O2 P_O P_N2 P_N P_NO', real=True, positive=True)
# P2_2_atm = 0.432992  # Assume this constant, previously used
#
# # Define equations
# eq1 = sym.Eq(185.6 * P_O2 - P_O**2, 0)  # Quadratic equation for oxygen
# eq2 = sym.Eq(0.01895 * P_N2 - P_N**2, 0)  # Quadratic equation for nitrogen
# eq3 = sym.Eq(0.314 * P_N * P_O - P_NO, 0)  # Equation for nitric oxide
# eq4 = sym.Eq((P_N2 + P_N + P_O2 + P_O + P_NO), 0.432992)  # Dalton's Law of partial pressures
# eq5 = sym.Eq((P_O2 + P_O + P_NO) / (P_N2 + P_N + P_NO), 0.2658)  # Ratio of oxygen to nitrogen
#
# # Attempt to solve the equations symbolically
# result = sym.solve([eq1, eq2, eq3, eq4, eq5], (P_O2, P_O, P_N2, P_N, P_NO), dict=True)
#
# # Check if solutions are within the desired range and output results
# if result:
#     print("Solutions found:")
#     valid_solutions = []
#     for sol in result:
#         # Check if all components of the solution are between 0 and 1
#         if all(0 < val <= 1 for val in sol.values()):
#             valid_solutions.append(sol)
#     if valid_solutions:
#         for sol in valid_solutions:
#             print(sol)
#     else:
#         print("No solutions within the range (0, 1].")
# else:
#     print("No solutions found or solve failed to converge.")

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
print(result)