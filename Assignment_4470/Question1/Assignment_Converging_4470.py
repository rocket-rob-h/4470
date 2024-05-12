import sympy as sym
import numpy as np
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

# import sympy as sym
#
# # Define symbols
# P_O2, P_O, P_N2, P_N, P_NO = sym.symbols('P_O2 P_O P_N2 P_N P_NO')
#
# # Equations based on the provided relationships
# eq1 = sym.Eq(7.52 * P_O2 + 3.76 * P_O - 2 * P_N2 - P_N + 2.76 * P_NO, 0)
# eq2 = sym.Eq(P_O2 + P_O + P_N2 + P_NO, 0.43)
# eq3 = sym.Eq(185.* P_O2 - P_O**2, 0)
# eq4 = sym.Eq(0.0234 * P_N2 - P_N**2, 0)
# eq5 = sym.Eq(0.227 * P_N * P_O - P_NO, 0)
# #
# # Attempt to solve the equations
# result = sym.solve([eq1, eq2, eq3, eq4, eq5], (P_O2, P_O, P_N2, P_N, P_NO))
# print(result)

#
#
# solutions = [
#     (61.3430039219891 + 1.73472347597681e-18*1j, -122.502867229610, 36.049387458986 + 8.67361737988404e-19*1j, -0.918452865715096, 25.5404758486347 + 8.56315551627477e-27*1j),
#     (7.65806823596473e-5, 0.136874753451702 + 1.6940658945086e-21*1j, 0.29563290991512 + 6.7762635780344e-21*1j, -0.0831733736962365, -0.00258424404918103 - 5.4871917742325e-32*1j),
#     (152.745058470859 + 3.46944695195361e-18*1j, -193.306883230554 - 3.46944695195361e-18*1j, 112.043527737579 + 1.73472347597681e-18*1j, 1.61920305985980, -71.0517029778846 - 8.56315551625591e-27*1j),
#     (0.000105253334933463 + 1.65436122510606e-24*1j, 0.160465497407145, 0.266552468224606 + 6.7762635780344e-21*1j, 0.0789767545323039 - 1.6940658945086e-21*1j, 0.00287678103331559 + 9.44461951394779e-32*1j)
# ]
#
# # Filter and print only positive real values
# filtered_solutions = []
# for solution in solutions:
#     filtered_solution = []
#     for value in solution:
#         # Check if the value is positive and real
#         if isinstance(value, complex):
#             if value.imag == 0 and value.real > 0:  # Check for real part positive and imaginary part zero
#                 filtered_solution.append(value.real)
#         elif value > 0:  # Handle non-complex numbers that are positive
#             filtered_solution.append(value)
#     if filtered_solution:
#         filtered_solutions.append(tuple(filtered_solution))
#
# # Output the filtered solutions
# print("Filtered positive real solutions:")
# for sol in filtered_solutions:
#     print(sol)

#
# P_O2, P_O, P_N2, P_N, P_NO = sym.symbols('P_O2 P_O P_N2 P_N P_NO')
#
# # Equations based on the provided relationships
# eq1 = sym.Eq(7.52 * P_O2 + 3.76 * P_O - 2 * P_N2 - P_N + 2.76 * P_NO, 0)
# eq2 = sym.Eq(2 * P_O2 + P_O + 2 * P_N2 + P_NO, 0.43299)
# eq3 = sym.Eq(185.6 * P_O2 - P_O**2, 0)
# eq4 = sym.Eq(0.01895 * P_N2 - P_N**2, 0)
# eq5 = sym.Eq(0.314 * P_N * P_O - P_NO, 0)
#
# # Attempt to solve the equations
# result = sym.solve([eq1, eq2, eq3, eq4, eq5], (P_O2, P_O, P_N2, P_N, P_NO))
# print(result)

PN2 = 0.26970243
PO2 = 0.00013719

PN = np.sqrt(0.01895*PN2)
PO = np.sqrt(185.64018*PO2)
PNO = 0.31424 * np.sqrt((185.64018 * 0.01895) * PO2 * PN2)


print('Partial Pressures')
print(f'PN2 is {PN2:.10f}')
print(f'PO2 is {PO2:.10f}')
print(f'PN is {PN:.10f}')
print(f'PO is {PO:.10f}')
print(f'PNO is {PNO:.10f}')

'''Mole Fractions'''


# Calculate total pressure
P_total = 0.51

# Calculate mole fractions
X_N2 = PN2 / P_total
X_O2 = PO2 / P_total
X_N = PN / P_total
X_O = PO / P_total
X_NO = PNO / P_total

# Print the results
print(f"Mole Fraction of N2: {X_N2:.8f}")
print(f"Mole Fraction of O2: {X_O2:.8f}")
print(f"Mole Fraction of N: {X_N:.8f}")
print(f"Mole Fraction of O: {X_O:.8f}")
print(f"Mole Fraction of NO: {X_NO:.8f}")


# Target mole fractions
mole_fractions = {
    "N": 0.16241,
    "NO": 0.0061796,
    "N2": 0.51083,
    "O": 0.31302,
    "O2": 0.00017813
}

# Calculate percentage differences
diff_N = abs((X_N - mole_fractions["N"]) / mole_fractions["N"]) * 100
diff_NO = abs((X_NO - mole_fractions["NO"]) / mole_fractions["NO"]) * 100
diff_N2 = abs((X_N2 - mole_fractions["N2"]) / mole_fractions["N2"]) * 100
diff_O = abs((X_O - mole_fractions["O"]) / mole_fractions["O"]) * 100
diff_O2 = abs((X_O2 - mole_fractions["O2"]) / mole_fractions["O2"]) * 100

# Print the results
print(f"Percentage Difference for N: {diff_N:.2f}%")
print(f"Percentage Difference for NO: {diff_NO:.2f}%")
print(f"Percentage Difference for N2: {diff_N2:.2f}%")
print(f"Percentage Difference for O: {diff_O:.2f}%")
print(f"Percentage Difference for O2: {diff_O2:.2f}%")