import numpy as np
from scipy.optimize import fsolve, newton
from Fluids_Library import *
from math import *

# Driver Helium properties
gamma_4 = 1.667
R_4 = 2077  # J/kg.K
P_4 = 3e6  # 3Mpa
cp_4 = 5192.6  # J/kg.k
T_4 = 300  # K
L_4 = 0.23  # m
U_4i = 0

# Test Gas Air properties
gamma_1 = 1.4
R_1 = 287  # J/kg.K
P_1 = 15e3  # 15kpa
cp_1 = 1005
T_1 = 300  # K

gamma_2 = gamma_1
a_1 = sound_speed(gamma_1, R_1, T_1 )
a_4 = sound_speed(gamma_4, R_4, T_4)



def solve_for_P2_P1(x):
    P2_P1 = x[0]
    term_inside_brackets = (gamma_4 - 1) * (a_1 / a_4) * (P2_P1 - 1) / np.sqrt(2 * gamma_1 * (2 * gamma_1 + (gamma_1 + 1) * (P2_P1 - 1)))
    return [(P_4 / P_1 )- P2_P1 * (1 - term_inside_brackets) ** (-2 * gamma_4 / (gamma_4 - 1))]

# Initial guess for P2_P1
P2_P1_initial_guess = [20.0]  # This is a guess and may need to be adjusted
# Use fsolve to find  root
P2_P1_solution = fsolve(solve_for_P2_P1, P2_P1_initial_guess)


print(f"Calculated P2/P1 ratio is: {P2_P1_solution[0]:.4f}")
# Rearrange  ratio to find Pressure of State 2 of accelerated test gas
P_2 = P_1 * P2_P1_solution[0]
P2_P1 = P2_P1_solution
print(f"P_2 is {P_2/1000:.2f} kPa")

M_1_value = find_mach_from_pressure_ratio(P2_P1, gamma_1, 2)
M_1 = M_1_value[0]
print(f"Mach number of P2/P1 is {M_1:.2f}")

U_s = Shock_velocity(a_1,gamma_1, P_2, P_1)
print(f"Shock Speed U_s {U_s:.2f} m/s ")

T2_T1 = T2_T1_ratio(M_1, gamma_1)
print(f"T2_T1 ratio is {T2_T1:.2f}")

T_2 = Top_ratio(T_1, T2_T1)
print(f"T_2 is {T_2:.2f} K")

U_2 = Post_shock_velocity(a_1, gamma_1, P_2, P_1)
print(f"U_2 velocity is {U_2:.2f} m/s")

a_2 = sound_speed(gamma_1,R_1,T_2)
print(f"a_2 sound speed is {a_2:.2f} m/s")


exit_r = (83.9 / 2)/1000
throat_r = (7 / 2 )/ 1000

A_Astar = (np.pi * (exit_r ** 2))/ ((np.pi * (throat_r ** 2)))
print(f"Exit area ratio is {A_Astar:.2f}")


M_2 = U_2 / a_2
print(f"Mach number 2 is {M_2:.2f}")

M_r = Mach_reflected(M_2, gamma_1)
print(f"Reflected shock mach number M_r is {M_r:.3f}")

T5_T2_Mach = T5_T2wopress(gamma_2, M_r)
print(f"T5_T2 ratio is {T5_T2_Mach:.4f}")

T_5 = Top_ratio(T5_T2_Mach, T_2)
print(f"T5 stagnation temperature is {T_5:.2f} K")

T5_T1 = T5_T1(T_5, T_1)
print(f"To Check Theoretical T5_T1 for  graph is {T5_T1:.2f}")

P5_P2 = shock_P2_P1(M_r, gamma_1)
print(f"P5_P1 ratio using ideal reflected is {P5_P2:.2f}")

P_5 = Top_ratio(P5_P2, P_2)

print(f"Theoretical P5 stagnation is {P_5 /1000:.2f} kPa")

P5_P1 = P_5 / P_1
print(f"To check the P5_P1 chart {P5_P1:.2f}")

# EQN Analytical
U_3 = U_2 #

a_3 = speed_sound3(U_4i, U_3, gamma_4, a_4)
print(f"Sound speed a_3 is {a_3:.2f}")

def EQN(U_2, a_3, gamma_4, P5_P2):
    gamma = gamma_4 + 1
    gamma_ = gamma_4 - 1
    numer = ((gamma/gamma_) - 1) * (P5_P2 - 1)
    denom = (1 + (gamma / gamma_)) * (1 + ( gamma / gamma_) * P5_P2)
    return U_2 - a_3 * numer / np.sqrt(denom)

Tailoring15kpa = EQN(U_2, a_3, gamma_4, P5_P2)

print(f"EQN of 15 kPa system is {Tailoring15kpa:.2f}")

print()
print("--Exit Conditions--" * 3)
print()

# Exit conditions
mach_from_A_Astar = mach_A_Astar(A_Astar, gamma_2, 10)
print(f"Mach from A_Astar is {mach_from_A_Astar:.3f}")

P5_Pe = stagnation_pressure_ratio(gamma_1, mach_from_A_Astar)
print(f"P5_Pe ratio is {P5_Pe:.10f} ratio ")

P_e = Bottom_ratio(P_5,P5_Pe)
print(f"Exit Pressure P_e is {P_e:.2f} Pa")

T5_Te = stagnation_temperature_ratio(gamma_1, mach_from_A_Astar)
print(f"Exit Temperature T5_Te ratio is {T5_Te:.3f} ratio")

T_e = Bottom_ratio(T_5,T5_Te)
print(f"Exit Temperature T_e is {T_e:.2f} K")


print()
print("--Plate Work --" * 3)
print()


theta_deg = np.radians(30)
beta_weak, beta_strong = solve_beta(mach_from_A_Astar, theta_deg)
beta_deg = np.degrees(beta_weak)
print(f"Weak shock beta: {beta_deg:.3f} degrees")
# print(f"Strong shock beta: {np.degrees(beta_strong):.3f} degrees")

beta_rad = beta_weak

M2_obl_15 = M2_obl(mach_from_A_Astar, beta_deg, theta_deg, g=1.4)
print(f"Mn_2 oblique {M2_obl_15:.3f}")

rho2_rho1_obl_15 = rho2_rho1_obl(mach_from_A_Astar, beta_deg, g=1.4)
print(f"rho2_rho1 oblique ratio is {rho2_rho1_obl_15:.2f}")

#density 1
rho1_obl = Ideal_gas_law_rho(P_e,R_1, T_e)
print(f"Density rho_1_nozzle, pre oblique is {rho1_obl:.5f} m^3/kg")

rho_2_obl = Top_ratio(rho1_obl, rho2_rho1_obl_15)
print(f"rho_2_oblique is {rho_2_obl:.5f} m^3/kg")

p2_p1_obl_15 = p2_p1_obl(mach_from_A_Astar, beta_rad, g=1.4)
print(f"P2_P1 across oblique shock is {p2_p1_obl_15:.3}")

P2 = P_e * p2_p1_obl_15
print(f"Plate pressure P2 is {P2:.2f} Pa")


T_2_T1_obl_15 = T2_T1_obl(mach_from_A_Astar, beta_rad, g=1.4)
print(f"The T2_T1_oblique ratio is {T_2_T1_obl_15:.2f}")

T_2_obl = Top_ratio(T_e,T_2_T1_obl_15)
print(f"The T_2_oblique from exit T_e is {T_2_obl:.2f} K")

print()

print(f"Beta weak is {beta_weak:.3f} radians")
print(f"Beta weak in degrees is {np.degrees(beta_weak):.2f}")

print()
print("Finished Analysis")