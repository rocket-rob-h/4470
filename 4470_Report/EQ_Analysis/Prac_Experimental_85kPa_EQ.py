import numpy as np
from scipy.optimize import fsolve, newton
from Fluids_Library import *

#use measured shock speed, no calibration constant.

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
P_1 = 85e3  # 15kpa
cp_1 = 1005
T_1 = 300  # K

gamma_2 = gamma_1
a_1 = sound_speed(gamma_1, R_1, T_1 )
a_4 = sound_speed(gamma_4, R_4, T_4)

probe_distance = 0.217 #m
probe_time1 = 0.01999 #s graph
probe_time2 = 0.020222 #s graph

U_s = (probe_distance/ (probe_time2-probe_time1))
print(f"Shock Speed U_s 85 kPa {U_s:.4f} m/s ")

M_s = Mach(U_s,a_1)
print(f"M_s Shock mach is {M_s:.4f}")

P2_P1 = shock_P2_P1(M_s, gamma_1)
print(f"The new {P2_P1} P2_P1 is based on M_s")

P_2 = Top_ratio(P_1,P2_P1)
print(f"P_2 from the M_s is {P_2/1000:.4f} kPa")

print()
print( " - Correct above this line - "*3)
print()


print( " - EQ Forced State 5 Values - "*3)

P_5_EQ = 3.289e6 #Mpa
gamma_5 = 1.3275 #Ratio
rho_5_EQ = 10.1 #m^3/kg
T_5_EQ = 1120 #K

# y5 = 1.3275
# P5 = 3.289 mpa
#T5 = 1120
# rho 5 = 10.1

Me_EQ = mach_A_Astar(A_Astar, gamma_5, 3)
print(f"Me_EQ from A_Astar at Nozzle exit EQ is {Me_EQ:.3f}")

T5_EQ_Te_EQ = stagnation_temperature_ratio(gamma_5, Me_EQ)
print(f"T5_EQ_Te_EQ ratio from Me_EQ is {T5_EQ_Te_EQ:.3f} ")

Te_EQ = Bottom_ratio(T_5_EQ,T5_EQ_Te_EQ)
print(f"Te_EQ Exit Temp at EQ is {Te_EQ:.3f} K")

P5_EQ_Pe_EQ = stagnation_pressure_ratio(gamma_5, Me_EQ)
print(f"P5_EQ_Pe_EQ: P5 on Pe ratio at EQ is {P5_EQ_Pe_EQ:.3f}")

Pe_EQ = Bottom_ratio(P_5_EQ, P5_EQ_Pe_EQ)
print(f"Exit Pressure Pe_EQ is {Pe_EQ:.3f} Pa")

rho5_rho_e_EQ = shock_density_ratio(Me_EQ, gamma_5)
print(f"rho5_rho_e_EQ: Density ratio p5/pe at EQ is {rho5_rho_e_EQ:.3f}")

rho_e_EQ = Bottom_ratio(rho_5_EQ,rho5_rho_e_EQ)
print(f"rho_e_EQ: Density at Exit EQ is {rho_e_EQ:.7f}")


print()
print( " - Plate Properties from EQ Forced State 5 Values - "*3)


theta_deg = np.radians(30)
beta_weak, beta_strong = solve_beta(Me_EQ, theta_deg)
beta_deg = np.degrees(beta_weak)
print(f"Weak shock beta: {beta_deg:.3f} degrees")
# print(f"Strong shock beta: {np.degrees(beta_strong):.3f} degrees")

beta_rad = beta_weak

M2_obl_15 = M2_obl(Me_EQ, beta_deg, theta_deg, g=1.4)
print(f"Mn_2 oblique {M2_obl_15:.3f}")

rho2_rho1_obl_15 = rho2_rho1_obl(Me_EQ, beta_deg, g=1.4)
print(f"rho2_rho1 oblique ratio is {rho2_rho1_obl_15:.4f}")

#density 1
rho1_obl = Ideal_gas_law_rho(Pe_EQ,R_1, Te_EQ)
print(f"Density rho_1_nozzle, pre oblique is {rho1_obl:.5f} m^3/kg")

rho_2_obl = Top_ratio(rho1_obl, rho2_rho1_obl_15)
print(f"rho_2_oblique is {rho_2_obl:.5f} m^3/kg")

p2_p1_obl_15 = p2_p1_obl(Me_EQ, beta_rad, g=1.4)
print(f"P2_P1 across oblique shock is {p2_p1_obl_15:.3}")

P2 = Pe_EQ * p2_p1_obl_15
print(f"Plate pressure P2 is {P2:.4f} Pa")


T_2_T1_obl_15 = T2_T1_obl(Me_EQ, beta_rad, g=1.4)
print(f"The T2_T1_oblique ratio is {T_2_T1_obl_15:.4f}")

T_2_obl = Top_ratio(Te_EQ,T_2_T1_obl_15)
print(f"The T_2_oblique from exit Te_EQ is {T_2_obl:.4f} K")

print()

print(f"Beta weak is {beta_weak:.3f} radians")
print(f"Beta weak in degrees is {np.degrees(beta_weak):.4f}")

print()
print("Finished Analysis")



print()
print( " - Don't use this data for experimental - "*3)
print()





#
#
# T2_T1 = T2_T1_ratio(M_s, gamma_1)
# print(f"T2_T1 ratio is {T2_T1:.4f}")
#
# T_2 = Top_ratio(T_1, T2_T1)
# print(f"T_2 is {T_2:.4f} K")
#
# U_2 = Post_shock_velocity(a_1, gamma_1, P_2, P_1)
# print(f"U_2 velocity is {U_2:.4f} m/s")
#
# a_2 = sound_speed(gamma_1,R_1,T_2)
# print(f"a_2 sound speed is {a_2:.4f} m/s")
#
#
# exit_r = (83.9 / 2)/1000
# throat_r = (7 / 2 )/ 1000
#
# A_Astar = (np.pi * (exit_r ** 2))/ ((np.pi * (throat_r ** 2)))
# print(f"Exit area ratio is {A_Astar:.4f}")
#
#
# M_2 = U_2 / a_2
# print(f"Mach number 2 is {M_2:.4f}")
#
# M_r = Mach_reflected(M_2, gamma_1)
# print(f"Reflected shock mach number M_r is {M_r:.3f}")
#
# T5_T2_Mach = T5_T2wopress(gamma_2, M_r)
# print(f"T5_T2 ratio is {T5_T2_Mach:.4f}")
#
# T_5 = Top_ratio(T5_T2_Mach, T_2)
# print(f"T5 stagnation temperature is {T_5:.4f} K")
#
# T5_T1 = T5_T1(T_5, T_1)
# print(f"To Check Theoretical T5_T1 for  graph is {T5_T1:.4f}")
#
# P5_P2 = shock_P2_P1(M_r, gamma_1)
# print(f"P5_P1 ratio using ideal reflected is {P5_P2:.4f}")
#
# P_5 = Top_ratio(P5_P2, P_2)
#
# print(f"Theoretical P5 stagnation is {P_5 /1000:.4f} kPa")
#
# P5_P1 = P_5 / P_1
# print(f"To check the P5_P1 chart {P5_P1:.4f}")
#
# # EQN Analytical
# U_3 = U_2 #
#
# a_3 = speed_sound3(U_4i, U_3, gamma_4, a_4)
# print(f"Sound speed a_3 is {a_3:.4f}")
#
# def EQN(U_2, a_3, gamma_4, P5_P2):
#     gamma = gamma_4 + 1
#     gamma_ = gamma_4 - 1
#     numer = ((gamma/gamma_) - 1) * (P5_P2 - 1)
#     denom = (1 + (gamma / gamma_)) * (1 + ( gamma / gamma_) * P5_P2)
#     return U_2 - a_3 * numer / np.sqrt(denom)
#
# Tailoring15kpa = EQN(U_2, a_3, gamma_4, P5_P2)
#
# print(f"EQN of  system is {Tailoring15kpa:.4f}")

# y5 = 1.3275
# P5 = 3.289 mpa
#T5 = 1120
# rho 5 = 10.1
