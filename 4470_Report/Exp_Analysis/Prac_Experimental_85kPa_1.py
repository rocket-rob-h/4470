import numpy as np
from scipy.optimize import fsolve, newton
from Fluids_Library import *

# Driver Helium properties
gamma_4 = 1.667
R_4 = 2077  # J/kg.K
P_4 = 3e6  # 3Mpa
cp_4 = 5192.6  # J/kg.k
T_4 = 300  # K
U_4i = 0

# Test Gas Air properties
gamma_1 = 1.4
R_1 = 287  # J/kg.K
P_1 = 85e3  # 15kpa
T_1 = 300  # K

gamma_2 = gamma_1
a_1 = sound_speed(gamma_1, R_1, T_1 )
a_4 = sound_speed(gamma_4, R_4, T_4)

probe_distance = 0.217 #m
probe_time1 = 0.01999 #s graph
probe_time2 = 0.020222 #s graph

U_s = (probe_distance/ (probe_time2-probe_time1))
print(f"Shock Speed U_s 85 kPa {U_s:.2f} m/s ")

M_s = Mach(U_s,a_1)
print(f"M_s Shock mach is {M_s:.2f}")

P2_P1 = shock_P2_P1(M_s, gamma_1)
print(f"The new {P2_P1} P2_P1 is based on M_s")

P_2 = Top_ratio(P_1,P2_P1)
print(f"P_2 from the M_s is {P_2/1000:.2f} kPa")

# The Pressure P5 at 85 Kpa from 2nd sensor is 3069.4625816075

print( " - EQ Forced State 5 Values - "*3)

P_5 = 3069.46258e3 #Mpa
gamma_5 = 1.3275 #Ratio
rho_5 = 10.1 #m^3/kg
T_5  = 1120 * (P_5/ 3.289e6) #K
print(f"The scaled T5 is {T_5:.2f} K")

Me = mach_A_Astar(A_Astar, gamma_5, 3)
print(f"Me from A_Astar at Nozzle exit Exp is {Me:.3f}")

T5_Te = stagnation_temperature_ratio(gamma_5, Me)
print(f"T5_Te ratio from Me is {T5_Te:.3f} ")

Te = Bottom_ratio(T_5,T5_Te)
print(f"Te Exit Temp at Exp is {Te:.3f} K")

P5_Pe = stagnation_pressure_ratio(gamma_5, Me)
print(f"P5_Pe: P5 on Pe ratio at Exp is {P5_Pe:.3f}")

Pe = Bottom_ratio(P_5, P5_Pe)
print(f"Exit Pressure Pe is {Pe:.3f} Pa")

rho5_rho_e = shock_density_ratio(Me, gamma_5)
print(f"rho5_rho_e: Density ratio p5/pe at Exp is {rho5_rho_e:.3f}")

rho_e = Bottom_ratio(rho_5,rho5_rho_e)
print(f"rho_e: Density at Exit Exp is {rho_e:.7f}")

print()
print( " - Plate Properties from Exp Forced State 5 Values - "*2)

theta_deg = np.radians(30)
beta_weak, beta_strong = solve_beta(Me, theta_deg)
beta_deg = np.degrees(beta_weak)
print(f"Weak shock beta: {beta_deg:.3f} degrees")
# print(f"Strong shock beta: {np.degrees(beta_strong):.3f} degrees")

beta_rad = beta_weak

M2_obl_15 = M2_obl(Me, beta_deg, theta_deg, g=1.4)
print(f"Mn_2 oblique {M2_obl_15:.3f}")

rho2_rho1_obl_15 = rho2_rho1_obl(Me, beta_deg, g=1.4)
print(f"rho2_rho1 oblique ratio is {rho2_rho1_obl_15:.4f}")

#density 1
rho1_obl = Ideal_gas_law_rho(Pe,R_1, Te)
print(f"Density rho_1_nozzle, pre oblique is {rho1_obl:.5f} m^3/kg")

rho_2_obl = Top_ratio(rho1_obl, rho2_rho1_obl_15)
print(f"rho_2_oblique is {rho_2_obl:.5f} m^3/kg")

p2_p1_obl_15 = p2_p1_obl(Me, beta_rad, g=1.4)
print(f"P2_P1 across oblique shock is {p2_p1_obl_15:.3f}")

P2 = Pe * p2_p1_obl_15
print(f"Plate pressure P2 is {P2:.4f} Pa")

T_2_T1_obl_15 = T2_T1_obl(Me, beta_rad, g=1.4)
print(f"The T2_T1_oblique ratio is {T_2_T1_obl_15:.4f}")

T_2_obl = Top_ratio(Te,T_2_T1_obl_15)
print(f"The T_2_oblique from exit Te is {T_2_obl:.4f} K")

print()

print(f"Beta weak is {beta_weak:.3f} radians")
print(f"Beta weak in degrees is {np.degrees(beta_weak):.4f}")

print()
print("Finished Analysis")


T2_T1 = T2_T1_ratio(M_s, gamma_1)
print(f"T2_T1 ratio is {T2_T1:.2f}")

T_2 = Top_ratio(T_1, T2_T1)
print(f"T_2 is {T_2:.2f} K")

# U_2 = Post_shock_velocity(a_1, gamma_1, P_2, P_1)

U_2 = Post_shock_velocity(a_1, gamma_1, P_2, P_1)
print(f"U_2 velocity is {U_2:.2f} m/s")
# U2 from Measured Us

a_2 = sound_speed(gamma_1,R_1,T_2)
print(f"a_2 sound speed is {a_2:.2f} m/s")

M_2 = U_2 / a_2
print(f"Mach number 2 is {M_2:.2f}")

M_r = Mach_reflected(M_2, gamma_1)
print(f"Reflected shock mach number M_r is {M_r:.3f}")


Exp_P5_P2 = P_5 / P_2
print(f"Experimental P5 from graph and P2 from Ms ratio is {Exp_P5_P2:2f}")

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

Tailoring85kpa = EQN(U_2, a_3, gamma_4, Exp_P5_P2)

print(f"EQN of  system is {Tailoring85kpa:.2f}")

