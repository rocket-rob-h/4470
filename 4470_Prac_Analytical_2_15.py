#Preamble

import numpy as np
from scipy.optimize import fsolve, newton
from Fluids_Library import *


#Initial Conditions
#Calculate Shock Pressures

# Driver Helium properties
gamma_4 = 1.667
R_4 = 2077  # J/kg.K
P_4 = 3e6  # 3Mpa
cp_4 = 5192.6  # J/kg.k
T_4 = 300  # K
L_4 = 0.23  # m
U_4i = 0

# Print statements for Driver Helium
print(f"Gamma for Helium (γ₄): {gamma_4} (dimensionless)")
print(f"Specific Gas Constant for Helium (R₄): {R_4} J/kg.K")
print(f"Pressure for Helium (P₄): {P_4/1e6} MPa")
print(f"Specific Heat Capacity at Constant Pressure for Helium (cp₄): {cp_4} J/kg.K")
print(f"Temperature for Helium (T₄): {T_4} K")
print(f"Length for Helium (L₄): {L_4} m")

# Test Gas Air properties
gamma_1 = 1.4
R_1 = 287  # J/kg.K
P_1 = 15e3  # 15kpa
cp_1 = 1005
T_1 = 300  # K

# Print statements for Test Gas Air
print(f"Gamma for Air (γ₁): {gamma_1} (dimensionless)")
print(f"Specific Gas Constant for Air (R₁): {R_1} J/kg.K")
print(f"Pressure for Air (P₁): {P_1/1000} kPa")
print(f"Specific Heat Capacity at Constant Pressure for Air (cp₁): {cp_1} J/kg.K")
print(f"Temperature for Air (T₁): {T_1} K")



gamma_2 = gamma_1

def sound_speed(gamma,R,Temperature):
    return np.sqrt(gamma * R * Temperature)

a_1 = sound_speed(gamma_1, R_1, T_1 )

a_4 = sound_speed(gamma_4, R_4, T_4)

P4_P1_known = P_4/P_1 # The given P4/P1 ratio

# Function to solve for P2/P1 using the equation provided
def solve_for_P2_P1(x):
    P2_P1 = x[0]
    term_inside_brackets = (gamma_4 - 1) * (a_1 / a_4) * (P2_P1 - 1) / np.sqrt(2 * gamma_1 * (2 * gamma_1 + (gamma_1 + 1) * (P2_P1 - 1)))
    return [P4_P1_known - P2_P1 * (1 - term_inside_brackets) ** (-2 * gamma_4 / (gamma_4 - 1))]

# Initial guess for P2_P1
P2_P1_initial_guess = [20.0]  # This is a guess and may need to be adjusted
# Use fsolve to find the root
P2_P1_solution = fsolve(solve_for_P2_P1, P2_P1_initial_guess)

print(f"The calculated P2/P1 ratio is: {P2_P1_solution[0]:.4f}")
# Rearrange the ratio to find Pressure of State 2 of accelerated test gas
P_2 = P_1 * P2_P1_solution[0]

#Now with this Pressure ratio, find the resulting Mach number.
#Compressible flow calculator says
print("---" * 5)
print("Data from Compressible flow calculator at P2/P1 = 19.896")

# OUTPUTS:
M_s =  4.14696203       # Upstream Mach number
Mb_2 =  0.43120352       # Downstream Mach number
Po2_Po1 =  0.12263917 # Total pressure ratio across the shock
P1_Po2 =  0.04422986  # Pressure ratio of state 1 to total pressure behind shock
rho2_rho1 =  4.64848391  # Density ratio across the shock
T2_T1 =  4.28028654    # Temperature ratio across the shock

print(f"The M_s is {M_s:.3f}")
print(f"The Mb_2 is { Mb_2:.3f}")
print(f"The rho2_rho1 is {rho2_rho1:.4f}")
print(f"The Po2_Po1 is {Po2_Po1:.4f}")
print(f"The T2_T1 is {T2_T1:.4f}")
print(f"The P1_Po2 is {P1_Po2:.4f}")


#The temperature of State 2 T2
T_2 = T_1 * T2_T1






def Ideal_gas_law_P(rho, R, Temperature):
    return rho * R * Temperature

def Ideal_gas_law_rho(Pressure,R, Temperature):
    return Pressure / (R * Temperature)

#Density initial at 15kpa Air
rho1 = Ideal_gas_law_rho(P_1, R_1, T_1)

#Density initial at 3MPA helium
rho4 = Ideal_gas_law_rho(P_4, R_4, T_4)


# Have mach number and sound speed, now can find velocity
# Remember reference frame issue

def Mach(velocity, sound_speed):
    return velocity / sound_speed

def Velocity(Mach, sound_speed):
    return Mach * sound_speed

#Shock Mach
def Shock_Mach(sound_speed, gamma, P2, P1):
    return sound_speed * np.sqrt ((gamma + 1) / (2 * gamma) * ((P2/P1)-1) + 1 )

U_s = Shock_Mach(a_1,gamma_1, P_2, P_1)


#U2 post shock speed from post shock mach number 0.7? Direct solve

# The velocity behind the shock wave.

def Post_shock_velocity(sound_speed, gamma, P2, P1):
    numerator = (2 * gamma) / (gamma + 1)
    denominator = (P2 / P1) + ((gamma - 1) / (gamma + 1))
    return (sound_speed / gamma) * ((P2 / P1) - 1) * np.sqrt(numerator / denominator)

U_2 = Post_shock_velocity(a_1, gamma_1, P_2, P_1)

a_2 = sound_speed(gamma_1,R_1,T_2)

#Checking with quick reference frame check with post shock mach number from Comp flow calc
U2 = U_s - Mb_2 * a_2

rho2 = Ideal_gas_law_rho(P_2,R_1,T_2)

# Wave Processing Time

#Reference frame for velocity changes.


def stagnation_pressure_ratio(gamma, Mach):
    return (1 + ((gamma - 1) / 2) * Mach ** 2) ** (-gamma / (gamma - 1))

def stagnation_temperature_ratio(gamma, Mach):
    return (1 + ((gamma - 1) / 2) * Mach**2) ** -1

def Top_ratio(Denominator, Ratio):
    return Denominator * Ratio

def Bottom_ratio(Numerator, Ratio):
    return Numerator / Ratio

exit_r = (83.9 / 2)/1000
throat_r = (7 / 2 )/ 1000

A_Astar = (np.pi * (exit_r ** 2))/ ((np.pi * (throat_r ** 2)))

print("--State 4--" * 7)
print(f"The sound speed for State 4 a_4 is {a_4:.2f} m/s")
print(f"The rho4 Density of Helium rho4 is {rho4:.5f} kg/m^3")

print("--State 1--" * 7)
print(f"The sound speed for State 1 a_1 is {a_1:.2f} m/s")
print(f"The rho1 Density of Air rho1 is {rho1:.5f} kg/m^3")
print(f"The shock velocity into the State 1 test gas U_s is {U_s:.2f} m/s" )


print("--State 2--" * 7)
print(f"The calculated P2/P1 ratio is: {P2_P1_solution[0]:.4f}")
print(f"The P_2 Pressure is {P_2/1000:.4f} kPa")
print(f"The T_2 Gas Temperature is {T_2:.2f} K")
print(f"The rho2 Density of Air is {rho2:.5f} kg/m^3")
print(f"The velocity behind the shock wave. Checking the reference frame velocity U2 - main way {U2:.2f} m/s")
print(f"The State 2 Post shock The velocity behind the shock wave short cut is  U_2 {U_2:.2f} m/s")
print(f"The sound speed for State 2 a_2 is {a_2:.2f} m/s")

print("-- Intermediate Steps --" *3)

#stagnation points
M_s_2 = Mach(U_s, a_1)
print(f"The mach number of the accelerated test gas is {M_s_2:.2f} ")

print("---" * 5)

#Post shock flow properties

#Tutorial 6

print(f"The exit area ratio is {A_Astar:.2f}")
print("Stagnation Regions")

#Stagnation regions
#P_e = design pressure of the nozzle?
# P_e = 280 # Pa off dial/given.

# Now Stagnation Pressure from Ratio
# Area ratio (A_A*): 143.660005

# Mach number (M): 7.52120812
# Mach angle (mu): 7.64052063 degrees
# Prandtl-Meyer angle (nu): 93.5377253 degrees
# Pressure ratio across shock (P_P0): 0.00015262
# Density ratio across shock (rho_rho0): 0.00187943
# Stagnation pressure ratio across shock (P*_P0*): 0.00028891
# Stagnation density ratio across shock (rho*_rho0*): 0.00296469
# Temperature ratio across shock (T_T0): 0.08121026
# Stagnation temperature ratio across shock (T*_T0*): 0.09745231


print("--State 5--" * 7)

P2_P1 = P2_P1_solution

def find_mach_from_pressure_ratio(pressure_ratio, gamma, initial_guess):
    """
    Find the upstream Mach number (M1) for a given pressure ratio (p2/p1) in a normal shock.

    :param pressure_ratio: The pressure ratio across the shock (p2/p1).
    :param gamma: Specific heat ratio (γ), for air it is typically 1.4.
    :param initial_guess: Initial guess for the Mach number.
    :return: The upstream Mach number (M1).
    """

    def f(M1):
        # Implicit equation for the upstream Mach number based on pressure ratio
        return 1 + (2 * gamma / (gamma + 1)) * (M1 ** 2 - 1) - pressure_ratio

    # Use the Newton-Raphson method to solve for M1
    M1 = newton(f, initial_guess)
    return M1

M_1_from_P2_P1 = find_mach_from_pressure_ratio(P2_P1, gamma_1, 4)
print(f" The mach number at this pressure ratio is {M_1_from_P2_P1}")

M_2 = U_2 / a_2
print(f"The Mach number 2 is {M_2:.2f}")

M_r = Mach_reflected(M_2, gamma_1)
print(f"The reflected shock mach number M_r is {M_r:.3f}")

T5_T2_Mach = T5_T2wopress(gamma_2, M_r)
print(f"The T5_T2 ratio is {T5_T2_Mach:.4f}")

T_5 = Top_ratio(T5_T2_Mach, T_2)
print(f"The T5 stagnation temperature is {T_5:.2f} K")

T5_T1 = T5_T1(T_5, T_1)
print(f"To Check the theoretical T5_T1 for the graph is {T5_T1:.2f}")

P5_P2 = shock_P2_P1(M_r, gamma_1)
print(f"FTP5_P1 ratio using ideal reflected is {P5_P2:.2f}")

P_5 = Top_ratio(P5_P2, P_2)

print(f"Theroretical P5 stagnation is {P_5 /1000:.2f} kPa")

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

print(f"The EQN of the system is {Tailoring15kpa:.2f}")

