#Preamble

import numpy as np
from scipy.optimize import fsolve


#Initial Conditions
#Calculate Shock Pressures

# Driver Helium properties
gamma_4 = 1.667
R_4 = 2077  # J/kg.K
P_4 = 3e6  # 3Mpa
cp_4 = 5192.6  # J/kg.k
T_4 = 300  # K
L_4 = 0.23  # m
U_4 = 0

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
print(f"The P_2 Pressure is {P_2:.2f} Pa")
print(f"The T_2 Gas Temperature is {T_2:.2f} K")
print(f"The rho2 Density of Air is {rho2:.5f} kg/m^3")
print(f"The velocity behind the shock wave. Checking the reference frame velocity U2 - main way {U2:.2f} m/s")
print(f"The State 2 Post shock The velocity behind the shock wave short cut is  U_2 {U_2:.2f} m/s")
print(f"The sound speed for State 2 a_2 is {a_2:.2f} m/s")

print("-- Intermediate Steps --" *3)

#stagnation points
# M_s_2 = Mach(U_s, a_1)
# print(f"The mach number of the accelerated test gas is {M_s_2:.2f} ")

print("---" * 5)

#Post shock flow properties

#Tutorial 6

print(f"The exit area ratio is {A_Astar:.2f}")
print("Stagnation Regions")

#Stagnation regions
#P_e = design pressure of the nozzle?
P_e = 280 # Pa off dial/given.

# Stagnation Pressure from Mach
Mach_exit = 7.52 # Mach number at 143.66 area ratio ??? slowed by reflection wave?

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

# Find temperature at the 280 pascal test section pressure. 15kpa pressure in test section?

# Altitudes (meters)
H1, H2 = 39600, 40000

# Pressures (millibars)
P1, P2 = 2.9308, 2.7752

# Temperatures (Kelvin)
T1, T2 = 249.93, 251.05

# Densities (kilograms per cubic meter)
rho1, rho2 = 0.02079, 0.01852

# Interpolated pressure (millibars)
P_target = 2.8

# Interpolation for altitude
alt_target = H1 + (H2 - H1) * (P_target - P1) / (P2 - P1)

# Interpolation for temperature
T_e = T1 + (T2 - T1) * (P_target - P1) / (P2 - P1)

# Interpolation for density
rho_target = rho1 + (rho2 - rho1) * (P_target - P1) / (P2 - P1)

# Output the interpolated values
print(f"Design altitude: {alt_target:.2f} meters")
print(f"Exit Temperature / T6: {T_e:.2f} Kelvin")
print(f"Exit Density/ rho6: {rho_target:.6f} kg/m^3")

print("--State 5--" * 7)

StagPRatio5 = stagnation_pressure_ratio(gamma_1, Mach_exit)
print(f"Stagnation Pressure ration at Nozzle Wall is {StagPRatio5:.10f} ")

StagTRatio5 = stagnation_temperature_ratio(gamma_1, Mach_exit)
print(f"Stagnation Temperature ration at Nozzle Wall is {StagTRatio5:.10f} ")

stagnation_pressure_nozzle5 = Bottom_ratio(P_e, StagPRatio5)
print(f"Stagnation Pressure at Nozzle Wall is {stagnation_pressure_nozzle5:.3f} ")

stagnation_temperature_nozzle5 = Bottom_ratio(T_1,StagTRatio5)
print(f"Stagnation Temperature at Nozzle Wall is {stagnation_temperature_nozzle5:.3f} K")
# Oxygen Dissociation temperature?

# These pressures and temperatures are the State 5 variables. Since we have the nozzle info and working back to the stagnation points, thats the state 5 reflected numbers.

# T_5 = stagnation_temperature_nozzle5
# P_5 = stagnation_pressure_nozzle5

# should make a state 6 or exit.


# Now to the primary shock graph since we now have a T5/T1 (inverse) and P5/P1 ratio.
T5_T1 = stagnation_temperature_nozzle5 / T_1
print(f"Stagnation T5/T1 for graph is {T5_T1:.2f} ")

# P5_P1 = 250 # at 298K
# M5 = 5.6 # at 298K


#it says the T1 = 298K though, not Te

#Not sure if this is needed - check the P5 calcs with stagnation or graph.
#a1s =
print("--State 3--" * 7)

# def sound_speed3(u4, a4, gamma4, u3):
#     return ((u4 + 2 * a4) / ((gamma4 - 1) - u3)) ** (gamma4 - 1)/ 2

def speed_sound3(U4, U3, gamma4, a_4):
    return ((U4 - U3) * (gamma4 - 1) / 2) + a_4


U_3 = U_2
a_3 = speed_sound3(U_4, U_3, gamma_4, a_4)
print(f"The sound speed of State 3 going backwards is {a_3:.2f} m/s")

# New stagnation values - stagnation pressure is correct then of 1.832 Mpa at end wall. Now temp must be T_e based.

ProperT5 = stagnation_temperature_ratio(gamma_1,4.147)

# probs dont even need T5, just go onto the taloring part and state 3 gas.

P5_P2 = stagnation_pressure_nozzle5 / P_2
print(f"The P5/P2 ratio for EQN is {P5_P2:.5f}")

def EQN(U_2, a_3, gamma_4, P5_P2):
    gamma = gamma_4 + 1
    gamma_ = gamma_4 - 1
    numer = ((gamma/gamma_) - 1) * (P5_P2 - 1)
    denom = (1 + (gamma / gamma_)) * (1 + ( gamma / gamma_) * P5_P2)
    return U_2 - a_3 * numer / np.sqrt(denom)

Tailoring15kpa = EQN(U_2, a_3, gamma_4, P5_P2)

print(f"The EQN of the system is {Tailoring15kpa}")
# = 254


# def calculate_EQN(u2, a3, gamma4, p5, p2):
#     """
#     Calculates the EQN based on the provided equation.
#
#     Parameters:
#     u2 : float
#         Velocity u2 from the equation.
#     a3 : float
#         Speed of sound a3 from the equation.
#     gamma4 : float
#         Ratio of specific heats gamma4 from the equation.
#     p5 : float
#         Pressure p5 from the equation.
#     p2 : float
#         Pressure p2 from the equation.
#
#     Returns:
#     float
#         The calculated value of EQN.
#     """
#     numerator = ((gamma4 + 1) / (gamma4 - 1) - 1) * ((p5 / p2) - 1)
#     denominator = np.sqrt((1 + (gamma4 + 1) / (gamma4 - 1)) * (1 + (gamma4 + 1) / (gamma4 - 1) * p5 / p2))
#     EQN = u2 - a3 * (numerator / denominator)
#     return EQN
#
# try2 = calculate_EQN(U_2, a_3, gamma_4, stagnation_pressure_nozzle5, P_2)
# print (try2)