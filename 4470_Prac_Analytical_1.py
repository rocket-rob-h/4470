#Preamble

import numpy as np
from scipy.optimize import fsolve


#Initial Conditions
#Calculate Shock Pressures

#Driver Helium
gamma_4 = 1.667
R_4 = 2077 #J/kg.K
P_4 = 3e6 #3Mpa
cp_4 = 5192.6 #J/kg.k
T_4 = 300 #K
L_4 = 0.23 #m



#Test Gas Air
gamma_1 = 1.4
R_1 = 287 #J/kg.K
P_1 = 15e3 #15kpa
cp_1 = 1005
T_1 = 300 #K
# L_1 =



def sound_speed(gamma,R,Temperature):
    return np.sqrt(gamma * R * Temperature)

a_1 = sound_speed(gamma_1, R_1, T_1 )
a_4 = sound_speed(gamma_4, R_4, T_4)

P4_P1_given = P_4/P_1 # The given P4/P1 ratio

# Function to be zeroed
def equation(P2_P1):
    term_inside_brackets = (gamma_4 - 1) / (a_1 / a_4) * (P2_P1 - 1) / np.sqrt(2 * gamma_1 * (2 * gamma_1 + (gamma_1 + 1) * (P2_P1 - 1)))
    return P4_P1_given - P2_P1 * (1 - term_inside_brackets) ** (-2 * gamma_4 / (gamma_4 - 1))

# Initial guess for P2_P1
P2_P1_initial_guess = 2.0  # This value may need to be adjusted based on the problem specifics

# Solve for P2/P1
P2_P1_solution = fsolve(equation, P2_P1_initial_guess)



P_2 = P_1 * P2_P1_solution[0]


#Now with this Pressure ratio, find the resulting Mach number.
#Compressible flow calculator says
print("---" * 5)
M_1 = 1.42815565
Mb_1 = 0.72814580
rho2_rho1 = 1.73841152
Po2_Po1 = 0.9508851
T2_T1 = 1.27294370
P1_Po2 = 0.31757134

print(f"The M_1 is {M_1:.2f}")
print(f"The Mb_1 is { Mb_1:.2f}")
print(f"The rho2_rho1 is {rho2_rho1:.2f}")
print(f"The Po2_Po1 is {Po2_Po1:.2f}")
print(f"The T2_T1 is {T2_T1:.2f}")
print(f"The P1_Po2 is {P1_Po2:.2f}")


T_2 = T_1 * T2_T1



#Density initial at 3MPA helium
#Density initial at 15kpa Air

def Ideal_gas_law_P(rho, R, Temperature):
    return rho * R * Temperature

def Ideal_gas_law_rho(Pressure,R, Temperature):
    return Pressure / (R * Temperature)

rho1 = Ideal_gas_law_rho(P_1, R_1, T_1)

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

def Post_shock_velocity(sound_speed, gamma, P2, P1):
    numerator = (2 * gamma) / (gamma + 1)
    denominator = (P2 / P1) + ((gamma - 1) / (gamma + 1))
    return (sound_speed / gamma) * ((P2 / P1) - 1) * np.sqrt(numerator / denominator)

U_2 = Post_shock_velocity(a_1, gamma_1, P_2, P_1)

a_2 = sound_speed(gamma_1,R_1,T_2)

#Checking with quick reference frame check with post shock mach number from Comp flow calc
U2 = U_s - Mb_1 * a_2

rho2 = Ideal_gas_law_rho(P_2,R_1,T_2)

# Wave Processing Time

#Reference frame for velocity changes.



print("--State 4--" * 7)
print(f"The sound speed for State 4 is {a_4:.2f} m/s")
print(f"The rho4 Density of Helium is {rho4:.5f} kg/m^3")

print("--State 1--" * 7)
print(f"The sound speed for State 1 is {a_1:.2f} m/s")
print(f"The rho1 Density of Air is {rho1:.5f} kg/m^3")
print(f"The shock velocity into the State 1 test gas is {U_s:.2f} m/s" )


print("--State 2--" * 7)
print(f"The calculated P2/P1 ratio is: {P2_P1_solution[0]:.4f}")
print(f"The P2 Pressure is {P_2:.2f} Pa")
print(f"The T2 Gas Temperature is {T_2:.2f} K")
print(f"The rho2 Density of Air is {rho2:.5f} kg/m^3")
print(f"Checking the reference frame velocity {U2:.2f} m/s")
print(f"The State 2 Post shock flow velocity is {U_2:.2f} m/s")
print(f"The sound speed for State 2 is {a_2:.2f} m/s")



print("--State 3--" * 7)









print("--State 5--" * 7)

#stagnation points
# M_1_2 = Mach(U_s, a_1)
# print(f"The mach number of the accelerated test gas is {M_1_2:.2f} ")

print("---" * 5)

#Post shock flow properties
