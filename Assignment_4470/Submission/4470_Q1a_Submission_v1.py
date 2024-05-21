# Import necessary libraries
import numpy as np
from scipy.optimize import fsolve

# Given parameters for the initial state
U_inf = 5200  # Free-stream velocity in meters per second (m/s)
T_inf = 1000  # Free-stream temperature in Kelvin (K)
P_inf = 500  # Free-stream pressure in Pascals (Pa)
Diameter = 0.075  # Diameter in meters (m)
mol_N2 = 0.79  # Mole fraction of Nitrogen (N2)
mol_O2 = 0.21  # Mole fraction of Oxygen (O2)
Mw_O2 = 32  # Molecular weight of Oxygen (O2) in grams per mole (g/mol)
Mw_N2 = 28  # Molecular weight of Nitrogen (N2) in grams per mole (g/mol)
Mw_air = 28.97  # Molecular weight of air in grams per mole (g/mol)
Ru = 8314  # Universal gas constant in Joules per kilomole per Kelvin (J/(kmol*K))

# Specific vibrational temperatures for the Kruger model
thetaV_O2 = 2270  # Vibrational temperature for O2 in Kelvin (K)
thetaT_O2 = 11390  # Translational temperature for O2 in Kelvin (K)
thetaV_N2 = 3390  # Vibrational temperature for N2 in Kelvin (K)

# Calculate the gas constant for a given molecular weight
def R_gas(Mw):
    return Ru / Mw

R_air = R_gas(Mw_air)
print(f'Gas constant for air (R) is {R_air:.2f} J/(kg*K)')

# Function to calculate density using the ideal gas law
def ideal_gas_law_density(P, R, T):
    return P / (R * T)

rho1 = ideal_gas_law_density(P_inf, 287, T_inf)
print(f'Initial density (rho1) is {rho1:.8f} kg/m^3')

# Functions to calculate post-shock conditions based on the given inputs and formulas
def P2(P1, rho1, U1, rho2, U2):
    return P1 + rho1 * U1**2 - rho2 * U2**2

def h2(h1, U1, U2):
    return h1 + (0.5 * U1**2) - (0.5 * U2**2)

def U2(rho1, U1, rho2):
    return (rho1 * U1) / rho2

def rho2(P2, Z2, R, T2):
    return P2 / (Z2 * R * T2)

# Given specific heat capacities and ratio at 1000K from Cengel Thermodynamics Textbook

Cp_1000K = 1142  # Specific heat capacity at constant pressure (Cp) in J/(kg*K)
Cv_1000K = 855  # Specific heat capacity at constant volume (Cv) in J/(kg*K)
gamma_1000 = 1.336  # Specific heat ratio (gamma)

# Calculate the initial enthalpy
h1 = Cp_1000K * T_inf
print(f'Initial enthalpy (h1) is {h1:.2f} J/kg')

# Function to calculate momentum pressure
def momentum_pressure(P1, rho1, U1):
    return P1 + (rho1 * U1**2)

PMOM = momentum_pressure(P_inf, rho1, U_inf)
print(f'Momentum pressure (PMOM) is {PMOM:.2f} Pa')

# Conversion to atmospheres
PMOM_atm = PMOM / 101325
print(f'Momentum pressure in atmospheres (PMOM) is {PMOM_atm:.6f} atm')

# Function to calculate total enthalpy
def total_enthalpy(h1, U1):
    return h1 + (0.5 * U1**2)

HTOT = total_enthalpy(h1, U_inf)
print(f'Total enthalpy (HTOT) is {HTOT / 1e6:.2f} MJ/kg')

# Iteration 1 for post-shock conditions
T2_1 = 5750  # Estimated post-shock temperature in Kelvin (K)
z1_1 = 1.313  # Estimated compressibility factor

# Calculate density after the first iteration
rho2_1 = rho2(PMOM, z1_1, R_air, T2_1)
print(f'Post-shock density after first iteration (rho2_1) is {rho2_1:.8f} kg/m^3')

# Calculate velocity after the first iteration
U2_1 = U2(rho1, U_inf, rho2_1)
print(f'Post-shock velocity after first iteration (U2_1) is {U2_1:.2f} m/s')

# Iteration 2 for refined post-shock conditions
P2_2 = PMOM - rho2_1 * U2_1**2
print(f'Post-shock pressure after second iteration (P2_2) is {P2_2:.2f} Pa')

P2_2_atm = P2_2 / 101325
print(f'Post-shock pressure in atmospheres after second iteration (P2_2_atm) is {P2_2_atm:.6f} atm')

# Total enthalpy after the second iteration
h2_2 = HTOT - (0.5 * U2_1**2)
print(f'Total enthalpy after second iteration (h2_2) is {h2_2 / 1e6:.2f} MJ/kg')

T2_2 = 5687  # Refined post-shock temperature in Kelvin (K)
z2 = 1.314  # Refined compressibility factor

# Calculate density after the second iteration
rho2_2 = rho2(P2_2, z2, R_air, T2_2)
print(f'Post-shock density after second iteration (rho2_2) is {rho2_2:.8f} kg/m^3')

# Calculate velocity after the second iteration
U2_2 = U2(rho2_1, U2_1, rho2_2)
print(f'Post-shock velocity after second iteration (U2_2) is {U2_2:.2f} m/s')

# Iteration 3 for further refined post-shock conditions
P2_3 = PMOM - (rho2_2 * U2_2**2)
print(f'Post-shock pressure after third iteration (P2_3) is {P2_3:.2f} Pa')

PMOM_3 = P2_3 / 101325
print(f'Post-shock pressure in atmospheres after third iteration (PMOM_3) is {PMOM_3:.6f} atm')

# Total enthalpy after the third iteration
h2_3 = HTOT - (0.5 * U2_2**2)
print(f'Total enthalpy after third iteration (h2_3) is {h2_3 / 1e6:.3f} MJ/kg')

# Convergence results
print('\nConvergence results:')
print(f'The post-shock density is {rho2_2:.7f} kg/m^3')
print(f'The post-shock pressure is {PMOM_3:.6f} Pa')
print(f'The post-shock velocity is {U2_2:.2f} m/s')
print(f'The post-shock temperature is {T2_2:.2f} K')
print(f'The sensible post-shock enthalpy is {h2_3:.2f} J/kg')

# Gibbs free energy calculations for dissociated composition at post-shock conditions
N_HoT_0K = 111543.3  # Enthalpy of formation for N at 0K
N_neg_FT_Ho_5687K = 263171.6  # Negative formation enthalpy for N at 5687K

N2_HoT_0K = -2072.3  # Enthalpy of formation for N2 at 0K
N2_neg_FT_Ho_5687K = 346000.8  # Negative formation enthalpy for N2 at 5687K

O2_HoT_0K = -2074.7  # Enthalpy of formation for O2 at 0K
O2_neg_FT_Ho_5687K = 371000.1  # Negative formation enthalpy for O2 at 5687K

O_HoT_0K = 57949.1  # Enthalpy of formation for O at 0K
O_neg_FT_Ho_5687K = 274002.6  # Negative formation enthalpy for O at 5687K

NO_HoT_0K = 19403.0  # Enthalpy of formation for NO at 0K
NO_neg_FT_Ho_5687K = 374003.2  # Negative formation enthalpy for NO at 5687K

joules = 4.184  # Conversion factor from calories to joules

# Gibbs free energy calculations in Joules
N2FT = N2_HoT_0K - N2_neg_FT_Ho_5687K
NFT = N_HoT_0K - N_neg_FT_Ho_5687K
OFT = O_HoT_0K - O_neg_FT_Ho_5687K
O2FT = O2_HoT_0K - O2_neg_FT_Ho_5687K
NOFT = NO_HoT_0K - NO_neg_FT_Ho_5687K

N2FT_J = N2FT * joules
NFT_J = NFT * joules
OFT_J = OFT * joules
O2FT_J = O2FT * joules
NOFT_J = NOFT * joules

# Gibbs free energy for dissociation reactions
GibbsNitrogen = (2 * NFT_J) - N2FT_J
GibbsOxygen = (2 * OFT_J) - O2FT_J
GibbsNO = NOFT_J - NFT_J - OFT_J

# Function to calculate Gibbs free energy equilibrium constant
def gibbs(DeltaGibbs, T):
    exponent = 8.314 * T
    return np.exp(-DeltaGibbs / exponent)

# Calculate equilibrium constants
KN = gibbs(GibbsNitrogen, T2_2)
KO = gibbs(GibbsOxygen, T2_2)
KNO = gibbs(GibbsNO, T2_2)

print(f'Equilibrium constant for N2 dissociation (KN) is {KN:.7f} atm')
print(f'Equilibrium constant for O2 dissociation (KO) is {KO:.7f} atm')
print(f'Equilibrium constant for NO formation (KNO) is {KNO:.7f} atm')
