# from scipy import integrate
# import numpy as np
#
# def integrate_function(f, a, b):
#     """
#     Compute the definite integral of f from a to b using SciPy's integration methods.
#
#     f : the function to integrate.
#     a : the lower limit of integration.
#     b : the upper limit of integration.
#
#     Returns the integral and the estimated error.
#     """
#     result, error = integrate.quad(f, a, b)
#     return result, error
#
# # Example usage:
# result, error = integrate_function(np.sin, 0, np.pi)
# print("Integral of sin(x) from 0 to pi:", result)
# print("Estimated error:", error)

from scipy import integrate
import numpy as np

# Function to return a constant value
def constant(x):
    return 8.47e3  # constant value

#FUNCTION LAYOUT
def integrate_function(f, a, b):
    """
    Compute the definite integral of f from a to b using SciPy's integration methods.

    f : the function to integrate.
    a : the lower limit of integration.
    b : the upper limit of integration.

    Returns the integral and the estimated error.
    """
    result, error = integrate.quad(f, a, b)
    return result, error

# Working Area
print("-"* 10)
result, error = integrate_function(constant, 0, 0.0506)
print(f"Integral between the lower and upper limit is {result:.2f} N/m" )
print("Estimated error:", error)


# def Drag(result,angle):
#     return result * np.sin(angle)
# drag1 = Drag(result,9)
# print(f"Drag at Point 1 is {drag1:.2f} N/m")
#
def Drag(result, angle_deg):
    angle_rad = np.radians(angle_deg)
    return result * np.sin(angle_rad)

angle_deg = 9  # angle in degrees
drag1 = Drag(result, angle_deg)
print(f"Drag at Point 1 is {drag1:.2f} N/m")

print("-"* 10)
def P2(x):

    return 24.84e3


ULP2 = 0.05 / np.sin(np.radians(18))
result2,error2 = integrate_function(P2,0,ULP2)
print(f"Upper limit P2 is {ULP2:.3f} m")
print(f"Integral at Point 2 is {result2:.2f} N/m")

def sutherlands_law(T, mu_0=1.789e-5, T_0=288, S=110):
    """
    Calculate the dynamic viscosity of air using Sutherland's Law.

    Parameters:
    T : float
        Temperature at which viscosity is calculated (in Kelvin).
    mu_0 : float, optional
        Reference viscosity (kg/m.s) at reference temperature T_0.
    T_0 : float, optional
        Reference temperature (in Kelvin).
    S : float, optional
        Sutherland's constant (in Kelvin).

    Returns:
    float
        Viscosity at temperature T (in kg/m.s).
    """
    return mu_0 * (T / T_0)**1.5 * (T_0 + S) / (T + S)


T1 = 350.95
T2 = 494.49
T3 = 251.6
T4 = 222.54 #ambient at that altitude
viscosity1 = sutherlands_law(T1)
viscosity2 = sutherlands_law(T2)
viscosity3 = sutherlands_law(T3)
viscosity4 = sutherlands_law(T4)
print(f"The dynamic viscosity of air at {T1} K is {viscosity1:.8f} kg/(m·s)")
print(f"The dynamic viscosity of air at {T2} K is {viscosity2:.8f} kg/(m·s)")
print(f"The dynamic viscosity of air at {T3} K is {viscosity3:.8f} kg/(m·s)")
print(f"The dynamic viscosity of air at {T4} K is {viscosity4:.8f} kg/(m·s)")

