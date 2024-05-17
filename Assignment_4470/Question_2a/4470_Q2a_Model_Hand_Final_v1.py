# Q2a Model and Hand calcs for Truncated Prism

import numpy as np
from stl import mesh

dynamic_pressure = 50000  # Pa
mach = 8
gamma = 1.289  # Mars CO2 atmosphere
R = 188.92  # CO2 in J/(kg·K)
T = 220  # Temperature
front_face_threshold = 0.99999  # Threshold to identify front-facing panels direction of flow
file_path = 'wk-07.stl'

# Functions
def CP_MAX(gamma):
    outer = 4/(gamma+1)
    inner = ((gamma+1)**2) / (4 * gamma)
    exp = gamma / (gamma - 1)
    return outer * (inner ** exp)

cp_max = CP_MAX(gamma)

def load_stl(file_path):
    return mesh.Mesh.from_file(file_path)

def calc_area(v0, v1, v2):
    """Area of a triangle from vertices."""
    return 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))

def calc_cp(theta, cp_max):
    return cp_max * np.sin(theta) ** 2

def calc_force(normal, area, cp, dynamic_pressure):
    """Force on triangle."""
    force_magnitude = dynamic_pressure * cp * area
    return force_magnitude * normal

def calc_sound_speed(gamma, R, T):
    """Sound speed."""
    return np.sqrt(gamma * R * T)

def calc_velocity(mach, sound_speed):
    """Velocity"""
    return mach * sound_speed

def calc_density(dynamic_pressure, velocity):
    """Density from dynamic pressure and velocity."""
    return 2 * dynamic_pressure / (velocity ** 2)

def calc_theta(normal, flow_direction):
    """Angle between the normal and the flow direction."""
    dot_product = np.dot(normal, flow_direction)
    magnitude_normal = np.linalg.norm(normal)
    magnitude_flow = np.linalg.norm(flow_direction)
    cos_theta = dot_product / (magnitude_normal * magnitude_flow)
    theta = np.arccos(cos_theta)
    return theta

def is_front_facing(normal, flow_direction, threshold):
    """Check if the normal is front-facing with respect to the flow direction to apply cpmax."""
    dot_product = np.dot(normal, -flow_direction)
    return dot_product > threshold

def calc_total_drag_force(file_path, flow_direction):
    """Compute the total drag force for a given STL file."""
    your_mesh = load_stl(file_path)
    normals = your_mesh.normals
    vertices = your_mesh.vectors
    total_drag_force = np.array([0.0, 0.0, 0.0])
    for i, normal in enumerate(normals):
        normal = normal / np.linalg.norm(normal)
        theta = calc_theta(normal, flow_direction)
        theta_degrees = np.degrees(theta)

        # Calculate the effective inclination angle for inclined panels between cpmax and sin(theta)^2
        if not np.isclose(abs(normal[0]), 1.0):
            inclination_angle = np.abs(theta_degrees - 90)
        else:
            inclination_angle = 0

        inclination_angle_radians = np.radians(inclination_angle)

        # Determine cp based on the inclination angle
        if is_front_facing(normal, flow_direction, front_face_threshold):
            cp = cp_max  # Front-facing panels
        else:
            cp = calc_cp(inclination_angle_radians, cp_max)  # Inclined panels

        v0, v1, v2 = vertices[i]
        area = calc_area(v0, v1, v2)
        # Calculate the force on the triangle
        force = calc_force(normal, area, cp, dynamic_pressure)
        # Sum the applied forces
        total_drag_force += force

        print(f"Triangle {i + 1}:")
        print(f"Area: {area:.4f} m^2")
        print(f"Inclination Angle away from 90 (degrees): {inclination_angle:.2f}")
        print(f"Cp: {cp:.4f}")
        print()

    drag_force_x = abs(total_drag_force[0])  # Positive drag force

    return drag_force_x


sound_speed = calc_sound_speed(gamma, R, T)
velocity = calc_velocity(mach, sound_speed)
density = calc_density(dynamic_pressure, velocity)

print(f"Sound speed: {sound_speed:.2f} m/s")
print(f"Free-stream velocity: {velocity:.2f} m/s")
print(f"Density: {density:.4f} kg/m^3")
print()
flow_direction = np.array([1, 0, 0])
total_drag_force = calc_total_drag_force(file_path, flow_direction)

'''Hand Calcs'''
def drag_force(Cd,q,S):
    return Cd * q * S

Face1 = 0.5
Face2 = 0.9
Area1 = 0.5**2
Area2 = (Face2**2) - Area1
cp_inclined_hand = calc_cp(np.radians(9),cp_max)
print(f'Inclined cp hand {cp_inclined_hand:.5f}')
print(f"Total drag force in x-axis: {total_drag_force:.1f} N")

def CD(cp1,A1,cp2,A2):
    num = (cp1 * A1) + (cp2 * A2)
    denom = A1+A2
    return num /denom

CD_Q2a_lel = CD(cp_max,Area1,cp_inclined_hand,Area2)
drag_calc = drag_force(CD_Q2a_lel,dynamic_pressure,Area1+Area2)

print(f'CP max is {cp_max:.5f}')
print('***'*5)
print('Hand calcs')
print(f'Area 1 is {Area1:.2f}')
print(f'Area 2 is {Area2:.2f}')
print(f'CD_Q2a is {CD_Q2a_lel:.4f}')
print('***'*5)
print(f'Total drag force from hand calculations {drag_calc:.2f} N')
print(f"Total drag force from stl file: {total_drag_force:.2f} N")

percentage_error = np.abs((total_drag_force - drag_calc) / total_drag_force) * 100
print(f"Percentage Error: {percentage_error:.6f}%")

