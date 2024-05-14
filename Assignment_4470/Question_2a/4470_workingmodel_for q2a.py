#
# import numpy as np
# from stl import mesh
#
# # Constants and Parameters
# cp_max = 1.877
# dynamic_pressure = 50000  # Pa (50 kPa)
# mach = 8
# gamma = 1.29  # Ratio of specific heats for Mars CO2 atmosphere
# R = 188.92  # Specific gas constant for CO2 in J/(kg·K)
# T = 220  # Temperature in Kelvin
# inclination_angle = np.radians(9)  # 9 degrees inclination in radians
#
# # Functions
#
# def load_stl(file_path):
#     """Load an STL file and return the mesh object."""
#     return mesh.Mesh.from_file(file_path)
#
# def calculate_area(v0, v1, v2):
#     """Calculate the area of a triangle given its vertices."""
#     return 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))
#
# def calculate_cp(theta, cp_max):
#     """Calculate the pressure coefficient for an inclined surface."""
#     return cp_max * np.sin(theta) ** 2
#
# def calculate_force(normal, area, cp, dynamic_pressure):
#     """Calculate the aerodynamic force on a triangle."""
#     force_magnitude = dynamic_pressure * cp * area
#     return force_magnitude * normal
#
# def compute_sound_speed(gamma, R, T):
#     """Compute the sound speed in the atmosphere."""
#     return np.sqrt(gamma * R * T)
#
# def compute_velocity(mach, sound_speed):
#     """Compute the velocity from Mach number and sound speed."""
#     return mach * sound_speed
#
# def compute_density(dynamic_pressure, velocity):
#     """Compute the density from dynamic pressure and velocity."""
#     return 2 * dynamic_pressure / (velocity ** 2)
#
# def calculate_theta(normal, flow_direction):
#     """Calculate the angle between the normal and the flow direction."""
#     dot_product = np.dot(normal, flow_direction)
#     magnitude_normal = np.linalg.norm(normal)
#     magnitude_flow = np.linalg.norm(flow_direction)
#     cos_theta = dot_product / (magnitude_normal * magnitude_flow)
#     theta = np.arccos(cos_theta)
#     return theta
#
# def compute_total_drag_force(file_path, flow_direction):
#     """Compute the total drag force for a given STL file."""
#     # Load the STL file
#     your_mesh = load_stl(file_path)
#
#     # Extract the vertices and normals
#     normals = your_mesh.normals
#     vertices = your_mesh.vectors
#
#     # Initialize total drag force
#     total_drag_force = np.array([0.0, 0.0, 0.0])
#
#     # Calculate the aerodynamic force on each triangle patch
#     for i, normal in enumerate(normals):
#         # Normalize the normal vector
#         normal = normal / np.linalg.norm(normal)
#
#         # Calculate the angle between the normal and the flow direction
#         theta = calculate_theta(normal, flow_direction)
#         theta_degrees = np.degrees(theta)
#
#         # Determine cp based on the triangle index and the angle
#         if i < 2:
#             cp = cp_max  # Front-facing panels
#         else:
#             cp = calculate_cp(inclination_angle, cp_max)  # Inclined panels
#
#         # Calculate area of the triangle
#         v0, v1, v2 = vertices[i]
#         area = calculate_area(v0, v1, v2)
#
#         # Calculate the force on the triangle
#         force = calculate_force(normal, area, cp, dynamic_pressure)
#
#         # Sum the forces
#         total_drag_force += force
#
#         # Print area, theta, and cp values for debugging
#         print(f"Triangle {i + 1}:")
#         print(f"Normal vector: {normal}")
#         print(f"Vertices: {v0}, {v1}, {v2}")
#         print(f"Area: {area:.4f}")
#         print(f"Theta (radians): {theta:.4f}, Theta (degrees): {theta_degrees:.2f}")
#         print(f"Cp: {cp:.4f}")
#         print(f"Force: {force}")
#         print()
#
#     # The drag force is the force in the x-direction
#     drag_force_x = abs(total_drag_force[0])  # Ensure positive drag force
#
#     return drag_force_x
#
# # Compute sound speed, velocity, and density
# sound_speed = compute_sound_speed(gamma, R, T)
# velocity = compute_velocity(mach, sound_speed)
# density = compute_density(dynamic_pressure, velocity)
#
# print(f"Sound speed: {sound_speed:.2f} m/s")
# print(f"Free-stream velocity: {velocity:.2f} m/s")
# print(f"Density: {density:.4f} kg/m^3")
# print()
#
# # Define the flow direction (example: along the x-axis)
# flow_direction = np.array([1, 0, 0])
#
# # Example usage
# file_path = 'wk-07.stl'
# total_drag_force = compute_total_drag_force(file_path, flow_direction)
# print(f"Total drag force in x-direction: {total_drag_force:.2f} N")

# Basic above.



import numpy as np
from stl import mesh

# Constants and Parameters
cp_max = 1.87733
dynamic_pressure = 50000  # Pa (50 kPa)
mach = 8
gamma = 1.2989  # Ratio of specific heats for Mars CO2 atmosphere
R = 188.92  # Specific gas constant for CO2 in J/(kg·K)
T = 220  # Temperature in Kelvin
front_face_threshold = 0.99999  # Threshold to identify front-facing panels

# Functions

def load_stl(file_path):
    """Load an STL file and return the mesh object."""
    return mesh.Mesh.from_file(file_path)

def calculate_area(v0, v1, v2):
    """Calculate the area of a triangle given its vertices."""
    return 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))

def calculate_cp(theta, cp_max):
    """Calculate the pressure coefficient for an inclined surface."""
    return cp_max * np.sin(theta) ** 2

def calculate_force(normal, area, cp, dynamic_pressure):
    """Calculate the aerodynamic force on a triangle."""
    force_magnitude = dynamic_pressure * cp * area
    return force_magnitude * normal

def compute_sound_speed(gamma, R, T):
    """Compute the sound speed in the atmosphere."""
    return np.sqrt(gamma * R * T)

def compute_velocity(mach, sound_speed):
    """Compute the velocity from Mach number and sound speed."""
    return mach * sound_speed

def compute_density(dynamic_pressure, velocity):
    """Compute the density from dynamic pressure and velocity."""
    return 2 * dynamic_pressure / (velocity ** 2)

def calculate_theta(normal, flow_direction):
    """Calculate the angle between the normal and the flow direction."""
    dot_product = np.dot(normal, flow_direction)
    magnitude_normal = np.linalg.norm(normal)
    magnitude_flow = np.linalg.norm(flow_direction)
    cos_theta = dot_product / (magnitude_normal * magnitude_flow)
    theta = np.arccos(cos_theta)
    return theta

def is_front_facing(normal, flow_direction, threshold):
    """Check if the normal is front-facing with respect to the flow direction."""
    dot_product = np.dot(normal, -flow_direction)
    return dot_product > threshold

def compute_total_drag_force(file_path, flow_direction):
    """Compute the total drag force for a given STL file."""
    # Load the STL file
    your_mesh = load_stl(file_path)

    # Extract the vertices and normals
    normals = your_mesh.normals
    vertices = your_mesh.vectors

    # Initialize total drag force
    total_drag_force = np.array([0.0, 0.0, 0.0])

    # Calculate the aerodynamic force on each triangle patch
    for i, normal in enumerate(normals):
        # Normalize the normal vector
        normal = normal / np.linalg.norm(normal)

        # Calculate the angle between the normal and the flow direction
        theta = calculate_theta(normal, flow_direction)
        theta_degrees = np.degrees(theta)

        # Calculate the effective inclination angle for inclined panels
        if not np.isclose(abs(normal[0]), 1.0):
            inclination_angle = np.abs(theta_degrees - 90)
        else:
            inclination_angle = 0

        inclination_angle_radians = np.radians(inclination_angle)

        # Determine cp based on the inclination angle
        if is_front_facing(normal, flow_direction, front_face_threshold):
            cp = cp_max  # Front-facing panels
        else:
            cp = calculate_cp(inclination_angle_radians, cp_max)  # Inclined panels

        # Calculate area of the triangle
        v0, v1, v2 = vertices[i]
        area = calculate_area(v0, v1, v2)

        # Calculate the force on the triangle
        force = calculate_force(normal, area, cp, dynamic_pressure)

        # Sum the forces
        total_drag_force += force

        # Print area, theta, and cp values for debugging
        print(f"Triangle {i + 1}:")
        # print(f"Normal vector: {normal}")
        # print(f"Vertices: {v0}, {v1}, {v2}")
        print(f"Area: {area:.4f} m^2")
        # print(f"Theta (degrees): {theta_degrees:.2f}")
        print(f"Inclination Angle (degrees): {inclination_angle:.2f}")
        print(f"Cp: {cp:.4f}")
        # print(f"Force: {force}")
        print()

    # The drag force is the force in the x-direction
    drag_force_x = abs(total_drag_force[0])  # Ensure positive drag force

    return drag_force_x

# Compute sound speed, velocity, and density
sound_speed = compute_sound_speed(gamma, R, T)
velocity = compute_velocity(mach, sound_speed)
density = compute_density(dynamic_pressure, velocity)

print(f"Sound speed: {sound_speed:.2f} m/s")
print(f"Free-stream velocity: {velocity:.2f} m/s")
print(f"Density: {density:.4f} kg/m^3")
print()

# Define the flow direction (example: along the x-axis)
flow_direction = np.array([1, 0, 0])

# Example usage
file_path = 'wk-07.stl'
total_drag_force = compute_total_drag_force(file_path, flow_direction)
print(f"Total drag force in x-direction: {total_drag_force:.1f} N")
print(f'The hand calculated drag force is 24745 N')


# advanced above