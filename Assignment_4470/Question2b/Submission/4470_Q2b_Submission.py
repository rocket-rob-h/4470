# import numpy as np
# from stl import mesh
#
# # Constants
# dynamic_pressure = 50000  # Pa (dynamic pressure)
# mach = 8  # Mach number
# gamma = 1.4  # Specific heat ratio
# R = 287  # Specific gas constant for air in J/(kg·K)
# T = 220  # Temperature in Kelvin
# Free_stream_pressure = 5  # Pa (free stream pressure)
# angle_of_attack_degrees = 7.52  # Angle of attack in degrees
# file_paths = ['top-s47661808.stl', 'bot-s47661808.stl']  # STL file paths
#
# # Function to calculate maximum coefficient of pressure
# def CP_MAX(gamma):
#     """Calculate maximum coefficient of pressure based on gamma."""
#     return (4 / (gamma + 1)) * ((gamma + 1) ** 2 / (4 * gamma)) ** (gamma / (gamma - 1))
#
# cp_max = CP_MAX(gamma)
#
# # Function to calculate the speed of sound
# def calc_sound_speed(gamma, R, T):
#     """Calculate the speed of sound."""
#     return np.sqrt(gamma * R * T)
#
# # Function to calculate the velocity
# def calc_velocity(mach, sound_speed):
#     """Calculate the velocity."""
#     return mach * sound_speed
#
# # Function to calculate the density from dynamic pressure and velocity
# def calc_density(dynamic_pressure, velocity):
#     """Calculate the density from dynamic pressure and velocity."""
#     return 2 * dynamic_pressure / (velocity ** 2)
#
# # Function to load the STL file and return the mesh object
# def load_stl(file_path):
#     """Load an STL file and return the mesh object."""
#     return mesh.Mesh.from_file(file_path)
#
# # Function to normalize a vector
# def normalize_vector(vector):
#     """Normalize a vector."""
#     return vector / np.linalg.norm(vector)
#
# # Function to calculate the area of a triangle from its vertices
# def calc_area(v0, v1, v2):
#     """Calculate the area of a triangle given its vertices."""
#     return 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))
#
# # Function to calculate sin(theta) using the velocity and normal vector
# def calc_sin_theta(velocity, normal):
#     """Calculate sin(theta) using the velocity and normal vector."""
#     return np.dot(normalize_vector(velocity), normal)
#
#
# # Function to rotate the flow direction by the angle of attack around the z-axis
# def rotate_flow(flow_dir, angle_degrees):
#     """Rotate the flow direction by the angle of attack around the z-axis."""
#     angle_radians = np.radians(angle_degrees)
#     # Rotation matrix for rotation around the z-axis
#     rot_matrix = np.array([
#         [np.cos(angle_radians), -np.sin(angle_radians), 0],
#         [np.sin(angle_radians), np.cos(angle_radians), 0],
#         [0, 0, 1]
#     ])
#     return np.dot(rot_matrix, flow_dir)
#
# # Function to calculate the force on a triangle
# def calc_force_for_triangle(normal, vertices, velocity, dynamic_pressure, static_pressure):
#     """Calculate the force on a triangle given the normal, vertices, velocity, dynamic pressure, and static pressure."""
#     # Normalize the normal vector
#     normal = normalize_vector(normal)
#
#     # Calculate the area of the triangle
#     area = calc_area(*vertices)
#
#     # Calculate sin(theta) where theta is the angle between the flow direction and the panel's normal vector
#     sin_theta = calc_sin_theta(velocity, normal)
#
#     # Determine the pressure coefficient (Cp) and pressure based on the orientation of the panel
#     if sin_theta < 0:
#         # If sin(theta) is negative, the panel is in the flow field
#         cp = cp_max * (sin_theta ** 2)
#         pressure = dynamic_pressure
#     else:
#         # If sin(theta) is positive, the panel is in the shadow region
#         cp = cp_max * (sin_theta ** 2)
#         pressure = static_pressure
#
#     # Calculate the force on the panel
#     # Force is negative due to the task sheet specifications
#     force = -pressure * cp * area * normal
#
#     # Normalize the velocity vector
#     velocity_normalized = normalize_vector(velocity)
#
#     # Calculate the drag force component (force in the direction of the flow)
#     drag = np.dot(force, velocity_normalized)
#
#     # Calculate the lift force component (force perpendicular to the flow direction)
#     # Assuming the lift direction is the y-axis component of the normalized velocity
#     lift = np.dot(force,
#                   np.array([velocity_normalized[1], velocity_normalized[0], 0]))  # Correcting the lift direction
#
#     return drag, lift, sin_theta, cp, pressure
#
# # lift is total force time second index
#
# # Function to compute the aerodynamic forces for a given list of STL files
# def compute_aerodynamic_forces(file_paths, angle_of_attack_degrees):
#     """Compute the total drag and lift forces for the given STL files."""
#     # Initialize total drag and lift forces for top and bottom panels
#     total_drag_top = 0
#     total_lift_top = 0
#     total_drag_bottom = 0
#     total_lift_bottom = 0
#
#     # Calculate sound speed, velocity, and density
#     sound_speed = calc_sound_speed(gamma, R, T)
#     velocity = calc_velocity(mach, sound_speed)
#     rho = calc_density(dynamic_pressure, velocity)
#
#     print()
#     print('*** Initial Conditions ***' * 3)
#     print()
#     print(f"At {T} Kelvin, Calculated Sound Speed: {sound_speed:.2f} m/s")
#     print(f"At Mach {mach}, Calculated Velocity: {velocity:.2f} m/s")
#     print(f"Therefore, the Air Density is: {rho:.5f} kg/m^3")
#     print(f"At {dynamic_pressure} Pa, dynamic pressure on the vehicle.")
#
#     # Rotate the flow direction based on the angle of attack
#     adjusted_velocity = rotate_flow(np.array([velocity, 0, 0]), angle_of_attack_degrees)
#     print(
#         f"The {angle_of_attack_degrees} degree angle of attack causes a rotated Velocity Components: {adjusted_velocity} around the Z-axis.")
#
#     # Loop through each STL file and compute the forces on each panel
#     for file_path in file_paths:
#         mesh_data = load_stl(file_path)
#         for normal, vertices in zip(mesh_data.normals, mesh_data.vectors):
#             drag, lift, sin_theta, cp, pressure = calc_force_for_triangle(normal, vertices, adjusted_velocity,
#                                                                           dynamic_pressure, Free_stream_pressure)
#
#             # Accumulate the drag and lift forces separately for top and bottom panels
#             if 'top' in file_path:
#                 total_drag_top += drag
#                 total_lift_top += lift  # Separating top and bottom files for analysis of Q2c
#             else:
#                 total_drag_bottom += drag
#                 total_lift_bottom += lift
#
#     print()
#     print('*** Results ***' * 5)
#     print()
#     print(
#         f"Top Panel - Total Drag Force: {total_drag_top:.2f} N, Total Lift Force: {total_lift_top:.2f} N, Lift/Drag Ratio: {total_lift_top / total_drag_top:.2f}")
#     print(
#         f"Bottom Panel - Total Drag Force: {total_drag_bottom:.2f} N, Total Lift Force: {total_lift_bottom:.2f} N, Lift/Drag Ratio: {total_lift_bottom / total_drag_bottom:.2f}")
#
#     # Calculate combined drag and lift forces
#     combined_drag = total_drag_top + total_drag_bottom
#     combined_lift = total_lift_bottom - total_lift_top
#
#     print()
#     print(
#         f"Combined - Total Drag Force: {combined_drag:.2f} N, Total Lift Force: {combined_lift:.2f} N, Lift/Drag Ratio: {combined_lift / combined_drag:.2f}")
#
#
# # Run the program to compute lift and drag
# compute_aerodynamic_forces(file_paths, angle_of_attack_degrees)

# import numpy as np
# from stl import mesh
#
# # Constants
# dynamic_pressure = 50000  # Pa (dynamic pressure)
# mach = 8  # Mach number
# gamma = 1.4  # Specific heat ratio
# R = 287  # Specific gas constant for air in J/(kg·K)
# T = 220  # Temperature in Kelvin
# Free_stream_pressure = 5  # Pa (free stream pressure)
# angle_of_attack_degrees = 7.52  # Angle of attack in degrees
# file_paths = ['top-s47661808.stl', 'bot-s47661808.stl']  # STL file paths
#
# # Function to calculate maximum coefficient of pressure
# def CP_MAX(gamma):
#     """Calculate maximum coefficient of pressure based on gamma."""
#     return (4 / (gamma + 1)) * ((gamma + 1) ** 2 / (4 * gamma)) ** (gamma / (gamma - 1))
#
# cp_max = CP_MAX(gamma)
#
# # Function to calculate the speed of sound
# def calc_sound_speed(gamma, R, T):
#     """Calculate the speed of sound."""
#     return np.sqrt(gamma * R * T)
#
# # Function to calculate the velocity
# def calc_velocity(mach, sound_speed):
#     """Calculate the velocity."""
#     return mach * sound_speed
#
# # Function to calculate the density from dynamic pressure and velocity
# def calc_density(dynamic_pressure, velocity):
#     """Calculate the density from dynamic pressure and velocity."""
#     return 2 * dynamic_pressure / (velocity ** 2)
#
# # Function to load the STL file and return the mesh object
# def load_stl(file_path):
#     """Load an STL file and return the mesh object."""
#     return mesh.Mesh.from_file(file_path)
#
# # Function to normalize a vector
# def normalize_vector(vector):
#     """Normalize a vector."""
#     return vector / np.linalg.norm(vector)
#
# # Function to calculate the area of a triangle from its vertices
# def calc_area(v0, v1, v2):
#     """Calculate the area of a triangle given its vertices."""
#     v0, v1, v2 = np.array(v0), np.array(v1), np.array(v2)
#     return 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))
#
# # Function to calculate sin(theta) using the velocity and normal vector
# def calc_sin_theta(velocity, normal):
#     """Calculate sin(theta) using the velocity and normal vector."""
#     return np.dot(normalize_vector(velocity), normal)
#
# # Function to rotate the flow direction by the angle of attack around the z-axis
# def rotate_flow(flow_dir, angle_degrees):
#     """Rotate the flow direction by the angle of attack around the z-axis."""
#     angle_radians = np.radians(angle_degrees)
#     # Rotation matrix for rotation around the z-axis
#     rot_matrix = np.array([
#         [np.cos(angle_radians), -np.sin(angle_radians), 0],
#         [np.sin(angle_radians), np.cos(angle_radians), 0],
#         [0, 0, 1]
#     ])
#     return np.dot(rot_matrix, flow_dir)
#
# # Function to calculate the force on a triangle
# def calc_force_for_triangle(normal, vertices, velocity, dynamic_pressure, static_pressure):
#     """Calculate the force on a triangle given the normal, vertices, velocity, dynamic pressure, and static pressure."""
#     # Normalize the normal vector
#     normal = normalize_vector(normal)
#
#     # Calculate the area of the triangle
#     area = calc_area(*vertices)
#
#     # Calculate sin(theta) where theta is the angle between the flow direction and the panel's normal vector
#     sin_theta = calc_sin_theta(velocity, normal)
#
#     # Determine the pressure coefficient (Cp) and pressure based on the orientation of the panel
#     if sin_theta < 0:
#         # If sin(theta) is negative, the panel is in the flow field
#         cp = cp_max * (sin_theta ** 2)
#         pressure = dynamic_pressure
#     else:
#         # If sin(theta) is positive, the panel is in the shadow region
#         cp = cp_max * (sin_theta ** 2)
#         pressure = static_pressure
#
#     # Calculate the force on the panel
#     # Force is negative due to the task sheet specifications
#     force = -pressure * cp * area * normal
#
#     # Normalize the velocity vector
#     velocity_normalized = normalize_vector(velocity)
#
#     # Calculate the drag force component (force in the direction of the flow)
#     drag = np.dot(force, velocity_normalized)
#
#     # Calculate the lift force component (force perpendicular to the flow direction)
#     # Assuming the lift direction is the y-axis component of the normalized velocity
#     lift = np.dot(force, np.array([velocity_normalized[1], velocity_normalized[0], 0]))  # Correcting the lift direction
#
#     return drag, lift, sin_theta, cp, pressure
#
# # Function to compute the aerodynamic forces for a given list of STL files
# def compute_aerodynamic_forces(file_paths, angle_of_attack_degrees):
#     """Compute the total drag and lift forces for the given STL files."""
#     # Initialize total drag and lift forces for top and bottom panels
#     total_drag_top = 0
#     total_lift_top = 0
#     total_drag_bottom = 0
#     total_lift_bottom = 0
#
#     # Calculate sound speed, velocity, and density
#     sound_speed = calc_sound_speed(gamma, R, T)
#     velocity = calc_velocity(mach, sound_speed)
#     rho = calc_density(dynamic_pressure, velocity)
#
#     print()
#     print('*** Initial Conditions ***' * 3)
#     print()
#     print(f"At {T} Kelvin, Calculated Sound Speed: {sound_speed:.2f} m/s")
#     print(f"At Mach {mach}, Calculated Velocity: {velocity:.2f} m/s")
#     print(f"Therefore, the Air Density is: {rho:.5f} kg/m^3")
#     print(f"At {dynamic_pressure} Pa, dynamic pressure on the vehicle.")
#
#     # Rotate the flow direction based on the angle of attack
#     adjusted_velocity = rotate_flow(np.array([velocity, 0, 0]), angle_of_attack_degrees)
#     print(f"The {angle_of_attack_degrees} degree angle of attack causes a rotated Velocity Components: {adjusted_velocity} around the Z-axis.")
#
#     # Loop through each STL file and compute the forces on each panel
#     for file_path in file_paths:
#         mesh_data = load_stl(file_path)
#         for i, (normal, vertices) in enumerate(zip(mesh_data.normals, mesh_data.vectors)):
#             print(f"Normal {i}: {normal}")  # Print the normals
#             drag, lift, sin_theta, cp, pressure = calc_force_for_triangle(normal, vertices, adjusted_velocity, dynamic_pressure, Free_stream_pressure)
#             if 'top' in file_path:
#                 total_drag_top += drag
#                 total_lift_top += lift  # Separating top and bottom files for analysis of Q2c
#             else:
#                 total_drag_bottom += drag
#                 total_lift_bottom += lift
#             print(f"Normal {i}: {normal}, Area: {calc_area(*vertices):.4f}, Sin(theta): {sin_theta:.4f}, Cp: {cp:.4f}, Pressure: {pressure:.2f} Pa, Drag: {drag:.2f}, Lift: {lift:.2f}")
#
#     print()
#     print('*** Results ***' * 5)
#     print()
#     print(f"Top Panel - Total Drag Force: {total_drag_top:.2f} N, Total Lift Force: {total_lift_top:.2f} N, Lift/Drag Ratio: {total_lift_top / total_drag_top:.2f}")
#     print(f"Bottom Panel - Total Drag Force: {total_drag_bottom:.2f} N, Total Lift Force: {total_lift_bottom:.2f} N, Lift/Drag Ratio: {total_lift_bottom / total_drag_bottom:.2f}")
#
#     # Calculate combined drag and lift forces
#     combined_drag = total_drag_top + total_drag_bottom
#     combined_lift = total_lift_bottom - total_lift_top
#
#     print()
#     print(f"Combined - Total Drag Force: {combined_drag:.2f} N, Total Lift Force: {combined_lift:.2f} N, Lift/Drag Ratio: {combined_lift / combined_drag:.2f}")
#
# # Run the program to compute lift and drag
# compute_aerodynamic_forces(file_paths, angle_of_attack_degrees)
#
# # Sound speed, velocity, and density calculations
# sound_speed = calc_sound_speed(gamma, R, T)
# velocity = calc_velocity(mach, sound_speed)
# rho = calc_density(dynamic_pressure, velocity)
# print(f"\nDensity: {rho:.4f} kg/m^3")

'''CORRECT NORMALS'''
# import numpy as np
#
# # Constants
# dynamic_pressure = 50000  # Pa (dynamic pressure)
# mach = 8  # Mach number
# gamma = 1.4  # Specific heat ratio
# R = 287  # Specific gas constant for air in J/(kg·K)
# T = 220  # Temperature in Kelvin
# Free_stream_pressure = 5  # Pa (free stream pressure)
# angle_of_attack_degrees = 7.52  # Angle of attack in degrees
# file_paths = ['top-s47661808.stl', 'bot-s47661808.stl']  # STL file paths
#
# # Function to calculate maximum coefficient of pressure
# def CP_MAX(gamma):
#     """Calculate maximum coefficient of pressure based on gamma."""
#     return (4 / (gamma + 1)) * ((gamma + 1) ** 2 / (4 * gamma)) ** (gamma / (gamma - 1))
#
# cp_max = CP_MAX(gamma)
#
# # Function to calculate the speed of sound
# def calc_sound_speed(gamma, R, T):
#     """Calculate the speed of sound."""
#     return np.sqrt(gamma * R * T)
#
# # Function to calculate the velocity
# def calc_velocity(mach, sound_speed):
#     """Calculate the velocity."""
#     return mach * sound_speed
#
# # Function to calculate the density from dynamic pressure and velocity
# def calc_density(dynamic_pressure, velocity):
#     """Calculate the density from dynamic pressure and velocity."""
#     return 2 * dynamic_pressure / (velocity ** 2)
#
# # Function to normalize a vector
# def normalize_vector(vector):
#     """Normalize a vector."""
#     return vector / np.linalg.norm(vector)
#
# # Function to calculate the area of a triangle from its vertices
# def calc_area(v0, v1, v2):
#     """Calculate the area of a triangle given its vertices."""
#     v0, v1, v2 = np.array(v0), np.array(v1), np.array(v2)
#     return 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))
#
# # Function to calculate sin(theta) using the velocity and normal vector
# def calc_sin_theta(velocity, normal):
#     """Calculate sin(theta) using the velocity and normal vector."""
#     return np.dot(normalize_vector(velocity), normal)
#
# # Function to rotate the flow direction by the angle of attack around the z-axis
# def rotate_flow(flow_dir, angle_degrees):
#     """Rotate the flow direction by the angle of attack around the z-axis."""
#     angle_radians = np.radians(angle_degrees)
#     # Rotation matrix for rotation around the z-axis
#     rot_matrix = np.array([
#         [np.cos(angle_radians), -np.sin(angle_radians), 0],
#         [np.sin(angle_radians), np.cos(angle_radians), 0],
#         [0, 0, 1]
#     ])
#     return np.dot(rot_matrix, flow_dir)
#
# # Function to parse an STL file to extract facet normals and vertices
# def parse_stl(file_path):
#     """Parse an STL file to extract facet normals and vertices."""
#     with open(file_path, 'r') as file:
#         lines = file.readlines()
#
#     normals = []
#     vertices = []
#     current_vertices = []
#     for line in lines:
#         parts = line.strip().split()
#         if parts[0] == 'facet' and parts[1] == 'normal':
#             # Parse the normal vector
#             normal = list(map(float, parts[2:]))
#             normals.append(normal)
#         elif parts[0] == 'vertex':
#             # Parse the vertex
#             vertex = list(map(float, parts[1:]))
#             current_vertices.append(vertex)
#         elif parts[0] == 'endloop':
#             vertices.append(current_vertices)
#             current_vertices = []
#
#     return normals, vertices
#
# # Function to calculate the force on a triangle
# def calc_force_for_triangle(normal, vertices, velocity, dynamic_pressure, static_pressure):
#     """Calculate the force on a triangle given the normal, vertices, velocity, dynamic pressure, and static pressure."""
#     # Normalize the normal vector
#     normal = normalize_vector(normal)
#
#     # Calculate the area of the triangle
#     area = calc_area(*vertices)
#
#     # Calculate sin(theta) where theta is the angle between the flow direction and the panel's normal vector
#     sin_theta = calc_sin_theta(velocity, normal)
#
#     # Determine the pressure coefficient (Cp) and pressure based on the orientation of the panel
#     if sin_theta < 0:
#         # If sin(theta) is negative, the panel is in the flow field
#         cp = cp_max * (sin_theta ** 2)
#         pressure = dynamic_pressure
#     else:
#         # If sin(theta) is positive, the panel is in the shadow region
#         cp = cp_max * (sin_theta ** 2)
#         pressure = static_pressure
#
#     # Calculate the force on the panel
#     # Force is negative due to the task sheet specifications
#     force = -pressure * cp * area * normal
#
#     # Normalize the velocity vector
#     velocity_normalized = normalize_vector(velocity)
#
#     # Calculate the drag force component (force in the direction of the flow)
#     drag = np.dot(force, velocity_normalized)
#
#     # Calculate the lift force component (force perpendicular to the flow direction)
#     # Assuming the lift direction is the y-axis component of the normalized velocity
#     lift = np.dot(force, np.array([velocity_normalized[1], velocity_normalized[0], 0]))  # Correcting the lift direction
#
#     return drag, lift, sin_theta, cp, pressure
#
# # Function to compute the aerodynamic forces for a given list of STL files
# def compute_aerodynamic_forces(file_paths, angle_of_attack_degrees):
#     """Compute the total drag and lift forces for the given STL files."""
#     # Initialize total drag and lift forces for top and bottom panels
#     total_drag_top = 0
#     total_lift_top = 0
#     total_drag_bottom = 0
#     total_lift_bottom = 0
#
#     # Calculate sound speed, velocity, and density
#     sound_speed = calc_sound_speed(gamma, R, T)
#     velocity = calc_velocity(mach, sound_speed)
#     rho = calc_density(dynamic_pressure, velocity)
#
#     print()
#     print('*** Initial Conditions ***' * 3)
#     print()
#     print(f"At {T} Kelvin, Calculated Sound Speed: {sound_speed:.2f} m/s")
#     print(f"At Mach {mach}, Calculated Velocity: {velocity:.2f} m/s")
#     print(f"Therefore, the Air Density is: {rho:.5f} kg/m^3")
#     print(f"At {dynamic_pressure} Pa, dynamic pressure on the vehicle.")
#
#     # Rotate the flow direction based on the angle of attack
#     adjusted_velocity = rotate_flow(np.array([velocity, 0, 0]), angle_of_attack_degrees)
#     print(f"The {angle_of_attack_degrees} degree angle of attack causes a rotated Velocity Components: {adjusted_velocity} around the Z-axis.")
#
#     # Loop through each STL file and compute the forces on each panel
#     for file_path in file_paths:
#         normals, vertices = parse_stl(file_path)  # Parse the STL file to get normals and vertices
#         for i, (normal, vertex_set) in enumerate(zip(normals, vertices)):
#             print(f"Normal {i}: {normal}")  # Print the normals
#             drag, lift, sin_theta, cp, pressure = calc_force_for_triangle(normal, vertex_set, adjusted_velocity, dynamic_pressure, Free_stream_pressure)
#             if 'top' in file_path:
#                 total_drag_top += drag
#                 total_lift_top += lift  # Separating top and bottom files for analysis of Q2c
#             else:
#                 total_drag_bottom += drag
#                 total_lift_bottom += lift
#             print(f"Normal {i}: {normal}, Area: {calc_area(*vertex_set):.4f}, Sin(theta): {sin_theta:.4f}, Cp: {cp:.4f}, Pressure: {pressure:.2f} Pa, Drag: {drag:.2f}, Lift: {lift:.2f}")
#
#     print()
#     print('*** Results ***' * 5)
#     print()
#     print(f"Top Panel - Total Drag Force: {total_drag_top:.2f} N, Total Lift Force: {total_lift_top:.2f} N, Lift/Drag Ratio: {total_lift_top / total_drag_top:.2f}")
#     print(f"Bottom Panel - Total Drag Force: {total_drag_bottom:.2f} N, Total Lift Force: {total_lift_bottom:.2f} N, Lift/Drag Ratio: {total_lift_bottom / total_drag_bottom:.2f}")
#
#     # Calculate combined drag and lift forces
#     combined_drag = total_drag_top + total_drag_bottom
#     combined_lift = total_lift_bottom - total_lift_top
#
#     print()
#     print(f"Combined - Total Drag Force: {combined_drag:.2f} N, Total Lift Force: {combined_lift:.2f} N, Lift/Drag Ratio: {combined_lift / combined_drag:.2f}")
#
# # Run the program to compute lift and drag
# compute_aerodynamic_forces(file_paths, angle_of_attack_degrees)
#
# # Sound speed, velocity, and density calculations
# sound_speed = calc_sound_speed(gamma, R, T)
# velocity = calc_velocity(mach, sound_speed)
# rho = calc_density(dynamic_pressure, velocity)
# print(f"\nDensity: {rho:.4f} kg/m^3")

'''CORRECT NORMALS'''

import numpy as np

# Constants
dynamic_pressure = 50000  # Pa (dynamic pressure)
mach = 8  # Mach number
gamma = 1.4  # Specific heat ratio
R = 287  # Specific gas constant for air in J/(kg·K)
T = 220  # Temperature in Kelvin
Free_stream_pressure = 5  # Pa (free stream pressure)
angle_of_attack_degrees = 7.52  # Angle of attack in degrees
file_paths = ['top-s47661808.stl', 'bot-s47661808.stl']  # STL file paths

# Function to calculate maximum coefficient of pressure
def CP_MAX(gamma):
    """Calculate maximum coefficient of pressure based on gamma."""
    return (4 / (gamma + 1)) * ((gamma + 1) ** 2 / (4 * gamma)) ** (gamma / (gamma - 1))

cp_max = CP_MAX(gamma)

# Function to calculate the speed of sound
def calc_sound_speed(gamma, R, T):
    """Calculate the speed of sound."""
    return np.sqrt(gamma * R * T)

# Function to calculate the velocity
def calc_velocity(mach, sound_speed):
    """Calculate the velocity."""
    return mach * sound_speed

# Function to calculate the density from dynamic pressure and velocity
def calc_density(dynamic_pressure, velocity):
    """Calculate the density from dynamic pressure and velocity."""
    return 2 * dynamic_pressure / (velocity ** 2)

# Function to normalize a vector
def normalize_vector(vector):
    """Normalize a vector."""
    return vector / np.linalg.norm(vector)

# Function to calculate the area of a triangle from its vertices
def calc_area(v0, v1, v2):
    """Calculate the area of a triangle given its vertices."""
    v0, v1, v2 = np.array(v0), np.array(v1), np.array(v2)
    return 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))

# Function to calculate sin(theta) using the velocity and normal vector
def calc_sin_theta(velocity, normal):
    """Calculate sin(theta) using the velocity and normal vector."""
    return np.dot(normalize_vector(velocity), normal)

# Function to rotate the flow direction by the angle of attack around the z-axis
def rotate_flow(flow_dir, angle_degrees):
    """Rotate the flow direction by the angle of attack around the z-axis."""
    angle_radians = np.radians(angle_degrees)
    # Rotation matrix for rotation around the z-axis
    rot_matrix = np.array([
        [np.cos(angle_radians), -np.sin(angle_radians), 0],
        [np.sin(angle_radians), np.cos(angle_radians), 0],
        [0, 0, 1]
    ])
    return np.dot(rot_matrix, flow_dir)

# Function to parse an STL file to extract facet normals and vertices
def parse_stl(file_path):
    """Parse an STL file to extract facet normals and vertices."""
    with open(file_path, 'r') as file:
        lines = file.readlines()

    normals = []
    vertices = []
    current_vertices = []
    for line in lines:
        parts = line.strip().split()
        if parts[0] == 'facet' and parts[1] == 'normal':
            # Parse the normal vector
            normal = list(map(float, parts[2:]))
            normals.append(normal)
        elif parts[0] == 'vertex':
            # Parse the vertex
            vertex = list(map(float, parts[1:]))
            current_vertices.append(vertex)
        elif parts[0] == 'endloop':
            vertices.append(current_vertices)
            current_vertices = []

    return normals, vertices

# Function to calculate the force on a triangle
def calc_force_for_triangle(normal, vertices, velocity, dynamic_pressure, static_pressure):
    """Calculate the force on a triangle given the normal, vertices, velocity, dynamic pressure, and static pressure."""
    # Normalize the normal vector
    normal = normalize_vector(normal)

    # Calculate the area of the triangle
    area = calc_area(*vertices)

    # Calculate sin(theta) where theta is the angle between the flow direction and the panel's normal vector
    sin_theta = calc_sin_theta(velocity, normal)

    # Determine the pressure coefficient (Cp) and pressure based on the orientation of the panel
    if sin_theta < 0:
        # If sin(theta) is negative, the panel is in the flow field
        cp = cp_max * (sin_theta ** 2)
        pressure = dynamic_pressure
    else:
        # If sin(theta) is positive, the panel is in the shadow region
        cp = cp_max * (sin_theta ** 2)
        pressure = static_pressure

    # Calculate the force on the panel
    # Force is negative due to the task sheet specifications
    force = -pressure * cp * area * normal

    # Normalize the velocity vector
    velocity_normalized = normalize_vector(velocity)

    # Calculate the drag force component (force in the direction of the flow)
    drag = np.dot(force, velocity_normalized)

    # Calculate the lift force component (force perpendicular to the flow direction)
    # Assuming the lift direction is the y-axis component of the normalized velocity
    lift = np.dot(force, np.array([velocity_normalized[1], velocity_normalized[0], 0]))  # Correcting the lift direction

    return drag, lift, sin_theta, cp, pressure

# Function to compute the aerodynamic forces for a given list of STL files
def compute_aerodynamic_forces(file_paths, angle_of_attack_degrees):
    """Compute the total drag and lift forces for the given STL files."""
    # Initialize total drag and lift forces for top and bottom panels
    total_drag_top = 0
    total_lift_top = 0
    total_drag_bottom = 0
    total_lift_bottom = 0

    # Calculate sound speed, velocity, and density
    sound_speed = calc_sound_speed(gamma, R, T)
    velocity = calc_velocity(mach, sound_speed)
    rho = calc_density(dynamic_pressure, velocity)

    print()
    print('*** Initial Conditions ***' * 3)
    print()
    print(f"At {T} Kelvin, Calculated Sound Speed: {sound_speed:.2f} m/s")
    print(f"At Mach {mach}, Calculated Velocity: {velocity:.2f} m/s")
    print(f"Therefore, the Air Density is: {rho:.5f} kg/m^3")
    print(f"At {dynamic_pressure} Pa, dynamic pressure on the vehicle.")

    # Rotate the flow direction based on the angle of attack
    adjusted_velocity = rotate_flow(np.array([velocity, 0, 0]), angle_of_attack_degrees)
    print(f"The {angle_of_attack_degrees} degree angle of attack causes a rotated Velocity Components: {adjusted_velocity} around the Z-axis.")

    # Loop through each STL file and compute the forces on each panel
    for file_path in file_paths:
        normals, vertices = parse_stl(file_path)  # Parse the STL file to get normals and vertices
        for i, (normal, vertex_set) in enumerate(zip(normals, vertices)):
            print(f"Normal {i}: {normal}")  # Print the normals
            drag, lift, sin_theta, cp, pressure = calc_force_for_triangle(normal, vertex_set, adjusted_velocity, dynamic_pressure, Free_stream_pressure)
            if 'top' in file_path:
                total_drag_top += drag
                total_lift_top += lift  # Separating top and bottom files for analysis of Q2c
            else:
                total_drag_bottom += drag
                total_lift_bottom += lift
            print(f"Normal {i}: {normal}, Area: {calc_area(*vertex_set):.4f}, Sin(theta): {sin_theta:.4f}, Cp: {cp:.4f}, Pressure: {pressure:.2f} Pa, Drag: {drag:.2f}, Lift: {lift:.2f}")

    print()
    print('*** Results ***' * 5)
    print()
    print(f"Top Panel - Total Drag Force: {total_drag_top:.2f} N, Total Lift Force: {total_lift_top:.2f} N, Lift/Drag Ratio: {total_lift_top / total_drag_top:.2f}")
    print(f"Bottom Panel - Total Drag Force: {total_drag_bottom:.2f} N, Total Lift Force: {total_lift_bottom:.2f} N, Lift/Drag Ratio: {total_lift_bottom / total_drag_bottom:.2f}")

    # Calculate combined drag and lift forces
    combined_drag = total_drag_top + total_drag_bottom
    combined_lift =  total_lift_bottom - total_lift_top

    print()
    print(f"Combined - Total Drag Force: {combined_drag:.2f} N, Total Lift Force: {combined_lift:.2f} N, Lift/Drag Ratio: {combined_lift / combined_drag:.2f}")

# Run the program to compute lift and drag
compute_aerodynamic_forces(file_paths, angle_of_attack_degrees)
#
# import numpy as np
#
# # Constants
# dynamic_pressure = 50000  # Pa (dynamic pressure)
# mach = 8  # Mach number
# gamma = 1.4  # Specific heat ratio
# R = 287  # Specific gas constant for air in J/(kg·K)
# T = 220  # Temperature in Kelvin
# Free_stream_pressure = 5  # Pa (free stream pressure)
# angle_of_attack_degrees = 7.52  # Angle of attack in degrees
# file_paths = ['top-s47661808.stl', 'bot-s47661808.stl']  # STL file paths
#
# # Function to calculate maximum coefficient of pressure
# def CP_MAX(gamma):
#     """Calculate maximum coefficient of pressure based on gamma."""
#     return (4 / (gamma + 1)) * ((gamma + 1) ** 2 / (4 * gamma)) ** (gamma / (gamma - 1))
#
# cp_max = CP_MAX(gamma)
#
# # Function to calculate the speed of sound
# def calc_sound_speed(gamma, R, T):
#     """Calculate the speed of sound."""
#     return np.sqrt(gamma * R * T)
#
# # Function to calculate the velocity
# def calc_velocity(mach, sound_speed):
#     """Calculate the velocity."""
#     return mach * sound_speed
#
# # Function to calculate the density from dynamic pressure and velocity
# def calc_density(dynamic_pressure, velocity):
#     """Calculate the density from dynamic pressure and velocity."""
#     return 2 * dynamic_pressure / (velocity ** 2)
#
# # Function to normalize a vector
# def normalize_vector(vector):
#     """Normalize a vector."""
#     return vector / np.linalg.norm(vector)
#
# # Function to calculate the area of a triangle from its vertices
# def calc_area(v0, v1, v2):
#     """Calculate the area of a triangle given its vertices."""
#     v0, v1, v2 = np.array(v0), np.array(v1), np.array(v2)
#     return 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))
#
# # Function to calculate sin(theta) using the velocity and normal vector
# def calc_sin_theta(velocity, normal):
#     """Calculate sin(theta) using the velocity and normal vector."""
#     return np.dot(normalize_vector(velocity), normal)
#
# # Function to rotate the flow direction by the angle of attack around the z-axis
# def rotate_flow(flow_dir, angle_degrees):
#     """Rotate the flow direction by the angle of attack around the z-axis."""
#     angle_radians = np.radians(angle_degrees)
#     # Rotation matrix for rotation around the z-axis
#     rot_matrix = np.array([
#         [np.cos(angle_radians), -np.sin(angle_radians), 0],
#         [np.sin(angle_radians), np.cos(angle_radians), 0],
#         [0, 0, 1]
#     ])
#     return np.dot(rot_matrix, flow_dir)
#
# # Function to parse an STL file to extract facet normals and vertices
# def parse_stl(file_path):
#     """Parse an STL file to extract facet normals and vertices."""
#     with open(file_path, 'r') as file:
#         lines = file.readlines()
#
#     normals = []
#     vertices = []
#     current_vertices = []
#     for line in lines:
#         parts = line.strip().split()
#         if parts[0] == 'facet' and parts[1] == 'normal':
#             # Parse the normal vector
#             normal = list(map(float, parts[2:]))
#             normals.append(normal)
#         elif parts[0] == 'vertex':
#             # Parse the vertex
#             vertex = list(map(float, parts[1:]))
#             current_vertices.append(vertex)
#         elif parts[0] == 'endloop':
#             vertices.append(current_vertices)
#             current_vertices = []
#
#     return normals, vertices
#
# # Function to calculate the force on a triangle
# def calc_force_for_triangle(normal, vertices, velocity, dynamic_pressure, static_pressure):
#     """Calculate the force on a triangle given the normal, vertices, velocity, dynamic pressure, and static pressure."""
#     # Normalize the normal vector
#     normal = normalize_vector(normal)
#
#     # Calculate the area of the triangle
#     area = calc_area(*vertices)
#
#     # Calculate sin(theta) where theta is the angle between the flow direction and the panel's normal vector
#     sin_theta = calc_sin_theta(velocity, normal)
#
#     # Determine the pressure coefficient (Cp) and pressure based on the orientation of the panel
#     if sin_theta < 0:
#         # If sin(theta) is negative, the panel is in the flow field
#         cp = cp_max * (sin_theta ** 2)
#         pressure = dynamic_pressure
#     else:
#         # If sin(theta) is positive, the panel is in the shadow region
#         cp = cp_max * (sin_theta ** 2)
#         pressure = static_pressure
#
#     # Calculate the force on the panel
#     force = -pressure * cp * area * normal
#
#     # Calculate drag and lift forces using the normal components
#     drag = force[0]  # Using x-component of the normal for drag
#     lift = force[1]  # Using y-component of the normal for lift
#
#     return drag, lift, sin_theta, cp, pressure
#
# # Function to compute the aerodynamic forces for a given list of STL files
# def compute_aerodynamic_forces(file_paths, angle_of_attack_degrees):
#     """Compute the total drag and lift forces for the given STL files."""
#     # Initialize total drag and lift forces for top and bottom panels
#     total_drag_top = 0
#     total_lift_top = 0
#     total_drag_bottom = 0
#     total_lift_bottom = 0
#
#     # Calculate sound speed, velocity, and density
#     sound_speed = calc_sound_speed(gamma, R, T)
#     velocity = calc_velocity(mach, sound_speed)
#     rho = calc_density(dynamic_pressure, velocity)
#
#     print()
#     print('*** Initial Conditions ***' * 3)
#     print()
#     print(f"At {T} Kelvin, Calculated Sound Speed: {sound_speed:.2f} m/s")
#     print(f"At Mach {mach}, Calculated Velocity: {velocity:.2f} m/s")
#     print(f"Therefore, the Air Density is: {rho:.5f} kg/m^3")
#     print(f"At {dynamic_pressure} Pa, dynamic pressure on the vehicle.")
#
#     # Rotate the flow direction based on the angle of attack
#     adjusted_velocity = rotate_flow(np.array([velocity, 0, 0]), angle_of_attack_degrees)
#     print(f"The {angle_of_attack_degrees} degree angle of attack causes a rotated Velocity Components: {adjusted_velocity} around the Z-axis.")
#
#     # Loop through each STL file and compute the forces on each panel
#     for file_path in file_paths:
#         normals, vertices = parse_stl(file_path)  # Parse the STL file to get normals and vertices
#         for i, (normal, vertex_set) in enumerate(zip(normals, vertices)):
#             print(f"Normal {i}: {normal}")  # Print the normals
#             drag, lift, sin_theta, cp, pressure = calc_force_for_triangle(normal, vertex_set, adjusted_velocity, dynamic_pressure, Free_stream_pressure)
#             if 'top' in file_path:
#                 total_drag_top += drag
#                 total_lift_top += lift  # Separating top and bottom files for analysis of Q2c
#             else:
#                 total_drag_bottom += drag
#                 total_lift_bottom += lift
#             print(f"Normal {i}: {normal}, Area: {calc_area(*vertex_set):.4f}, Sin(theta): {sin_theta:.4f}, Cp: {cp:.4f}, Pressure: {pressure:.2f} Pa, Drag: {drag:.2f}, Lift: {lift:.2f}")
#
#     print()
#     print('*** Results ***' * 5)
#     print()
#     print(f"Top Panel - Total Drag Force: {total_drag_top:.2f} N, Total Lift Force: {total_lift_top:.2f} N, Lift/Drag Ratio: {total_lift_top / total_drag_top:.2f}")
#     print(f"Bottom Panel - Total Drag Force: {total_drag_bottom:.2f} N, Total Lift Force: {total_lift_bottom:.2f} N, Lift/Drag Ratio: {total_lift_bottom / total_drag_bottom:.2f}")
#
#     # Calculate combined drag and lift forces
#     combined_drag = total_drag_top + total_drag_bottom
#     combined_lift =  total_lift_bottom - total_lift_top
#
#     print()
#     print(f"Combined - Total Drag Force: {combined_drag:.2f} N, Total Lift Force: {combined_lift:.2f} N, Lift/Drag Ratio: {combined_lift / combined_drag:.2f}")
#
# # Run the program to compute lift and drag
# compute_aerodynamic_forces(file_paths, angle_of_attack_degrees)



# import numpy as np
#
# # Constants
# dynamic_pressure = 50000  # Pa (dynamic pressure)
# mach = 8  # Mach number
# gamma = 1.4  # Specific heat ratio
# R = 287  # Specific gas constant for air in J/(kg·K)
# T = 220  # Temperature in Kelvin
# Free_stream_pressure = 5  # Pa (free stream pressure)
# angle_of_attack_degrees = 7.52  # Angle of attack in degrees
# file_paths = ['top-s47661808.stl', 'bot-s47661808.stl']  # STL file paths
#
# # Function to calculate maximum coefficient of pressure
# def CP_MAX(gamma):
#     """Calculate maximum coefficient of pressure based on gamma."""
#     return (4 / (gamma + 1)) * ((gamma + 1) ** 2 / (4 * gamma)) ** (gamma / (gamma - 1))
#
# cp_max = CP_MAX(gamma)
#
# # Function to calculate the speed of sound
# def calc_sound_speed(gamma, R, T):
#     """Calculate the speed of sound."""
#     return np.sqrt(gamma * R * T)
#
# # Function to calculate the velocity
# def calc_velocity(mach, sound_speed):
#     """Calculate the velocity."""
#     return mach * sound_speed
#
# # Function to calculate the density from dynamic pressure and velocity
# def calc_density(dynamic_pressure, velocity):
#     """Calculate the density from dynamic pressure and velocity."""
#     return 2 * dynamic_pressure / (velocity ** 2)
#
# # Function to normalize a vector
# def normalize_vector(vector):
#     """Normalize a vector."""
#     return vector / np.linalg.norm(vector)
#
# # Function to calculate the area of a triangle from its vertices
# def calc_area(v0, v1, v2):
#     """Calculate the area of a triangle given its vertices."""
#     v0, v1, v2 = np.array(v0), np.array(v1), np.array(v2)
#     return 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))
#
# # Function to calculate sin(theta) using the velocity and normal vector
# def calc_sin_theta(velocity, normal):
#     """Calculate sin(theta) using the velocity and normal vector."""
#     return np.dot(normalize_vector(velocity), normal)
#
# # Function to rotate the flow direction by the angle of attack around the z-axis
# def rotate_flow(flow_dir, angle_degrees):
#     """Rotate the flow direction by the angle of attack around the z-axis."""
#     angle_radians = np.radians(angle_degrees)
#     # Rotation matrix for rotation around the z-axis
#     rot_matrix = np.array([
#         [np.cos(angle_radians), -np.sin(angle_radians), 0],
#         [np.sin(angle_radians), np.cos(angle_radians), 0],
#         [0, 0, 1]
#     ])
#     return np.dot(rot_matrix, flow_dir)
#
# # Function to parse an STL file to extract facet normals and vertices
# def parse_stl(file_path):
#     """Parse an STL file to extract facet normals and vertices."""
#     with open(file_path, 'r') as file:
#         lines = file.readlines()
#
#     normals = []
#     vertices = []
#     current_vertices = []
#     for line in lines:
#         parts = line.strip().split()
#         if parts[0] == 'facet' and parts[1] == 'normal':
#             # Parse the normal vector
#             normal = list(map(float, parts[2:]))
#             normals.append(normal)
#         elif parts[0] == 'vertex':
#             # Parse the vertex
#             vertex = list(map(float, parts[1:]))
#             current_vertices.append(vertex)
#         elif parts[0] == 'endloop':
#             vertices.append(current_vertices)
#             current_vertices = []
#
#     return normals, vertices
#
# # Function to calculate the force on a triangle
# def calc_force_for_triangle(normal, vertices, velocity, dynamic_pressure, static_pressure):
#     """Calculate the force on a triangle given the normal, vertices, velocity, dynamic pressure, and static pressure."""
#     # Normalize the normal vector
#     normal = normalize_vector(normal)
#
#     # Calculate the area of the triangle
#     area = calc_area(*vertices)
#
#     # Calculate sin(theta) where theta is the angle between the flow direction and the panel's normal vector
#     sin_theta = calc_sin_theta(velocity, normal)
#
#     # Determine the pressure coefficient (Cp) and pressure based on the orientation of the panel
#     if sin_theta < 0:
#         # If sin(theta) is negative, the panel is in the flow field
#         cp = cp_max * (sin_theta ** 2)
#         pressure = dynamic_pressure
#     else:
#         # If sin(theta) is positive, the panel is in the shadow region
#         cp = cp_max * (sin_theta ** 2)
#         pressure = static_pressure
#
#     # Calculate the force on the panel
#     force = -pressure * cp * area * normal
#
#     # Decompose the force into drag and lift components
#     drag = np.dot(force, normalize_vector(velocity))  # Component in the direction of the flow
#     lift = np.dot(force, np.cross(-normalize_vector(velocity), [0, 0, 1]))  # Perpendicular component in the plane
#
#     return drag, lift, sin_theta, cp, pressure
#
# # Function to compute the aerodynamic forces for a given list of STL files
# def compute_aerodynamic_forces(file_paths, angle_of_attack_degrees):
#     """Compute the total drag and lift forces for the given STL files."""
#     # Initialize total drag and lift forces for top and bottom panels
#     total_drag_top = 0
#     total_lift_top = 0
#     total_drag_bottom = 0
#     total_lift_bottom = 0
#
#     # Calculate sound speed, velocity, and density
#     sound_speed = calc_sound_speed(gamma, R, T)
#     velocity = calc_velocity(mach, sound_speed)
#     rho = calc_density(dynamic_pressure, velocity)
#
#     print()
#     print('*** Initial Conditions ***' * 3)
#     print()
#     print(f"At {T} Kelvin, Calculated Sound Speed: {sound_speed:.2f} m/s")
#     print(f"At Mach {mach}, Calculated Velocity: {velocity:.2f} m/s")
#     print(f"Therefore, the Air Density is: {rho:.5f} kg/m^3")
#     print(f"At {dynamic_pressure} Pa, dynamic pressure on the vehicle.")
#
#     # Rotate the flow direction based on the angle of attack
#     adjusted_velocity = rotate_flow(np.array([velocity, 0, 0]), angle_of_attack_degrees)
#     print(f"The {angle_of_attack_degrees} degree angle of attack causes a rotated Velocity Components: {adjusted_velocity} around the Z-axis.")
#
#     # Loop through each STL file and compute the forces on each panel
#     for file_path in file_paths:
#         normals, vertices = parse_stl(file_path)  # Parse the STL file to get normals and vertices
#         for i, (normal, vertex_set) in enumerate(zip(normals, vertices)):
#             print(f"Normal {i}: {normal}")  # Print the normals
#             drag, lift, sin_theta, cp, pressure = calc_force_for_triangle(normal, vertex_set, adjusted_velocity, dynamic_pressure, Free_stream_pressure)
#             if 'top' in file_path:
#                 total_drag_top += drag
#                 total_lift_top += lift  # Separating top and bottom files for analysis of Q2c
#             else:
#                 total_drag_bottom += drag
#                 total_lift_bottom += lift
#             print(f"Normal {i}: {normal}, Area: {calc_area(*vertex_set):.4f}, Sin(theta): {sin_theta:.4f}, Cp: {cp:.4f}, Pressure: {pressure:.2f} Pa, Drag: {drag:.2f}, Lift: {lift:.2f}")
#
#     print()
#     print('*** Results ***' * 5)
#     print()
#     print(f"Top Panel - Total Drag Force: {total_drag_top:.2f} N, Total Lift Force: {total_lift_top:.2f} N, Lift/Drag Ratio: {total_lift_top / total_drag_top:.2f}")
#     print(f"Bottom Panel - Total Drag Force: {total_drag_bottom:.2f} N, Total Lift Force: {total_lift_bottom:.2f} N, Lift/Drag Ratio: {total_lift_bottom / total_drag_bottom:.2f}")
#
#     # Calculate combined drag and lift forces
#     combined_drag = total_drag_top + total_drag_bottom
#     combined_lift = total_lift_top - total_lift_bottom
#
#     print()
#     print(f"Combined - Total Drag Force: {combined_drag:.2f} N, Total Lift Force: {combined_lift:.2f} N, Lift/Drag Ratio: {combined_lift / combined_drag:.2f}")
#
# # Run the program to compute lift and drag
# compute_aerodynamic_forces(file_paths, angle_of_attack_degrees)


# import numpy as np
#
# # Constants
# dynamic_pressure = 50000  # Pa (dynamic pressure)
# mach = 8  # Mach number
# gamma = 1.4  # Specific heat ratio
# R = 287  # Specific gas constant for air in J/(kg·K)
# T = 220  # Temperature in Kelvin
# Free_stream_pressure = 5  # Pa (free stream pressure)
# angle_of_attack_degrees = 7.52 # Angle of attack in degrees
# file_paths = ['top-s47661808.stl', 'bot-s47661808.stl']  # STL file paths
#
# # Function to calculate maximum coefficient of pressure
# def CP_MAX(gamma):
#     """Calculate maximum coefficient of pressure based on gamma."""
#     return (4 / (gamma + 1)) * ((gamma + 1) ** 2 / (4 * gamma)) ** (gamma / (gamma - 1))
#
# cp_max = CP_MAX(gamma)
#
# # Function to calculate the speed of sound
# def calc_sound_speed(gamma, R, T):
#     """Calculate the speed of sound."""
#     return np.sqrt(gamma * R * T)
#
# # Function to calculate the velocity
# def calc_velocity(mach, sound_speed):
#     """Calculate the velocity."""
#     return mach * sound_speed
#
# # Function to calculate the density from dynamic pressure and velocity
# def calc_density(dynamic_pressure, velocity):
#     """Calculate the density from dynamic pressure and velocity."""
#     return 2 * dynamic_pressure / (velocity ** 2)
#
# # Function to normalize a vector
# def normalize_vector(vector):
#     """Normalize a vector."""
#     return vector / np.linalg.norm(vector)
#
# # Function to calculate the area of a triangle from its vertices
# def calc_area(v0, v1, v2):
#     """Calculate the area of a triangle given its vertices."""
#     v0, v1, v2 = np.array(v0), np.array(v1), np.array(v2)
#     return 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))
#
# # Function to calculate sin(theta) using the velocity and normal vector
# def calc_sin_theta(velocity, normal):
#     """Calculate sin(theta) using the velocity and normal vector."""
#     return np.dot(normalize_vector(velocity), normal)
#
# # Function to rotate the flow direction by the angle of attack around the z-axis
# def rotate_flow(flow_dir, angle_degrees):
#     """Rotate the flow direction by the angle of attack around the z-axis."""
#     angle_radians = np.radians(angle_degrees)
#     # Rotation matrix for rotation around the z-axis
#     rot_matrix = np.array([
#         [np.cos(angle_radians), -np.sin(angle_radians), 0],
#         [np.sin(angle_radians), np.cos(angle_radians), 0],
#         [0, 0, 1]
#     ])
#     return np.dot(rot_matrix, flow_dir)
#
# # Function to parse an STL file to extract facet normals and vertices
# def parse_stl(file_path):
#     """Parse an STL file to extract facet normals and vertices."""
#     with open(file_path, 'r') as file:
#         lines = file.readlines()
#
#     normals = []
#     vertices = []
#     current_vertices = []
#     for line in lines:
#         parts = line.strip().split()
#         if parts[0] == 'facet' and parts[1] == 'normal':
#             # Parse the normal vector
#             normal = list(map(float, parts[2:]))
#             normals.append(normal)
#         elif parts[0] == 'vertex':
#             # Parse the vertex
#             vertex = list(map(float, parts[1:]))
#             current_vertices.append(vertex)
#         elif parts[0] == 'endloop':
#             vertices.append(current_vertices)
#             current_vertices = []
#
#     return normals, vertices
#
# # Function to calculate the force on a triangle
# def calc_force_for_triangle(normal, vertices, velocity, dynamic_pressure, static_pressure):
#     """Calculate the force on a triangle given the normal, vertices, velocity, dynamic pressure, and static pressure."""
#     # Normalize the normal vector
#     normal = normalize_vector(normal)
#
#     # Calculate the area of the triangle
#     area = calc_area(*vertices)
#
#     # Calculate sin(theta) where theta is the angle between the flow direction and the panel's normal vector
#     sin_theta = calc_sin_theta(velocity, normal)
#
#     # Determine the pressure coefficient (Cp) and pressure based on the orientation of the panel
#     if sin_theta < 0:
#         # If sin(theta) is negative, the panel is in the flow field
#         cp = cp_max * (sin_theta ** 2)
#         pressure = dynamic_pressure
#     else:
#         # If sin(theta) is positive, the panel is in the shadow region
#         cp = cp_max * (sin_theta ** 2)
#         pressure = static_pressure
#
#     # Calculate the force on the panel
#     force = -pressure * cp * area * normal
#
#     # Calculate drag and lift forces using the normal components
#     drag = force[0]  # Using x-component of the normal for drag
#     lift = force[1]   # Using y-component of the normal for lift
#
#     return drag, lift, sin_theta, cp, pressure
#
# # Function to compute the aerodynamic forces for a given list of STL files
# def compute_aerodynamic_forces(file_paths, angle_of_attack_degrees):
#     """Compute the total drag and lift forces for the given STL files."""
#     # Initialize total drag and lift forces for top and bottom panels
#     total_drag_top = 0
#     total_lift_top = 0
#     total_drag_bottom = 0
#     total_lift_bottom = 0
#
#     # Calculate sound speed, velocity, and density
#     sound_speed = calc_sound_speed(gamma, R, T)
#     velocity = calc_velocity(mach, sound_speed)
#     rho = calc_density(dynamic_pressure, velocity)
#
#     print()
#     print('*** Initial Conditions ***' * 3)
#     print()
#     print(f"At {T} Kelvin, Calculated Sound Speed: {sound_speed:.2f} m/s")
#     print(f"At Mach {mach}, Calculated Velocity: {velocity:.2f} m/s")
#     print(f"Therefore, the Air Density is: {rho:.5f} kg/m^3")
#     print(f"At {dynamic_pressure} Pa, dynamic pressure on the vehicle.")
#
#     # Rotate the flow direction based on the angle of attack
#     adjusted_velocity = rotate_flow(np.array([velocity, 0, 0]), angle_of_attack_degrees)
#     print(f"The {angle_of_attack_degrees} degree angle of attack causes a rotated Velocity Components: {adjusted_velocity} around the Z-axis.")
#
#     # Loop through each STL file and compute the forces on each panel
#     for file_path in file_paths:
#         normals, vertices = parse_stl(file_path)  # Parse the STL file to get normals and vertices
#         for i, (normal, vertex_set) in enumerate(zip(normals, vertices)):
#             print(f"Normal {i}: {normal}")  # Print the normals
#             drag, lift, sin_theta, cp, pressure = calc_force_for_triangle(normal, vertex_set, adjusted_velocity, dynamic_pressure, Free_stream_pressure)
#             if 'top' in file_path:
#                 total_drag_top += drag
#                 total_lift_top += lift  # Separating top and bottom files for analysis of Q2c
#             else:
#                 total_drag_bottom += drag
#                 total_lift_bottom += lift
#             print(f"Normal {i}: {normal}, Area: {calc_area(*vertex_set):.4f}, Sin(theta): {sin_theta:.4f}, Cp: {cp:.4f}, Pressure: {pressure:.2f} Pa, Drag: {drag:.2f}, Lift: {lift:.2f}")
#
#     print()
#     print('*** Results ***' * 5)
#     print()
#     print(f"Top Panel - Total Drag Force: {total_drag_top:.2f} N, Total Lift Force: {total_lift_top:.2f} N, Lift/Drag Ratio: {total_lift_top / total_drag_top:.2f}")
#     print(f"Bottom Panel - Total Drag Force: {total_drag_bottom:.2f} N, Total Lift Force: {total_lift_bottom:.2f} N, Lift/Drag Ratio: {total_lift_bottom / total_drag_bottom:.2f}")
#
#     # Calculate combined drag and lift forces
#     combined_drag = total_drag_top + total_drag_bottom
#     combined_lift =  total_lift_bottom - total_lift_top
#
#     print()
#     print(f"Combined - Total Drag Force: {combined_drag:.2f} N, Total Lift Force: {combined_lift:.2f} N, Lift/Drag Ratio: {combined_lift / combined_drag:.2f}")
#
# # Run the program to compute lift and drag
# compute_aerodynamic_forces(file_paths, angle_of_attack_degrees)


#
# #
# import numpy as np
#
# # Constants
# dynamic_pressure = 50000  # Pa (dynamic pressure)
# mach = 8  # Mach number
# gamma = 1.4  # Specific heat ratio
# R = 287  # Specific gas constant for air in J/(kg·K)
# T = 220  # Temperature in Kelvin
# Free_stream_pressure = 5  # Pa (free stream pressure)
# angle_of_attack_degrees = 7.52  # Angle of attack in degrees
# file_paths = ['top-s47661808.stl', 'bot-s47661808.stl']  # STL file paths
#
# # Function to calculate maximum coefficient of pressure
# def CP_MAX(gamma):
#     """Calculate maximum coefficient of pressure based on gamma."""
#     return (4 / (gamma + 1)) * ((gamma + 1) ** 2 / (4 * gamma)) ** (gamma / (gamma - 1))
#
# cp_max = CP_MAX(gamma)
#
# # Function to calculate the speed of sound
# def calc_sound_speed(gamma, R, T):
#     """Calculate the speed of sound."""
#     return np.sqrt(gamma * R * T)
#
# # Function to calculate the velocity
# def calc_velocity(mach, sound_speed):
#     """Calculate the velocity."""
#     return mach * sound_speed
#
# # Function to calculate the density from dynamic pressure and velocity
# def calc_density(dynamic_pressure, velocity):
#     """Calculate the density from dynamic pressure and velocity."""
#     return 2 * dynamic_pressure / (velocity ** 2)
#
# # Function to normalize a vector
# def normalize_vector(vector):
#     """Normalize a vector."""
#     return vector / np.linalg.norm(vector)
#
# # Function to calculate the area of a triangle from its vertices
# def calc_area(v0, v1, v2):
#     """Calculate the area of a triangle given its vertices."""
#     v0 = np.array(v0)
#     v1 = np.array(v1)
#     v2 = np.array(v2)
#     return 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))
#
# # Function to calculate sin(theta) using the velocity and normal vector
# def calc_sin_theta(velocity, normal):
#     """Calculate sin(theta) using the velocity and normal vector."""
#     return np.dot(normalize_vector(velocity), normal)
#
# # Function to rotate the flow direction by the angle of attack around the z-axis
# def rotate_flow(flow_dir, angle_degrees):
#     """Rotate the flow direction by the angle of attack around the z-axis."""
#     angle_radians = np.radians(angle_degrees)
#     # Rotation matrix for rotation around the z-axis
#     rot_matrix = np.array([
#         [np.cos(angle_radians), -np.sin(angle_radians), 0],
#         [np.sin(angle_radians), np.cos(angle_radians), 0],
#         [0, 0, 1]
#     ])
#     return np.dot(rot_matrix, flow_dir)
#
# # Function to calculate the force on a triangle
# def calc_force_for_triangle(normal, vertices, velocity, dynamic_pressure, static_pressure):
#     """Calculate the force on a triangle given the normal, vertices, velocity, dynamic pressure, and static pressure."""
#     # Normalize the normal vector
#     normal = normalize_vector(normal)
#
#     # Calculate the area of the triangle
#     area = calc_area(*vertices)
#
#     # Calculate sin(theta) where theta is the angle between the flow direction and the panel's normal vector
#     sin_theta = calc_sin_theta(velocity, normal)
#
#     # Determine the pressure coefficient (Cp) and pressure based on the orientation of the panel
#     if sin_theta < 0:
#         # If sin(theta) is negative, the panel is in the flow field
#         cp = cp_max * (sin_theta ** 2)
#         pressure = dynamic_pressure
#     else:
#         # If sin(theta) is positive, the panel is in the shadow region
#         cp = cp_max * (sin_theta ** 2)
#         pressure = static_pressure
#
#     # Calculate the force on the panel
#     # Force is negative due to the task sheet specifications
#     force = -pressure * cp * area * normal
#
#     # Normalize the velocity vector
#     velocity_normalized = normalize_vector(velocity)
#
#     # Calculate the drag force component (force in the direction of the flow)
#     drag = np.dot(force, velocity_normalized)
#
#     # Calculate the lift force component (force perpendicular to the flow direction)
#     # Assuming the lift direction is the y-axis component of the normalized velocity
#     lift = np.dot(force, np.array([velocity_normalized[1], velocity_normalized[0], 0]))  # Correcting the lift direction
#
#     return drag, lift, sin_theta, cp, pressure
#
# # Function to parse an STL file to extract facet normals and vertices
# def parse_stl(file_path):
#     """Parse an STL file to extract facet normals and vertices."""
#     with open(file_path, 'r') as file:
#         lines = file.readlines()
#
#     normals = []
#     vertices = []
#     current_vertices = []
#     for line in lines:
#         parts = line.strip().split()
#         if parts[0] == 'facet' and parts[1] == 'normal':
#             # Parse the normal vector
#             normal = list(map(float, parts[2:]))
#             normals.append(normal)
#         elif parts[0] == 'vertex':
#             # Parse the vertex
#             vertex = list(map(float, parts[1:]))
#             current_vertices.append(vertex)
#         elif parts[0] == 'endloop':
#             vertices.append(current_vertices)
#             current_vertices = []
#
#     return normals, vertices
#
# # Function to compute the aerodynamic forces for a given list of STL files
# def compute_aerodynamic_forces(file_paths, angle_of_attack_degrees):
#     """Compute the total drag and lift forces for the given STL files."""
#     # Initialize total drag and lift forces for top and bottom panels
#     total_drag_top = 0
#     total_lift_top = 0
#     total_drag_bottom = 0
#     total_lift_bottom = 0
#
#     # Calculate sound speed, velocity, and density
#     sound_speed = calc_sound_speed(gamma, R, T)
#     velocity = calc_velocity(mach, sound_speed)
#     rho = calc_density(dynamic_pressure, velocity)
#
#     print()
#     print('*** Initial Conditions ***' * 3)
#     print()
#     print(f"At {T} Kelvin, Calculated Sound Speed: {sound_speed:.2f} m/s")
#     print(f"At Mach {mach}, Calculated Velocity: {velocity:.2f} m/s")
#     print(f"Therefore, the Air Density is: {rho:.5f} kg/m^3")
#     print(f"At {dynamic_pressure} Pa, dynamic pressure on the vehicle.")
#
#     # Rotate the flow direction based on the angle of attack
#     adjusted_velocity = rotate_flow(np.array([velocity, 0, 0]), angle_of_attack_degrees)
#     print(f"The {angle_of_attack_degrees} degree angle of attack causes a rotated Velocity Components: {adjusted_velocity} around the Z-axis.")
#
#     # Loop through each STL file and compute the forces on each panel
#     for file_path in file_paths:
#         normals, vertices = parse_stl(file_path)
#         for normal, vertex_set in zip(normals, vertices):
#             drag, lift, sin_theta, cp, pressure = calc_force_for_triangle(normal, vertex_set, adjusted_velocity, dynamic_pressure, Free_stream_pressure)
#             if 'top' in file_path:
#                 total_drag_top += drag
#                 total_lift_top += lift  # Separating top and bottom files for analysis of Q2c
#             else:
#                 total_drag_bottom += drag
#                 total_lift_bottom += lift
#             print(f"Normal: {normal}, Area: {calc_area(*vertex_set):.4f}, Sin(theta): {sin_theta:.4f}, Cp: {cp:.4f}, Pressure: {pressure:.2f} Pa, Drag: {drag:.2f}, Lift: {lift:.2f}")
#
#     print()
#     print('*** Results ***' * 5)
#     print()
#     print(f"Top Panel - Total Drag Force: {total_drag_top:.2f} N, Total Lift Force: {total_lift_top:.2f} N, Lift/Drag Ratio: {total_lift_top / total_drag_top:.2f}")
#     print(f"Bottom Panel - Total Drag Force: {total_drag_bottom:.2f} N, Total Lift Force: {total_lift_bottom:.2f} N, Lift/Drag Ratio: {total_lift_bottom / total_drag_bottom:.2f}")
#
#     # Calculate combined drag and lift forces
#     combined_drag = total_drag_top + total_drag_bottom
#     combined_lift = total_lift_bottom - total_lift_top
#
#     print()
#     print(f"Combined - Total Drag Force: {combined_drag:.2f} N, Total Lift Force: {combined_lift:.2f} N, Lift/Drag Ratio: {combined_lift / combined_drag:.2f}")
#
# # Run the program to compute lift and drag
# compute_aerodynamic_forces(file_paths, angle_of_attack_degrees)
#
# # Sound speed, velocity, and density calculations
# sound_speed = calc_sound_speed(gamma, R, T)
# velocity = calc_velocity(mach, sound_speed)
# rho = calc_density(dynamic_pressure, velocity)
# print(f"\nDensity: {rho:.4f} kg/m^3")
