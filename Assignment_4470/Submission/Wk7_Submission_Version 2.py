# import numpy as np
# from stl import mesh
#
# # Constants
# cp_max = 1.8393710511306676
# dynamic_pressure = 50000  # Pa
# mach = 8
# gamma = 1.4
# R = 287
# T = 220
# Free_stream_pressure = 1116  # Pa
# angle_of_attack_degrees = 0 # degrees
#
# def calc_sound_speed(gamma, R, T):
#     """Calculate the speed of sound."""
#     return np.sqrt(gamma * R * T)
#
# def calc_velocity(mach, sound_speed):
#     """Calculate the velocity."""
#     return mach * sound_speed
#
# def calc_density(dynamic_pressure, velocity):
#     """Calculate the density from dynamic pressure and velocity."""
#     return 2 * dynamic_pressure / (velocity ** 2)
#
# def load_stl(file_path):
#     return mesh.Mesh.from_file(file_path)
#
# def normalize_vector(vector):
#     return vector / np.linalg.norm(vector)
#
# def calc_area(v0, v1, v2):
#     return 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))
#
# def calc_sin_theta(velocity, normal):
#     return np.dot(normalize_vector(velocity), normal)
#
# def rotate_flow(flow_dir, angle_degrees):
#     angle_radians = np.radians(angle_degrees)
#     rot_matrix = np.array([
#         [np.cos(angle_radians), -np.sin(angle_radians), 0],
#         [np.sin(angle_radians), np.cos(angle_radians), 0],
#         [0, 0, 1]
#     ])
#     return np.dot(rot_matrix, flow_dir)
#
# def calc_force_for_triangle(normal, vertices, velocity, dynamic_pressure, static_pressure):
#     normal = normalize_vector(normal)
#     area = calc_area(*vertices)
#     sin_theta = calc_sin_theta(velocity, normal)
#     if sin_theta < 0:
#         cp = cp_max * (sin_theta**2)
#         pressure = dynamic_pressure
#     else:
#         cp = cp_max * (sin_theta**2)  # Apply static pressure in the shadow region
#         pressure = static_pressure
#     force = -pressure * cp * area * normal
#     velocity_normalized = normalize_vector(velocity)
#     drag = np.dot(force, velocity_normalized)
#     lift = np.dot(force, np.array([-velocity_normalized[1], velocity_normalized[0], 0]))
#     return drag, lift, sin_theta, cp, pressure
#
# def compute_aerodynamic_forces(file_path, angle_of_attack_degrees):
#     total_drag = 0.0
#     total_lift = 0.0
#     sound_speed = calc_sound_speed(gamma, R, T)
#     velocity = calc_velocity(mach, sound_speed)
#     print(f"Calculated Sound Speed: {sound_speed:.2f} m/s")
#     print(f"Calculated Velocity: {velocity:.2f} m/s")
#     adjusted_velocity = rotate_flow(np.array([velocity, 0, 0]), angle_of_attack_degrees)
#     print(f"Rotated Velocity Vector: {adjusted_velocity}")
#     mesh_data = load_stl(file_path)
#     for normal, vertices in zip(mesh_data.normals, mesh_data.vectors):
#         drag, lift, sin_theta, cp, pressure = calc_force_for_triangle(normal, vertices, adjusted_velocity, dynamic_pressure, Free_stream_pressure)
#         total_drag += drag
#         total_lift += lift
#         print(normal)
#         print(f"Normal: {normal}, Area: {calc_area(*vertices):.4f}, Sin(theta): {sin_theta:.4f}, Cp: {cp:.4f}, Pressure: {pressure:.2f} Pa, Drag: {drag:.2f}, Lift: {lift:.2f}")
#     print(f"Total Drag Force: {total_drag:.2f} N, Total Lift Force: {total_lift:.2f} N, Lift-to-Drag Ratio: {total_lift/total_drag:.2f}")
#
# # Main computation
# file_path = 'wk-07.stl'
# compute_aerodynamic_forces(file_path, angle_of_attack_degrees)
#
# # Sound speed, velocity, and density calculations
# def density(q, U):
#     return 2 * q / U**2
#
# sound_speed = calc_sound_speed(gamma, R, T)
# velocity = calc_velocity(mach, sound_speed)
# rho = density(dynamic_pressure, velocity)
# print(f"\nDensity: {rho:.4f} kg/m^3")
#
# print()
# # Function to load the STL file and print normals
# def load_stl_and_print_normals(file_path):
#     """Load an STL file and print the normals."""
#     mesh_data = mesh.Mesh.from_file(file_path)
#     for normal in mesh_data.normals:
#         print(f"Normal: {normal}")
#
# Code file execution
# file_path = 'wk-07.stl'
# load_stl_and_print_normals(file_path)



#
# import numpy as np
# from stl import mesh
#
# # Constants
# cp_max = 1.8393710511306676
# dynamic_pressure = 50000  # Pa
# mach = 8
# gamma = 1.4
# R = 287
# T = 220
# Free_stream_pressure = 5  # Pa
# angle_of_attack_degrees = 0  # degrees
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
# # Function to load the STL file and return the normals and vertices
# def load_stl(file_path):
#     """Load an STL file and return the normals and vertices."""
#     stl_mesh = mesh.Mesh.from_file(file_path)
#     normals = np.copy(stl_mesh.normals)
#     vertices = np.copy(stl_mesh.vectors)
#     return normals, vertices
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
#     lift = np.dot(force, np.array([-velocity_normalized[1], velocity_normalized[0], 0]))
#
#     return drag, lift, sin_theta, cp, pressure
#
# # Function to compute the aerodynamic forces for a given STL file
# def compute_aerodynamic_forces(file_path, angle_of_attack_degrees):
#     """Compute the total drag and lift forces for the given STL file."""
#     total_drag = 0.0
#     total_lift = 0.0
#     sound_speed = calc_sound_speed(gamma, R, T)
#     velocity = calc_velocity(mach, sound_speed)
#     print(f"Calculated Sound Speed: {sound_speed:.2f} m/s")
#     print(f"Calculated Velocity: {velocity:.2f} m/s")
#     adjusted_velocity = rotate_flow(np.array([velocity, 0, 0]), angle_of_attack_degrees)
#     print(f"Rotated Velocity Vector: {adjusted_velocity}")
#     normals, vertices = load_stl(file_path)
#     for normal, vertex_set in zip(normals, vertices):
#         drag, lift, sin_theta, cp, pressure = calc_force_for_triangle(normal, vertex_set, adjusted_velocity, dynamic_pressure, Free_stream_pressure)
#         total_drag += drag
#         total_lift += lift
#         print(f"Normal: {normal}, Area: {calc_area(*vertex_set):.4f}, Sin(theta): {sin_theta:.4f}, Cp: {cp:.4f}, Pressure: {pressure:.2f} Pa, Drag: {drag:.2f}, Lift: {lift:.2f}")
#     print(f"Total Drag Force: {total_drag:.2f} N, Total Lift Force: {total_lift:.2f} N, Lift-to-Drag Ratio: {total_lift/total_drag:.2f}")
#
# # Main computation
# file_path = 'wk-07.stl'
# compute_aerodynamic_forces(file_path, angle_of_attack_degrees)
#
# # Sound speed, velocity, and density calculations
# def density(q, U):
#     """Calculate the density from dynamic pressure and velocity."""
#     return 2 * q / U ** 2
#
# sound_speed = calc_sound_speed(gamma, R, T)
# velocity = calc_velocity(mach, sound_speed)
# rho = density(dynamic_pressure, velocity)
# print(f"\nDensity: {rho:.4f} kg/m^3")
#
# print()
# # Function to load the STL file and print normals
# def load_stl_and_print_normals(file_path):
#     """Load an STL file and print the normals."""
#     stl_mesh = mesh.Mesh.from_file(file_path)
#     for normal in stl_mesh.normals:
#         print(f"Normal: {normal}")
#
# # Main computation
# file_path = 'wk-07.stl'
# load_stl_and_print_normals(file_path)
#
#
# import numpy as np
#
# # Constants
# cp_max = 1.8393710511306676
# dynamic_pressure = 50000  # Pa
# mach = 8
# gamma = 1.4
# R = 287
# T = 220
# Free_stream_pressure = 5  # Pa
# angle_of_attack_degrees = 0  # degrees
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
#         pressure = static_pressure # is this causing an issue should this be smaller so -pressure * cp * area * normal is a smaller number or something??>>?>>>?>???
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
#     lift = np.dot(force, np.array([velocity_normalized[1], velocity_normalized[0], 0]))
#
#     return drag, lift, sin_theta, cp, pressure
#
# # Function to compute the aerodynamic forces for given facet normals and vertices
# def compute_aerodynamic_forces(normals, vertices, angle_of_attack_degrees):
#     """Compute the total drag and lift forces for the given facet normals and vertices."""
#     total_drag = 0.0
#     total_lift = 0.0
#     sound_speed = calc_sound_speed(gamma, R, T)
#     velocity = calc_velocity(mach, sound_speed)
#     print(f"Calculated Sound Speed: {sound_speed:.2f} m/s")
#     print(f"Calculated Velocity: {velocity:.2f} m/s")
#     adjusted_velocity = rotate_flow(np.array([velocity, 0, 0]), angle_of_attack_degrees)
#     print(f"Rotated Velocity Vector: {adjusted_velocity}")
#     for normal, vertex_set in zip(normals, vertices):
#         drag, lift, sin_theta, cp, pressure = calc_force_for_triangle(normal, vertex_set, adjusted_velocity, dynamic_pressure, Free_stream_pressure)
#         total_drag += drag
#         total_lift += lift
#         print(f"Normal: {normal}, Area: {calc_area(*vertex_set):.4f}, Sin(theta): {sin_theta:.4f}, Cp: {cp:.4f}, Pressure: {pressure:.2f} Pa, Drag: {drag:.2f}, Lift: {lift:.2f}")
#     print(f"Total Drag Force: {total_drag:.2f} N, Total Lift Force: {total_lift:.2f} N, Lift-to-Drag Ratio: {total_lift/total_drag:.2f}")
#
# # Facet normals and vertices from the provided text file
# normals = [
#     [-1.000000, 0.000000, 0.000000],
#     [-1.000000, -0.000000, 0.000000],
#     [-0.156434, -0.000000, -0.987688],
#     [-0.156434, 0.000000, -0.987688],
#     [-0.156434, 0.987688, 0.000000],
#     [-0.156434, 0.987688, 0.000000],
#     [-0.156434, 0.000000, 0.987688],
#     [-0.156434, 0.000000, 0.987688],
#     [-0.156434, -0.987688, 0.000000],
#     [-0.156434, -0.987688, 0.000000]
# ]
#
# vertices = [
#     [[0.000000, -0.250000, -0.250000], [0.000000, 0.250000, 0.250000], [0.000000, 0.250000, -0.250000]],
#     [[0.000000, -0.250000, -0.250000], [0.000000, -0.250000, 0.250000], [0.000000, 0.250000, 0.250000]],
#     [[1.262750, -0.450000, -0.450000], [0.000000, -0.250000, -0.250000], [1.262750, 0.450000, -0.450000]],
#     [[0.000000, -0.250000, -0.250000], [0.000000, 0.250000, -0.250000], [1.262750, 0.450000, -0.450000]],
#     [[0.000000, 0.250000, -0.250000], [0.000000, 0.250000, 0.250000], [1.262750, 0.450000, -0.450000]],
#     [[0.000000, 0.250000, 0.250000], [1.262750, 0.450000, 0.450000], [1.262750, 0.450000, -0.450000]],
#     [[0.000000, -0.250000, 0.250000], [1.262750, 0.450000, 0.450000], [0.000000, 0.250000, 0.250000]],
#     [[0.000000, -0.250000, 0.250000], [1.262750, -0.450000, 0.450000], [1.262750, 0.450000, 0.450000]],
#     [[1.262750, -0.450000, -0.450000], [1.262750, -0.450000, 0.450000], [0.000000, -0.250000, 0.250000]],
#     [[1.262750, -0.450000, -0.450000], [0.000000, -0.250000, 0.250000], [0.000000, -0.250000, -0.250000]]
# ]
#
# # Main computation
# compute_aerodynamic_forces(normals, vertices, angle_of_attack_degrees)
#
# # Sound speed, velocity, and density calculations
# def density(q, U):
#     """Calculate the density from dynamic pressure and velocity."""
#     return 2 * q / U ** 2
#
# sound_speed = calc_sound_speed(gamma, R, T)
# velocity = calc_velocity(mach, sound_speed)
# rho = density(dynamic_pressure, velocity)
# print(f"\nDensity: {rho:.4f} kg/m^3")


import numpy as np

# Constants
cp_max = 1.8393710511306676
dynamic_pressure = 50000  # Pa
mach = 8
gamma = 1.4
R = 287
T = 220
Free_stream_pressure = 5  # Pa
angle_of_attack_degrees = 9  # degrees


def calc_sound_speed(gamma, R, T):
    """Calculate the speed of sound."""
    return np.sqrt(gamma * R * T)


def calc_velocity(mach, sound_speed):
    """Calculate the velocity."""
    return mach * sound_speed


def calc_density(dynamic_pressure, velocity):
    """Calculate the density from dynamic pressure and velocity."""
    return 2 * dynamic_pressure / (velocity ** 2)


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


def normalize_vector(vector):
    return vector / np.linalg.norm(vector)


def calc_area(v0, v1, v2):
    return 0.5 * np.linalg.norm(np.cross(np.array(v1) - np.array(v0), np.array(v2) - np.array(v0)))


def calc_sin_theta(velocity, normal):
    return np.dot(normalize_vector(velocity), normalize_vector(normal))


def rotate_flow(flow_dir, angle_degrees):
    angle_radians = np.radians(angle_degrees)
    rot_matrix = np.array([
        [np.cos(angle_radians), -np.sin(angle_radians), 0],
        [np.sin(angle_radians), np.cos(angle_radians), 0],
        [0, 0, 1]
    ])
    return np.dot(rot_matrix, flow_dir)


def calc_force_for_triangle(normal, vertices, velocity, dynamic_pressure, static_pressure):
    normal = normalize_vector(normal)
    area = calc_area(*vertices)
    sin_theta = calc_sin_theta(velocity, normal)
    if sin_theta < 0:
        cp = cp_max * (sin_theta ** 2)
        pressure = dynamic_pressure
    else:
        cp = cp_max * (sin_theta ** 2)
        pressure = static_pressure
    force = -pressure * cp * area * normal
    velocity_normalized = normalize_vector(velocity)
    drag = np.dot(force, velocity_normalized)
    lift = np.dot(force, np.array([-velocity_normalized[1], velocity_normalized[0], 0]))
    return drag, lift, sin_theta, cp, pressure


def compute_aerodynamic_forces(file_path, angle_of_attack_degrees):
    total_drag = 0.0
    total_lift = 0.0
    sound_speed = calc_sound_speed(gamma, R, T)
    velocity = calc_velocity(mach, sound_speed)
    adjusted_velocity = rotate_flow(np.array([velocity, 0, 0]), angle_of_attack_degrees)
    normals, vertices = parse_stl(file_path)
    for i, (normal, vertex_set) in enumerate(zip(normals, vertices)):
        drag, lift, sin_theta, cp, pressure = calc_force_for_triangle(normal, vertex_set, adjusted_velocity,
                                                                      dynamic_pressure, Free_stream_pressure)
        total_drag += drag
        total_lift += lift
        print(
            f"Normal {i}: {normal}, Area: {calc_area(*vertex_set):.4f}, Sin(theta): {sin_theta:.4f}, Cp: {cp:.4f}, Pressure: {pressure:.2f} Pa, Drag: {drag:.2f}, Lift: {lift:.2f}")
    print(
        f"Total Drag Force: {total_drag:.2f} N, Total Lift Force: {total_lift:.2f} N, Lift-to-Drag Ratio: {total_lift / total_drag:.2f}")


# Main computation
file_path = 'wk-07.stl'
compute_aerodynamic_forces(file_path, angle_of_attack_degrees)

# Sound speed, velocity, and density calculations
sound_speed = calc_sound_speed(gamma, R, T)
velocity = calc_velocity(mach, sound_speed)
rho = calc_density(dynamic_pressure, velocity)
print(f"\nDensity: {rho:.4f} kg/m^3")
