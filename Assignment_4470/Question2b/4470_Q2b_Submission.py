
import numpy as np
from stl import mesh

# Constants
dynamic_pressure = 50000  # Pa
mach = 8
gamma = 1.4
R = 287
T = 220
Free_stream_pressure = 10  # Pa
angle_of_attack_degrees = 7.52 # degrees
file_paths = ['top-s47661808.stl', 'bot-s47661808.stl']

def CP_MAX(gamma):
    """Calculate maximum coefficient of pressure based on gamma."""
    return (4 / (gamma + 1)) * ((gamma + 1)**2 / (4 * gamma))**(gamma / (gamma - 1))

cp_max = CP_MAX(gamma)

def calc_sound_speed(gamma, R, T):
    """Calculate the speed of sound."""
    return np.sqrt(gamma * R * T)

def calc_velocity(mach, sound_speed):
    """Calculate the velocity."""
    return mach * sound_speed

def calc_density(dynamic_pressure, velocity):
    """Calculate the density from dynamic pressure and velocity."""
    return 2 * dynamic_pressure / (velocity ** 2)

def load_stl(file_path):
    return mesh.Mesh.from_file(file_path)

def normalize_vector(vector):
    return vector / np.linalg.norm(vector)

def calc_area(v0, v1, v2):
    return 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))

def calc_sin_theta(velocity, normal):
    return np.dot(normalize_vector(velocity), normal)

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
        cp = cp_max * (sin_theta**2)
        pressure = dynamic_pressure  # apply dynamic pressure to panels in flow field
    else:
        cp = cp_max * (sin_theta**2)  # Apply static pressure in panels in shadow to balance forces
        pressure = static_pressure
    force = -pressure * cp * area * normal  # negative pressure as per task sheet
    velocity_normalized = normalize_vector(velocity)
    drag = np.dot(force, velocity_normalized)
    lift = np.dot(force, np.array([velocity_normalized[1], velocity_normalized[0], 0]))  # check this if velocity_normalized[1] is a negative or not to account for top panel lift
    return drag, lift, sin_theta, cp, pressure

def compute_aerodynamic_forces(file_paths, angle_of_attack_degrees):
    '''Set Initial Lift and Drag to 0 before loop'''
    total_drag_top = 0
    total_lift_top = 0
    total_drag_bottom = 0
    total_lift_bottom = 0
    sound_speed = calc_sound_speed(gamma, R, T)
    velocity = calc_velocity(mach, sound_speed)
    rho = calc_density(dynamic_pressure, velocity)
    print()
    print('*** Initial Conditions ***' * 3)
    print()
    print(f"At {T} in Kelvin, Calculated Sound Speed: {sound_speed:.2f} m/s")
    print(f"At Mach {mach}, Calculated Velocity: {velocity:.2f} m/s")
    rho = calc_density(dynamic_pressure, velocity)
    print(f"Therefore the Air Density is : {rho:.5f} kg/m^3")
    print(f'At {dynamic_pressure} Pa, dynamic pressure on the vehicle. ')
    adjusted_velocity = rotate_flow(np.array([velocity, 0, 0]), angle_of_attack_degrees)
    print(f"The {angle_of_attack_degrees} degree angle of attack causes a rotated Velocity Components : {adjusted_velocity} around Z axis.")

    for file_path in file_paths:
        mesh_data = load_stl(file_path)
        for normal, vertices in zip(mesh_data.normals, mesh_data.vectors):
            drag, lift, sin_theta, cp, pressure = calc_force_for_triangle(normal, vertices, adjusted_velocity, dynamic_pressure, Free_stream_pressure)
            if 'top' in file_path:
                total_drag_top += drag
                total_lift_top += lift # Sperating top and bottom files for analysis of Q2c
            else:
                total_drag_bottom += drag
                total_lift_bottom += lift
            # print(f" Sin(theta): {sin_theta:.4f}, Cp: {cp:.4f}, Pressure: {pressure:.2f} Pa, Drag: {drag:.2f}, Lift: {lift:.2f}")

    # Sound speed, velocity, and density calculations
    print()
    print('*** Results ***'*5)
    print()
    print(f"Top Panel - Total Drag Force: {total_drag_top:.2f} N, Total Lift Force: {total_lift_top:.2f} N, Lift/Drag Ratio: {total_lift_top/total_drag_top:.2f}")
    print(f"Bottom Panel - Total Drag Force: {total_drag_bottom:.2f} N, Total Lift Force: {total_lift_bottom:.2f} N, Lift/Drag Ratio: {total_lift_bottom/total_drag_bottom:.2f}")
    combined_drag = total_drag_top + total_drag_bottom
    combined_lift = total_lift_top + total_lift_bottom
    print()
    print(f"Combined - Total Drag Force: {combined_drag:.2f} N, Total Lift Force: {combined_lift:.2f} N, Lift/Drag Ratio: {combined_lift/combined_drag:.2f}")


# Run Lift/Drag Program
compute_aerodynamic_forces(file_paths, angle_of_attack_degrees)