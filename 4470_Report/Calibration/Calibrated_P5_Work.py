# #Calibration
# Average voltage for P5 Experimental


import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Load the data from an LVM file using numpy
kpa15_filename = 'kpa15_data.lvm'
X_Value, p_st1, p_st2, p_plate, TC1 = np.loadtxt(kpa15_filename, delimiter=',', skiprows=24, usecols=(0, 1, 2, 3, 4), unpack=True)

# Define the time range you are interested in
time_range = (0.019645, 0.01980)  # Adjust this to the time range you are interested in

# Filter the data to get only the values within the specified time range
filtered_indices = (X_Value >= time_range[0]) & (X_Value <= time_range[1])
# filtered_p_st1 = p_st1[filtered_indices]
filtered_p_st2 = p_st2[filtered_indices]

# Calculate the averages of the filtered values
# average_p_st1 = np.mean(filtered_p_st1) if filtered_p_st1.size > 0 else float('nan')
average_p_st2 = np.mean(filtered_p_st2) if filtered_p_st2.size > 0 else float('nan')

# Output the results
# print(f"Average voltage for p_st1 within the time range {time_range} is: {average_p_st1:.12f}")
print(f"Average voltage for p_st2 within the time range {time_range} is: {average_p_st2:.12f}")

# Plotting for visualization
plt.figure(figsize=(10, 5))
plt.xlim(left=0.019350, right=0.01965)
plt.ylim(bottom=-0.05, top=0.4)
# plt.plot(X_Value, p_st1, label='Pressure Sensor 1 (Volts)')
plt.plot(X_Value, p_st2, label='Pressure Sensor 2 (Volts)')
plt.axvline(x=time_range[0], color='green', linestyle='--', label='Start of Time Range')
plt.axvline(x=time_range[1], color='red', linestyle='--', label='End of Time Range')
plt.xlabel('Time (s)')
plt.ylabel('Voltage (Volts)')
plt.title('Pressure 5 Calibration Zone 15 kPa')
plt.legend()
plt.show()



# Plotting zoomed in P5 at top


# Path to the .lvm file
kpa15_filename = 'kpa15_data.lvm'

# Load the data
X_Value, p_st1, p_st2, p_plate, TC1 = np.loadtxt(kpa15_filename, delimiter=',', skiprows=24, usecols=(0, 1, 2, 3, 4), unpack=True)

# Interpolation functions for p_st1 and p_st2
interp_p_st1 = interp1d(X_Value, p_st1, kind='linear', fill_value="extrapolate")
interp_p_st2 = interp1d(X_Value, p_st2, kind='linear', fill_value="extrapolate")

def get_pressure_values_at_x(x_location):
    """Return pressure values at a specific x location."""
    p_st1_at_x = interp_p_st1(x_location)
    p_st2_at_x = interp_p_st2(x_location)
    return p_st1_at_x, p_st2_at_x

# Example usage
x_location = 0.0197  # Example x value
p_st1_value, p_st2_value = get_pressure_values_at_x(x_location)
print(f'Pressure Sensor 1 at x={x_location}: {p_st1_value:.3f} Volts')
print(f'Pressure Sensor 2 at x={x_location}: {p_st2_value:.3f} Volts')



# Plotting
plt.figure()
# plt.plot(X_Value, p_st1, label='Pressure Sensor 1 (Volts)')
plt.plot(X_Value, p_st2, color= 'orange',label='Pressure Sensor 2 (Volts)')

plt.xlabel('Time (s)', fontsize=10)
plt.ylabel('Pressure (Volts)', fontsize=10)
plt.title('Shock velocity determination', weight='bold', fontsize=10)
plt.legend(fontsize=8)
plt.xlim(left=0.0196, right=0.02)
# plt.ylim(bottom=-0.05, top=0.7)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.subplots_adjust(top=0.943, bottom=0.095, left=0.13, right=0.945, hspace=0.200, wspace=0.200)
# plt.savefig("stag finder2.png", format='png', dpi=500)
plt.show()


# Average voltage for P5 Experimental
