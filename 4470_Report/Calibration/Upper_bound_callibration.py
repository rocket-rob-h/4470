import numpy as np
import matplotlib.pyplot as plt

# Load the data from an LVM file using numpy
kpa15_filename = 'kpa15_data.lvm'
X_Value, p_st1, p_st2, p_plate, TC1 = np.loadtxt(kpa15_filename, delimiter=',', skiprows=24, usecols=(0, 1, 2, 3, 4), unpack=True)

# Define the time range you are interested in
time_range = (0.019370, 0.01959)  # Adjust this to the time range you are interested in

# Filter the data to get only the values within the specified time range
filtered_indices = (X_Value >= time_range[0]) & (X_Value <= time_range[1])
filtered_p_st1 = p_st1[filtered_indices]
filtered_p_st2 = p_st2[filtered_indices]

# Calculate the averages of the filtered values
average_p_st1 = np.mean(filtered_p_st1) if filtered_p_st1.size > 0 else float('nan')
average_p_st2 = np.mean(filtered_p_st2) if filtered_p_st2.size > 0 else float('nan')

# Output the results
print(f"Average voltage for p_st1 within the time range {time_range} is: {average_p_st1:.12f}")
print(f"Average voltage for p_st2 within the time range {time_range} is: {average_p_st2:.12f}")

# Plotting for visualization
plt.figure(figsize=(10, 5))
plt.xlim(left=0.019350, right=0.01965)
plt.ylim(bottom=-0.05, top=0.4)
plt.plot(X_Value, p_st1, label='Pressure Sensor 1 (Volts)')
plt.plot(X_Value, p_st2, label='Pressure Sensor 2 (Volts)')
plt.axvline(x=time_range[0], color='green', linestyle='--', label='Start of Time Range')
plt.axvline(x=time_range[1], color='red', linestyle='--', label='End of Time Range')
plt.xlabel('Time (s)')
plt.ylabel('Voltage (Volts)')
plt.title('Pressure Two Calibration Zone 15 kPa')
plt.legend()
plt.show()
