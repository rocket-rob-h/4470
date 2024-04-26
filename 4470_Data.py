# import numpy as np
# import matplotlib.pyplot as plt
#
# # Path to the .lvm file
# kpa15_filename = 'kpa15_data.lvm'
#
#
# # Load the data, specifying columns to read (skipping the 'Comment' column)
# # Indices: 0 = X_Value, 1 = p_st1, 2 = p_st2, 3 = p_plate, 4 = TC1
# X_Value, p_st1, p_st2, p_plate, TC1 = np.loadtxt(kpa15_filename, delimiter=',', skiprows=24, usecols=(0, 1, 2, 3, 4), unpack=True)
#
#
# # plt.xlim(X_Value[100], X_Value[200])
#
# # Plotting example: plot first sensor pressure data vs time
# plt.figure()
# plt.plot(X_Value, p_st1, label='Pressure Sensor 1 (Volts)')
# plt.plot(X_Value, p_st2, label = 'Pressure Sensor 2(Volts)')
# plt.xlabel('Time (s)')
# plt.ylabel('Pressure (Volts)')
# plt.title('Pressure 1+2 Shock speed zoom')
# plt.legend()
#
# plt.xlim(left=0.01915, right=0.0201)
#
# plt.savefig("Pressure 1+2 Shock speed zoom", format='png')
# plt.show()
#
#
#

import numpy as np
import matplotlib.pyplot as plt

# Path to the .lvm file
kpa15_filename = 'kpa15_data.lvm'

# Load the data
X_Value, p_st1, p_st2, p_plate, TC1 = np.loadtxt(kpa15_filename, delimiter=',', skiprows=24, usecols=(0, 1, 2, 3, 4), unpack=True)

# Calculate the absolute derivative (approximate)
delta_p_st1 = np.abs(np.diff(p_st1))  # Absolute change in p_st1
delta_p_st2 = np.abs(np.diff(p_st2))  # Absolute change in p_st2

# Define a threshold for what you consider a 'large change'
threshold = 0.08  # This threshold needs tuning based on your data characteristics

# Find indices where the change exceeds the threshold
significant_changes_st1 = np.where(delta_p_st1 > threshold)[0]
significant_changes_st2 = np.where(delta_p_st2 > threshold)[0]

# Output the times at which significant changes occur for both sensors
times_with_significant_change_st1 = X_Value[significant_changes_st1]
times_with_significant_change_st2 = X_Value[significant_changes_st2]

# Print the results
print("Times with significant changes in Sensor 1:", times_with_significant_change_st1)
print("Times with significant changes in Sensor 2:", times_with_significant_change_st2)


# Plotting
plt.figure()
plt.plot(X_Value, p_st1, label='Pressure Sensor 1 (Volts)')
plt.plot(X_Value, p_st2, label='Pressure Sensor 2 (Volts)')
plt.scatter(X_Value[significant_changes_st1], p_st1[significant_changes_st1], color='red', label='Significant Change Sensor 1')
plt.scatter(X_Value[significant_changes_st2], p_st2[significant_changes_st2], color='blue', label='Significant Change Sensor 2')
plt.xlabel('Time (s)')
plt.ylabel('Pressure (Volts)')
plt.title('Pressure 1+2 Shock speed zoom')
plt.legend()

# Set x-limits to focus on a specific time interval
plt.xlim(left=0.01915, right=0.0201)

# Save and show the plot
plt.savefig("Pressure 1+2 Shock speed zoom.png", format='png')
plt.show()



