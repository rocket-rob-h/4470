import numpy as np
import matplotlib.pyplot as plt

# Path to the .lvm file
kpa15_filename = 'kpa15_data.lvm'

# Load the data
X_Value, p_st1, p_st2, p_plate, TC1 = np.loadtxt(kpa15_filename, delimiter=',', skiprows=24, usecols=(0, 1, 2, 3, 4), unpack=True)

# Calculate the absolute derivative (approximate)
# delta_p_st1 = np.abs(np.diff(p_st1))  # Absolute change in p_st1
# delta_p_st2 = np.abs(np.diff(p_st2))  # Absolute change in p_st2

# # Define a threshold for what you consider a 'large change'
# threshold = 0.08  # This threshold needs tuning based on your data characteristics
#
# # Find indices where the change exceeds the threshold
# significant_changes_st1 = np.where(delta_p_st1 > threshold)[0]
# significant_changes_st2 = np.where(delta_p_st2 > threshold)[0]
#
# # Output the times at which significant changes occur for both sensors
# times_with_significant_change_st1 = X_Value[significant_changes_st1]
# times_with_significant_change_st2 = X_Value[significant_changes_st2]

# # Print the results
# print("Times with significant changes in Sensor 1:", times_with_significant_change_st1)
# print("Times with significant changes in Sensor 2:", times_with_significant_change_st2)


# Plotting
plt.figure()
plt.plot(X_Value, p_plate, label='Plate Sensor (Volts)')
# Add a horizontal line at y = 0
plt.axhline(0, color='black', linewidth=2, label='Resting')  # You can adjust the color and linewidth as needed

# plt.scatter(X_Value[significant_changes_st1], p_st1[significant_changes_st1], color='red', label='Change Sensor 1')
# plt.scatter(X_Value[significant_changes_st2], p_st2[significant_changes_st2], color='blue', label='Change Sensor 2')

# # Annotate each point with its corresponding time value
# for i in significant_changes_st1:
#     plt.annotate(f'{X_Value[i]:.6f}', (X_Value[i], p_st1[i]), textcoords="offset points", xytext=(35,12), ha='center', fontsize=12)
#
# for i in significant_changes_st2:
#     plt.annotate(f'{X_Value[i]:.6f}', (X_Value[i], p_st2[i]), textcoords="offset points", xytext=(35,12), ha='center', fontsize=12)


#Plotting
plt.xlabel('Time (s)', fontsize=10)
plt.ylabel('Pressure (Volts)', fontsize=10)
plt.title('Plate Pressure', weight='bold', fontsize=10)
plt.legend(fontsize=8)


# Set x and y limits to focus on specific intervals
# plt.xlim(left=0.01915, right=0.0194)
# plt.ylim(bottom=-0.05, top=0.4)

# Increase the fontsize of the axis number values
plt.tick_params(axis='both', which='major', labelsize=8)
# Adjust the padding and the spacing to match the screenshot
plt.subplots_adjust(top=0.943, bottom=0.095, left=0.13, right=0.945, hspace=0.200, wspace=0.200)

# Save and show the plot

plt.savefig("Plate Pressure.png", format='png', dpi=500)
plt.show()