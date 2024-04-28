import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

# Load the data from an LVM file using numpy
kpa15_filename = 'kpa15_data.lvm'
X_Value, p_st1, p_st2, p_plate, TC1 = np.loadtxt(kpa15_filename, delimiter=',', skiprows=24, usecols=(0, 1, 2, 3, 4), unpack=True)

# # Calculate the absolute derivative (approximate) to identify changes
# delta_p_st1 = np.abs(np.diff(p_st1))  # Absolute change in p_st1
# delta_p_st2 = np.abs(np.diff(p_st2))  # Absolute change in p_st2
#
# # Define a threshold for what you consider a 'significant change'
# threshold = 0.06  # Adjust this threshold based on your data characteristics
#
# # Find indices where the change exceeds the threshold
# significant_changes_st1 = np.where(delta_p_st1 > threshold)[0]
# significant_changes_st2 = np.where(delta_p_st2 > threshold)[0]
#
# # Assuming the first significant change is what you are interested in
# first_significant_change_st1 = significant_changes_st1[0] if significant_changes_st1.size > 0 else len(p_st1)
# first_significant_change_st2 = significant_changes_st2[0] if significant_changes_st2.size > 0 else len(p_st2)
#
# # Calculate averages from the start to the first significant change
# average_p_st1 = np.mean(p_st1[:first_significant_change_st1 + 1])  # Include the point at the change
# average_p_st2 = np.mean(p_st2[:first_significant_change_st2 + 1])  # Include the point at the change
#
# # Print the results
# print("Average p_st1 from start to first significant change:", average_p_st1)
# print("Average p_st2 from start to first significant change:", average_p_st2)
#
# # Plotting for visualization
# plt.figure()
# plt.plot(X_Value, p_st1, label='Pressure Sensor 1 (Volts)')
# plt.plot(X_Value, p_st2, label='Pressure Sensor 2 (Volts)')
# plt.scatter(X_Value[significant_changes_st1], p_st1[significant_changes_st1], color='red', label='Change Sensor 1')
# plt.scatter(X_Value[significant_changes_st2], p_st2[significant_changes_st2], color='blue', label='Change Sensor 2')
# plt.xlabel('Time (s)', fontsize=10)
# plt.ylabel('Pressure (Volts)', fontsize=10)
# plt.title('Shock velocity determination', weight='bold', fontsize=10)
# plt.legend(fontsize=8)
# plt.tick_params(axis='both', which='major', labelsize=8)
# plt.subplots_adjust(top=0.943, bottom=0.095, left=0.13, right=0.945, hspace=0.200, wspace=0.200)
#
#
# #Plotting
# plt.xlabel('Time (s)', fontsize=10)
# plt.ylabel('Pressure (Volts)', fontsize=10)
# plt.title('Shock velocity determination', weight='bold', fontsize=10)
# plt.legend(fontsize=8)
#
#
# # Set x and y limits to focus on specific intervals
# # plt.xlim(left=0.01915, right=0.0194)
# # plt.ylim(bottom=-0.05, top=0.4)
#
# # Increase the fontsize of the axis number values
# plt.tick_params(axis='both', which='major', labelsize=8)
# # Adjust the padding and the spacing to match the screenshot
# plt.subplots_adjust(top=0.943, bottom=0.095, left=0.13, right=0.945, hspace=0.200, wspace=0.200)
#
# # Save and show the plot
# #
# # plt.savefig("Calibration_15kPa.png", format='png', dpi=500)
# plt.show()




# # Symbols for voltage and pressure
# V = sp.symbols('V')
# m, b = sp.symbols('m b')
#
# # Initial and final voltage and pressure points
# V0 = sp.Rational(0.0002330144852021676)
# P0 = sp.Rational(85)
# V1 = sp.Rational(0.043)
# P1 = sp.Rational(245.59)
#
# # Calculate the slope (m) and intercept (b)
# m = (P1 - P0) / (V1 - V0)
# b = P0 - m * V0
#
# # Define the linear equation for pressure as a function of voltage
# P = m * V + b
#
# # Display the symbolic equation
# print("Symbolic equation for pressure as a function of voltage:")
# sp.pprint(P, use_unicode=True)


# Define the voltage ranges for st1 and st2
st1_voltage_range = (0.2, 0.24)
st2_voltage_range = (0.2, 0.23)

# Filter the data to get only the values within the specified voltage ranges
filtered_p_st1 = p_st1[(p_st1 >= st1_voltage_range[0]) & (p_st1 <= st1_voltage_range[1])]
filtered_p_st2 = p_st2[(p_st2 >= st2_voltage_range[0]) & (p_st2 <= st2_voltage_range[1])]

# Calculate the averages of the filtered values
average_p_st1 = np.mean(filtered_p_st1) if filtered_p_st1.size > 0 else float('nan')
average_p_st2 = np.mean(filtered_p_st2) if filtered_p_st2.size > 0 else float('nan')

# Output the results
print(f"Average value for p_st1 within the voltage range {st1_voltage_range} is: {average_p_st1:.12f}")
print(f"Average value for p_st2 within the voltage range {st2_voltage_range} is: {average_p_st2:.12f}")