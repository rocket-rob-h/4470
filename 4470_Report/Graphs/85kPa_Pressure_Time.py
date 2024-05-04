
import numpy as np
import matplotlib.pyplot as plt

# Path to the .lvm file
kpa15_filename = 'kpa85_data.lvm'

# Load the data
X_Value, p_st1, p_st2, p_plate, TC1 = np.loadtxt(kpa15_filename, delimiter=',', skiprows=24, usecols=(0, 1, 2, 3, 4), unpack=True)

# Resting pressure P1 in kPa
P1 = 85  # Resting pressure is already in kPa

# Calibration functions
def Calib85_1(V):
    return P1 + 1340.381407 * V  # Calibration for p_st1

def Calib85_2(V):
    return P1 + 1277.449 * V     # Calibration for p_st2

# Apply calibration to convert voltage to pressure in kPa
pressure_st1 = Calib85_1(p_st1)
pressure_st2 = Calib85_2(p_st2)

# Plotting the calibrated pressure
plt.figure()
plt.plot(X_Value, pressure_st1, label='Pressure Sensor 1 (kPa)')
plt.plot(X_Value, pressure_st2, label='Pressure Sensor 2 (kPa)')

plt.xlabel('Time (s)', fontsize=10)
plt.ylabel('Pressure (kPa)', fontsize=10)
plt.title('Calibrated Pressure Distribution at 85 kPa', weight='bold', fontsize=10)
plt.legend(fontsize=8)

# Set x and y limits to focus on specific intervals
plt.xlim(left=0.017, right=0.035)
# plt.ylim(bottom=-0.05, top=0.7)

# Increase the fontsize of the axis number values
plt.tick_params(axis='both', which='major', labelsize=8)
# Adjust the padding and the spacing to match the screenshot
plt.subplots_adjust(top=0.943, bottom=0.095, left=0.13, right=0.945, hspace=0.200, wspace=0.200)


plt.savefig("85kPa_Pressure_Time.png", format='png', dpi=500)
plt.show()
