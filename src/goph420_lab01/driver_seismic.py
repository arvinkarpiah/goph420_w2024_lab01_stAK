#This script does the following:

#1. Loads the velocity data from the file s_wave_data.txt.
#2. Plots the raw data and prints the value T used in the integrations.
#3. Estimates the integral in equation (14) from Lab01 using different sampling intervals and 
#   both trapezoid rule and Simpson's rules.
#4. Prints the result for each case
#5. Calculates and plots the convergence of the integral estimates using trapezoid rule and 
#   Simpson's rules on a log-log plot.
#6. Saves the two output figures in the directory figures.
#(i) Raw data_Q7
#(ii) Convergence plot_Q7

#Import dependencies and functions
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/Users/arvinkarpiah/Desktop/GOPH420/Lab_01/Repository/goph420_w2024_lab01_stAK/src')
import integration
from integration import integrate_newton

#Load seismic data
data_path = '/Users/arvinkarpiah/Desktop/GOPH420/Lab_01/Repository/goph420_w2024_lab01_stAK/data/s_wave_data.txt'
SeismicData = np.loadtxt(data_path)
Time_Raw = SeismicData[:,0]
Velocity_Raw = SeismicData[:,1]

#Analysis for an appropriate period of event, T to use
max_abs_vel = max(abs(Velocity_Raw))

# Initialize an empty list to store indices
bin_list = []

# Determine last point in the series that satisfies below condition
for i in range(len(Velocity_Raw)):
    if abs(Velocity_Raw[i]) > (0.5/100) * max_abs_vel:
        bin_list.append(i)  

last_point = max(bin_list) 

#Plot raw data
plt.figure()
plt.plot(Time_Raw,Velocity_Raw,'b-')
plt.axvline(x=Time_Raw[last_point], color='r', linestyle='-')
plt.title('S-wave velocity vs arrival time ')
plt.xlabel('Time(s)')
plt.ylabel('Velocity(mm/s)')
plt.text(Time_Raw[last_point], max(Velocity_Raw), 'T=6.77s', verticalalignment='bottom', horizontalalignment='right')
plt.savefig('/Users/arvinkarpiah/Desktop/GOPH420/Lab_01/Repository/goph420_w2024_lab01_stAK/figures/Raw data_Q7.png')
plt.savefig('/Users/arvinkarpiah/Desktop/GOPH420/Lab_01/Repository/goph420_w2024_lab01_stAK/Raw data_Q7.png')

#s will control sampling interval to be used
s = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256])

#Period in seconds
T = Time_Raw[last_point]
print('T=', T)

#Initialize empty lists to store variables
I_TrapeziumRule = np.zeros(len(s), dtype=float)  # Array to store results for Trapezium rule
I_SimpsonsRule = np.zeros(len(s), dtype=float)  # Array to store results for Simpsons rule
Sampling = np.zeros(len(s), dtype=float)  # Array to store sampling intervals
epsilon_a_trap = np.zeros(len(s)-1, dtype=float) # Array to store approximate relative error (trapezium)
epsilon_a_simp = np.zeros(len(s)-1, dtype=float) # Array to store approximate relative error (Simpsons)

# Calculate integrals using Trapezium and Simpsons rules for different sampling intervals
for k in range(len(s)):
    Time = Time_Raw[0:last_point+1:s[k]]
    print(len(Time-1))
    Sampling[k] = Time[1] - Time[0]  # Sampling interval in seconds
    Velocity_squared = (Velocity_Raw[0:last_point+1:s[k]]) ** 2  # Squared velocity in (mm/s)^2

    I_TrapeziumRule[k] = integration.integrate_newton(Time, Velocity_squared, 'trap')
    I_SimpsonsRule[k] = integration.integrate_newton(Time, Velocity_squared, 'simp')

    print('Average squared velocity during the event using Trapezium rule is', (I_TrapeziumRule[k])/T, '(mm/s)^2 for',
           'Sampling interval=', Sampling[k], 's')
    print('Average squared velocity during the event using Simpsons rules is', (I_SimpsonsRule[k])/T, '(mm/s)^2 for',
           'Sampling interval=', Sampling[k], 's')
     

 # Calculate approximate relative errors
for kk in range(len(s)-1):
     
     epsilon_a_trap[kk] = abs((I_TrapeziumRule[kk] -  I_TrapeziumRule[kk+1])/(I_TrapeziumRule[kk]))
     epsilon_a_simp[kk] = abs((I_SimpsonsRule[kk] -  I_SimpsonsRule[kk+1])/(I_SimpsonsRule[kk]))

#Plot convergence plot and save it
plt.figure()
plt.loglog(Sampling[0:-1],epsilon_a_trap,'-ro',label='Trapezium Rule')
plt.loglog(Sampling[0:-1],epsilon_a_simp,'-bo',label='SimpsonsRule')

#custom_x_labels = [0.01, 0.02, 0.04, 0.08, 0.16, 0.32,0.64, 1.28, 2.56]
#plt.xticks(Sampling, custom_x_labels)
plt.xlabel('Sampling interval(s)')
plt.ylabel('Approximate Relative Error')
plt.title('Approximate relative error vs sampling intervals')
plt.legend()
plt.savefig('/Users/arvinkarpiah/Desktop/GOPH420/Lab_01/Repository/goph420_w2024_lab01_stAK/figures/Convergence plot_Q7.png')
plt.savefig('/Users/arvinkarpiah/Desktop/GOPH420/Lab_01/Repository/goph420_w2024_lab01_stAK/Convergence plot_Q7.png')
#plt.show()

       





