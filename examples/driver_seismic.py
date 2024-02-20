#Loads the velocity data from the file s_wave_data.txt. You may want to use numpy.loadtxt() to do this.
#Plots the raw data and indicates the value T used in the integrations.
#Estimates the integral in equation (14) using different sampling intervals and both trapezoid rule and Simpson's rules.
#Calculates and plots the convergence of the integral estimates using trapezoid rule and Simpson's rules on a log-log plot. Y
#You may want to use the function matplotlib.pyplot.loglog() which has similar functionality to matplotlib.pyplot.plot(). 
#Make sure to include labels and a legend indicating which convergence curve is for trapezoid rule and which is for Simpson's rules.
#Saves the output figures in the directory figures/ in the root of your project folder. 
#You may want to use matplotlib.pyplot.savefig() to save your figures as image files.


#Import dependencies
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/Users/arvinkarpiah/Desktop/GOPH420/Lab_01/Repository/goph420_w2024_lab01_stAK/src')
from goph420_lab01 import integration
from goph420_lab01.integration import integrate_newton

#Load seismic data
data_path = '/Users/arvinkarpiah/Desktop/GOPH420/Lab_01/Repository/goph420_w2024_lab01_stAK/examples/s_wave_data.txt'
SeismicData = np.loadtxt(data_path)
Time_Raw = SeismicData[:,0]
Velocity_Raw = SeismicData[:,1]

#Analysis for an appropriate period of event, T to use
max_abs_vel = max(abs(Velocity_Raw))

bin_list = []  # Initialize an empty list to store indices

for i in range(len(Velocity_Raw)):
    if Velocity_Raw[i] > (0.5/100) * max_abs_vel:
        bin_list.append(i)  

last_point = max(bin_list) #Last point in the series that satisfies above condition

#Plot raw data
plt.figure()
plt.plot(Time_Raw,Velocity_Raw,'b-')
plt.axvline(x=Time_Raw[last_point], color='r', linestyle='-')
plt.title('S-wave velocity vs arrival time ')
plt.xlabel('Time(s)')
plt.ylabel('Velocity(mm/s)')
plt.text(Time_Raw[last_point], max(Velocity_Raw), 'T=6.46s', verticalalignment='bottom', horizontalalignment='right')
plt.show()
plt.savefig('/Users/arvinkarpiah/Desktop/GOPH420/Lab_01/Repository/goph420_w2024_lab01_stAK/figures/Raw data.png')

# Compute integration for different sampling intervals and rules

s = np.array([1, 2, 4, 8, 16, 32])  # Sampling intervals

T = Time_Raw[last_point]  # Period in seconds
print('T=', T)

I_TrapeziumRule = np.zeros(len(s), dtype=float)  # Array to store results for Trapezium rule
I_SimpsonsRule = np.zeros(len(s), dtype=float)  # Array to store results for Simpsons rule
Sampling = np.zeros(len(s), dtype=float)  # Array to store sampling intervals
epsilon_a_trap = np.zeros(len(s)-1, dtype=float) # Array to store approximate relative error (trapezium)
epsilon_a_simp = np.zeros(len(s)-1, dtype=float) # Array to store approximate relative error (Simpsons)

for k in range(len(s)):
    Time = Time_Raw[0:last_point+1:s[k]]
    Sampling[k] = Time[1] - Time[0]  # Sampling interval in seconds
    Velocity_squared = (Velocity_Raw[0:last_point+1:s[k]]) ** 2  # Squared velocity in (mm/s)^2

    # Calculate integrals using Trapezium and Simpsons rules
    I_TrapeziumRule[k] = integration.integrate_newton(Time, Velocity_squared, 'trap')
    I_SimpsonsRule[k] = integration.integrate_newton(Time, Velocity_squared, 'simp')


    print('Average squared velocity during the event using Trapezium rule is', (I_TrapeziumRule[k])/T, '(mm/s)^2 for',
           'Sampling interval=', Sampling[k], 's')
    print('Average squared velocity during the event using Simpsons rules is', (I_SimpsonsRule[k])/T, '(mm/s)^2 for',
           'Sampling interval=', Sampling[k], 's')
 

for kk in range(len(s)-1):
     
     # Calculate approximate relative errors
     epsilon_a_trap[kk] = abs((I_TrapeziumRule[kk+1] -  I_TrapeziumRule[kk])/(I_TrapeziumRule[kk+1]))
     epsilon_a_simp[kk] = abs((I_SimpsonsRule[kk+1] -  I_SimpsonsRule[kk])/(I_SimpsonsRule[kk+1]))

plt.figure()
plt.loglog(Sampling[1:],epsilon_a_trap,'-ro',label='Trapezium Rule')
plt.loglog(Sampling[1:],epsilon_a_simp,'-bo',label='SimpsonsRule')
plt.xlabel('Sampling interval(s)')
plt.ylabel('Approximate Relative Error')
plt.title('Approximate relative error vs sampling intervals')
plt.legend()
plt.savefig('/Users/arvinkarpiah/Desktop/GOPH420/Lab_01/Repository/goph420_w2024_lab01_stAK/figures/Convergence.png')
plt.show()



    
    





