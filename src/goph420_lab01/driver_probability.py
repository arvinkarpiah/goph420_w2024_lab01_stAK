#This script does the following:

#1. Computes the probability of a seismic event with magnitude greater than 4.0 
#   using Gauss-Legendre quadrature method 
#2. Estimates the probability that the true value is between 10.25m and 10.35 m using 
#   Gauss-Legendre quadrature method
#3. Plots the convergence of the above estimates using increasing number of integration points
#4. Saves the output figures XXXXX in the directory figures 

#Import dependencies and functions
import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import derivative
import random
import math
# import sys
# sys.path.append('/Users/arvinkarpiah/Desktop/GOPH420/Lab_01/Repository/goph420_w2024_lab01_stAK/src')
# from goph420_lab01 import integration
# from goph420_lab01.integration import integrate_gauss
import integration
from integration import integrate_gauss

####################### Question 8(i)

#Define mean and standard deviation
mu = 1.5
std = 0.5

#Define normal probability density function
def Phi(z): 
    result = (1/(math.sqrt(2*math.pi))*math.exp(-0.5*(z**2)))
    return result

#Transform inputs to standard normal variate z

X = 4
Z = (X-mu)/(std) 

#Compute integration with different number of points

xnpts = np.array([1,2,3,4,5])

Probability_GL1 = np.zeros(len(xnpts), dtype=float)

for k in range(len(xnpts)):
 
    Probability_GL1[k] = 0.5 - integration.integrate_gauss(Phi,[0,Z],xnpts[k])

    print('Probability of a seismic event with magnitude greater than 4.0 computed using Gauss-Legendre Quadrature with', xnpts[k], 'points is', Probability_GL1[k])

# Generate a random number between the limts to evaluate f**2N at
chi = random.uniform(0,Z)
TruncErr = np.zeros(len(xnpts), dtype=float)

for i in range(len(xnpts)):
    T3 = derivative(Phi, chi, n=2*xnpts[i], order=2*xnpts[i]+1)
    T1 = (Z-0)**(2*xnpts[i]+1)
    T2 = (math.factorial(xnpts[i]))**4
    T4 = (2*xnpts[i]+1)
    T5 = (math.factorial(2*xnpts[i]))**3
    
    TruncErr[i] = abs((T1*T2*T3)/(T4*T5))
                                    
#Plot convergence plot and save it
plt.figure()
custom_x_labels = [1, 2, 3, 4, 5]
plt.loglog(xnpts,TruncErr,'-ro')
plt.xlabel('Number of points used in integration')
plt.xticks(xnpts, custom_x_labels)
plt.ylabel('Truncation Error')
plt.title('Truncation error vs number of points')
plt.savefig('/Users/arvinkarpiah/Desktop/GOPH420/Lab_01/Repository/goph420_w2024_lab01_stAK/figures/Convergence_GaussLegendre_Q8i.png')
plt.savefig('/Users/arvinkarpiah/Desktop/GOPH420/Lab_01/Repository/goph420_w2024_lab01_stAK/Convergence_GaussLegendre_Q8i.png')

####################### Question 8(ii)

#Define mean and standard deviation
mu=10.28
std=0.05

X1 = 10.35
Z1 = (X1-mu)/std

X2 = 10.25
Z2 = (X2-mu)/std

#Compute integration with different number of points

Probability_GL2 = np.zeros(len(xnpts), dtype=float)

for k in range(len(xnpts)):

    Probability_GL2[k] = integration.integrate_gauss(Phi,[0,Z1],xnpts[k]) +  integration.integrate_gauss(Phi,[0,-1*Z2],xnpts[k])

    print('Probability of true value being between 10.25 and 10.35 computed using Gauss-Legendre Quadrature with', xnpts[k], 'points is',Probability_GL2[k])


#Compute trunctation error

# Generate a random number between the limts to evaluate f**2N at               
chi = random.uniform(Z1,Z2)

for i in range(len(xnpts)):
    T3 = derivative(Phi, chi, n=2*xnpts[i], order=2*xnpts[i]+1)
    T1 = (Z2-Z1)**(2*xnpts[i]+1)
    T2 = (math.factorial(xnpts[i]))**4
    T4 = (2*xnpts[i]+1)
    T5 = (math.factorial(2*xnpts[i]))**3
    
    TruncErr[i] = abs((T1*T2*T3)/(T4*T5))

                                   
#Plot convergence plot and save it
plt.figure()
plt.loglog(xnpts,TruncErr,'-bo')
plt.xlabel('Number of points used in integration')
plt.xticks(xnpts, custom_x_labels)
plt.ylabel('Truncation Error')
plt.title('Truncation error vs number of points')
plt.legend()
plt.savefig('/Users/arvinkarpiah/Desktop/GOPH420/Lab_01/Repository/goph420_w2024_lab01_stAK/figures/Convergence_GaussLegendre_Q8ii.png')
plt.savefig('/Users/arvinkarpiah/Desktop/GOPH420/Lab_01/Repository/goph420_w2024_lab01_stAK/Convergence_GaussLegendre_Q8ii.png')
#plt.show()


