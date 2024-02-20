#Computes the probability of a seismic event with magnitude greater than 4.0 assuming that magnitudes are normally distributed 
#with an annual mean of 1.5 and a standard deviation of 0.5. 
#Use Gauss-Legendre quadrature to do this. 
#You may also need to implement equation (17) as a function to do this.
#For a measured distance with mean L = 10.28 m and standard error ΔL=0.05 m, estimate the probability 
#that the true value is between 10.25m and 10.35 m. Use Gauss-Legendre quadrature to do this.
#Plots the convergence of the above estimates using increasing number of integration points in a log-log space.
#Saves the output figures in the directory figures/ in the root of your project folder. 
#You may want to use matplotlib.pyplot.savefig() to save your figures as image files.


#Import dependencies
import numpy as np
import matplotlib.pyplot as plt
import math
import sys
sys.path.append('/Users/arvinkarpiah/Desktop/GOPH420/Lab_01/Repository/goph420_w2024_lab01_stAK/src')
from goph420_lab01 import integration
from goph420_lab01.integration import integrate_gauss

#Question 8(i)

#Define mean and standard deviation
mu = 1.5
std = 0.5

#Define standard normal probability density function

def Phi(x,dx):  
    result = (1/math.sqrt(2*math.pi))*math.exp(-0.5*((x-mu)/(std))**2) * dx

#Define magnitude to be investigated, Z and convert to X

Z = 4
X = Z*std + mu

I_2pts = integration.integrate_gauss(Phi,[X,0],2)

#Probability_2pts = 0.5 - I_2pts


