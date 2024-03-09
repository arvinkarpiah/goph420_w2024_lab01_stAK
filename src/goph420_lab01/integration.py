#This script contains two functions that are can be called from other scripts:

# 1. integrate_newton: Performs integration of functions using the Trapezium and Simpsons' Rules.
# 2. integrate_gauss: Performs integration using the Gauss-Legendre method.

import math
from math import sqrt

####### integrate_newton

def integrate_newton(x, f, alg):
    if len(x) != len(f):
        raise ValueError("Input lists must have the same length")
    
    alg = alg.strip().lower()  # Convert to lowercase and remove leading/trailing whitespace
    
    valid_algs = ['trap', 'simp']
    if alg not in valid_algs:
        raise ValueError(f"Invalid integration rule '{alg}'. Use 'trap' or 'simp'.")

    a = x[0]  # First time value
    b = x[-1]  # Last time value
    f0 = f[0]  # First velocity value
    fN = f[-1]  # Last velocity value
    N = len(x) - 1  # Number of segments

    if alg == 'trap':
        total_sum = sum(f[1:-1])  # sum of all elements except the first and last

        I = ((b - a) / (2 * N)) * (f0 + 2 * total_sum + fN)

        
    elif alg == 'simp':
        total_sum_1 = sum(f[1:-1:2])  # sum of odd elements except the first and last
        total_sum_2 = sum(f[2:-1:2])  # sum of even elements except the first and last

        total_sum_3 = sum(f[1:-4:2])  # sum of odd elements except the first and last
        total_sum_4 = sum(f[2:-4:2])  # sum of even elements except the first and last

        if len(x) % 2 == 1:  # If number of data points are odd
            I = ((b - a) / (3 * N)) * (f0 + 4 * total_sum_1 + 2 * total_sum_2 + fN)

        else:  # If number of data points are even
            I_2g = ((x[-4] - x[0]) / (3 * (N - 3))) * (f0 + 4 * total_sum_3 + 2 * total_sum_4 + f[-4]) #Simpson's 1/3 Rule
            I_3 = (1 / 8) * (x[-1] - x[-4]) * (f[-4] + 3 * f[-3] + 3 * f[-2] + f[-1]) #Simpson's 3/8 Rule
            I = I_2g + I_3

    return float(I)          

####### integrate_gauss

def integrate_gauss(f, lims, npts):
    if len(lims) != 2:
        raise ValueError("lims must contain exactly two values")
    
    try:
        a, b = map(float, lims)
    except (ValueError, TypeError):
        raise ValueError("lims must be convertible to float")

    valid_npts = [1, 2, 3, 4, 5]
    if npts not in valid_npts:
        raise ValueError(f"Invalid value for npts. Must be one of the following values: {valid_npts}")

    if callable(f):
        # Change limits
        A = (b + a) / 2
        B = (b - a) / 2

        if npts == 1:
            c_0 = 2.0
            x_0 = 0.0

            f_0 = f(A + B * x_0)
                
            I = B * c_0 * f_0
       

        elif npts == 2:
            c_0 = 1
            c_1 = 1
            x_0 = -1/sqrt(3)
            x_1 = 1/sqrt(3)
    

            f_0=f(A + B * x_0)
            f_1=f(A + B * x_1)

            I = B*(c_0 * f_0 + c_1 * f_1)
                  

        elif npts == 3:

            c_0 = 5/9
            c_1 = 8/9
            c_2 = 5/9
            x_0 = -sqrt(3/5)
            x_1 = 0
            x_2 = sqrt(3/5)

            f_0=f(A + B * x_0)
            f_1=f(A + B * x_1)
            f_2=f(A + B * x_2)
                
            I = B*(c_0 * f_0 + c_1 * f_1 + c_2 * f_2)

        if npts == 4:

            c_0 = (18-sqrt(30))/36
            c_1 = (18+sqrt(30))/36
            c_2 = (18+sqrt(30))/36
            c_3 = (18-sqrt(30))/36
            x_0 = -sqrt(3/7+(2/7)*sqrt(6/5))
            x_1 = -sqrt(3/7-(2/7)*sqrt(6/5))
            x_2 = sqrt(3/7-(2/7)*sqrt(6/5))
            x_3 = sqrt(3/7+(2/7)*sqrt(6/5))

            f_0=f(A + B * x_0)
            f_1=f(A + B * x_1)
            f_2=f(A + B * x_2)
            f_3=f(A + B * x_3)
    
            I = B*(c_0 * f_0 + c_1 * f_1 + c_2 * f_2 + c_3 * f_3)


        if npts == 5:


            c_0 = (322+13*sqrt(70))/900
            c_1 = (322-13*sqrt(70))/900
            c_2 = 128/225
            c_3 = (322-13*sqrt(70))/900
            c_4 = (322+13*sqrt(70))/900
            x_0 = -(1/3)*sqrt(5-2*sqrt(10/7))
            x_1 = -(1/3)*sqrt(5+2*sqrt(10/7))
            x_2 = 0.0
            x_3 = (1/3)*sqrt(5+2*sqrt(10/7))
            x_4 = (1/3)*sqrt(5-2*sqrt(10/7))

            f_0=f(A + B * x_0)
            f_1=f(A + B * x_1)
            f_2=f(A + B * x_2)
            f_3=f(A + B * x_3)
            f_4=f(A + B * x_4)
    
            I = B*(c_0 * f_0 + c_1 * f_1 + c_2 * f_2 + c_3 * f_3 + c_4 * f_4)

        return float(I)
    else:
        raise TypeError("f is not a callable")

