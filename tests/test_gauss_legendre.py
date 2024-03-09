#Tnhis script does the following:

#1. Creates a function that generates polynomials up to order of 9
#2. Tests the integrate_gauss function for polynomials up to order of 9
#3. Prints the result as a way to double check.

import unittest
import random
import math
import numpy as np
from scipy.integrate import fixed_quad
import sys
sys.path.append('/Users/arvinkarpiah/Desktop/GOPH420/Lab_01/Repository/goph420_w2024_lab01_stAK/src')
from goph420_lab01 import integration
from goph420_lab01.integration import integrate_gauss

# Create a function that generates polynomials of order n in the pattern x^n + x^n-1 + x^n-2 +....+x^(n-n)
def polynomial_generator(x, n):
        result = np.zeros_like(x)
        for i in range(n + 1):
            result += x ** i
        return result
    
    # Define integral function
def get_integral_function(j):
    return lambda x: polynomial_generator(x, j)

def F(t): 
    
    result = 1-math.exp((-12.5/68.1)*t)
    return result

# Generate bounds of integration
a=0
b=2
# Perform test
class TestIntegrationGauss(unittest.TestCase):
    def test_integrate_gauss_polynomial(self):
            for j in range(10):
                integral = get_integral_function(j)
                exact_integral, _ = fixed_quad(integral, a, b)
                computed_integral = integrate_gauss(integral,[a,b],5) 
                self.assertAlmostEqual(computed_integral, exact_integral,delta=1e-3)
                print('for n=',j,'exact_integral is',exact_integral,'and computed integral is',computed_integral)

    def test_integrate_gauss_example_from_book(self):
         E = (9.8*68.1/12.5)*(integrate_gauss(F,[0,10],5))
         self.assertAlmostEqual(E, 289.4351,delta=1e-3)
         

            
if __name__ == '__main__':
    unittest.main()


