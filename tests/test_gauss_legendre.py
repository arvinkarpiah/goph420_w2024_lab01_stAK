#Test integrate_gauss() by designing unit tests in which data is generated from polynomials up to order n = 9 and 
#showing that the function can get the exact result with an appropriate option for npts. You may also try to write 
#unit tests for polynomial and non-polynomial functions comparing results with scipy.integrate.fixed_quad() which has the 
#parameter n for the number of integration points, but this is not mandatory. Implement any tests in tests/test_gauss_legendre.py 
#as described previously.

import unittest
import numpy as np
import sys
sys.path.append('/Users/arvinkarpiah/Desktop/GOPH420/Lab_01/Repository/goph420_w2024_lab01_stAK/src')
from goph420_lab01 import integration
from goph420_lab01.integration import integrate_gauss

# Data to test

# Dataset with even number of points to test for trapezium and Simpson's rule

def f(x, dx):
    result = (0.2 + 25 * x - 200 * x ** 2 + 675 * x ** 3 - 900 * x ** 4 + 400 * x ** 5) * dx
    return result

import unittest

class TestIntegration(unittest.TestCase):
     
    # Testing integrate_gauss for npts=2
    def test_integrategauss(self):
        exact_integral =   0.90 #1.822578
        self.assertAlmostEqual(integrate_gauss(f, [0, 0.8], 2), exact_integral, delta=1e-4)
       

if __name__ == '__main__':
    unittest.main()