
import unittest
import numpy as np
import sys
import math
sys.path.append('/Users/arvinkarpiah/Desktop/GOPH420/Lab_01/Repository/goph420_w2024_lab01_stAK/src')
from goph420_lab01 import integration
from goph420_lab01.integration import integrate_newton

# Data to test

# Dataset with even number of points to test for trapezium and Simpson's rule
x = np.linspace(-10.5, 10.5, num=601)
f = 2*x**2 - 8*x

# Dataset with odd number of points to test for trapezium and Simpson's rule
a = np.linspace(-10.5, 10.5, num=600)
b = -8*a -20

# Dataset with odd number of points to test for trapezium and Simpson's rule
c = np.linspace(0, (math.pi)/2, num=3)
d = (np.cos(c))*(np.sin(c))


class TestIntegration(unittest.TestCase):
     
     # Testing integrate_newton for both trapezium and Simpson's Rule for odd number of points
     def test_integratenewton_odd(self):
        exact_integral = float(((2/3)*(x[-1]**3)+4*x[-1]**2)-((2/3)*(x[0]**3)+4*x[0]**2))
        self.assertAlmostEqual(integrate_newton(x, f, 'trap'), exact_integral,delta=1e-4)
        self.assertAlmostEqual(integrate_newton(x, f, 'simp'), exact_integral,delta=1e-4)
     
    # Testing integrate_newton for both trapezium and Simpson's Rule for even number of points
     def test_integratenewton_even(self):
        exact_integral = float(((-4)*(a[-1]**2)-20*a[-1])-((-4)*(a[0]**2)-20*a[0]))
        self.assertAlmostEqual(integrate_newton(a, b, 'trap'), exact_integral,delta=1e-4)
        self.assertAlmostEqual(integrate_newton(a, b, 'simp'), exact_integral,delta=1e-4)

    # Testing integrate_newton for both trapezium and Simpson's Rule for odd number of points
     def test_integratenewton_GOPH420quiz1(self):
        exact_integral = 0.5
        #self.assertAlmostEqual(integrate_newton(c, d, 'trap'), exact_integral,delta=1e-4)
        self.assertAlmostEqual(integrate_newton(c, d, 'simp'), exact_integral,delta=1e-4)


if __name__ == '__main__':
    unittest.main()