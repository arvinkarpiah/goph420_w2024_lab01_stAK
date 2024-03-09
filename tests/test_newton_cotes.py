#This script does the following:

# 1.Tests the integrate_newton function with different combinations of
# a)Odd and even number of data points
# b)Trapezium and Simpsons' Rules.
# 2.Prints the result as a way to double check.

import unittest
import numpy as np
import sys
import math
sys.path.append('/Users/arvinkarpiah/Desktop/GOPH420/Lab_01/Repository/goph420_w2024_lab01_stAK/src')
from goph420_lab01 import integration
from goph420_lab01.integration import integrate_newton

# Dataset with even number of points to test for trapezium and Simpson's rule
x = np.linspace(-10.5, 10.5, num=600)
f = 2*x**2 - 8*x

# Dataset with odd number of points to test for trapezium and Simpson's rule
a = np.linspace(-10.5, 10.5, num=601)
b = -8*a -20


class TestIntegrationNewton(unittest.TestCase):
     
    # Testing integrate_newton for both trapezium and Simpson's Rule for even number of points
    def test_integratenewton_even(self):
        exact_integral_even = (2/3)*(x[-1]**3) - 4*x[-1]**2 - ((2/3)*(x[0]**3) - 4*x[0]**2)
        computed_integral_trap_even = integrate_newton(x, f, 'trap')
        computed_integral_simp_even = integrate_newton(x, f, 'simp')
        self.assertAlmostEqual(computed_integral_trap_even, exact_integral_even, delta=1e-2)
        self.assertAlmostEqual(computed_integral_trap_even, exact_integral_even, delta=1e-2)
        print('for even num.of points with Trapezium Rule','exact_integral is',exact_integral_even,'and computed integral is',computed_integral_trap_even)
        print('for even num.of points with Simpsons Rules','exact_integral is',exact_integral_even,'and computed integral is',computed_integral_simp_even)

    #Testing integrate_newton for both trapezium and Simpson's Rule for odd number of points
    def test_integratenewton_odd(self):
       exact_integral_odd = (-4)*(a[-1]**2) - 20*a[-1] - ((-4)*(a[0]**2) - 20*a[0])
       computed_integral_trap_odd = integrate_newton(a, b, 'trap')
       computed_integral_simp_odd = integrate_newton(a, b, 'simp')
       self.assertAlmostEqual(computed_integral_trap_odd, exact_integral_odd, delta=1e-6)
       self.assertAlmostEqual(computed_integral_simp_odd, exact_integral_odd, delta=1e-6)
       print('for even num.of points with Trapezium Rule','exact_integral is',exact_integral_odd,'and computed integral is',computed_integral_trap_odd)
       print('for even num.of points with Simpsons Rules','exact_integral is',exact_integral_odd,'and computed integral is',computed_integral_simp_odd)


if __name__ == '__main__':
    unittest.main()
