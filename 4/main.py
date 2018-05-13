#! /usr/bin/env python3

import numpy as np
from math import sqrt, sin, pi, cos, e, exp
import scipy.integrate as integrate
#import matplotlib.pyplot as plt

f = lambda x: exp(x / 2)


'''
Chebyshev polynomials of the first kind. 
@n is order of polynom.
'''
def chebyshev(x, n):
    return 1 if n == 0 else x if n == 1 else 2.0 * x * chebyshev(x, n-1) - chebyshev(x, n - 2)

def legendre(x, n):
    return np.polynomial.legendre.legval(x, n)


'''
Solve Fredholm Integral Equations of the Second Kind 
by the collocation method.
'''
def fredgolm_I_solve(a, b, f, n):
    x = [(1 + cos(k * pi / (n - 1))) / 2.0 for k in range(n)]
    y = [f(s) for s in x]
    M = np.zeros((n,n))

    for i in range(n):
        for j in range(n):
            integr = integrate.quad(lambda t: (1. / (sqrt(1. + x[i] + t)) * chebyshev(2.0 * t - 1.0, j)), a, b) 
            M[i][j] = chebyshev(2.0 * x[i] - 1.0, j) - integr[0]

    A = np.linalg.solve(M, y)
    y = []
    x = np.linspace(a, b, n)
    for i in range(n):
        c = 0
        for j in range(n):
            c += A[j] * chebyshev(x[i], j)
        y.append(f(x[i]) - c)
        print(x[i], y[i])

    return x, y



def main():
    a, b = [0, 1]
    n    = int(input('Enter n: '))
    x, y = fredgolm_I_solve(a, b, f, n)
    '''
    t = 0.5
    for i in range(n):
    '''


if __name__ == '__main__':
	main()
