#! /usr/bin/env python3
import numpy as np
from math import sqrt, sin, pi, cos
import scipy.integrate as integrate
#import matplotlib.pyplot as plt

f = lambda x: x+sqrt(1.0+x)

'''
Chebyshev polynomials of the first kind. 
@n is order of polynom.
'''
def chebyshev(x, n):
    return 1 if n == 0 else x if n == 1 else 2.0 * x * chebyshev(x, n-1) - chebyshev(x, n - 2)


'''
Solve Fredholm Integral Equations of the Second Kind 
by the collocation method.
'''
def fredgolm(a, b, f, n):
    x = [(1+cos(k*pi/(n-1)))/2.0 for k in range(n)]
    y = [f(s) for s in x]
    M = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            integr = integrate.quad(lambda t: ((1.0/2.0*x[i]*t 
                + 2*sqrt(1.0 + 1.0/12.0*(x[i]+t)))*chebyshev(2.0*t-1.0,j)), a, b) 
            M[i][j] = chebyshev(2.0*x[i]-1.0,j) - integr[0]
    A = np.linalg.solve(M,y)
    y = []
    x = np.linspace(a,b,n)
    for i in range(n):
        c = 0
        for j in range(n):
            c += A[j]*chebyshev(x[i],j)
            y.append(f(x[i])-c)
        print(x[i], y[i])
    return x, y


def main():
    a, b = [0, 1]
    n = int(input('Enter n: '))
    x, y = fredgolm(a,b,f,n)
    print(x, y)


if __name__ == '__main__':
	main()
