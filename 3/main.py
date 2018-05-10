#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import random as rand
from numpy.linalg import norm
from numpy.linalg import inv
import math
import functools
import fractions

gs = int(input('Enter s: '))
ga = int(input('Enter a: '))
gb = int(input('Enter b: '))

PRECISION = 1e-10
TYPE_MAX  = 0
TYPE_MIN  = 1


'''
Set value for my matrix
'''
def matrix_init(s, a, b):
    m = [[0]*s for _ in range(s)]
    for i in range(s):
        for k in range(s):
            if i == k:
                m[i][i] = i * (s + 1 - i) * (4 / (s + 1) ** 2) * a
            elif abs(i - k) > 1: 
                m[i][k] = 0.1 / ((i - k) ** 2 + 1)
            else:
                m[i][k] = b 
                m[k][i] = b 
    return np.mat(m)

m = matrix_init(gs, ga, gb)
print(m)
w, v = np.linalg.eig(m)
print("Maximum abs value:" + str(max(np.abs((w)))))
print("Minimum abs value:" + str(min(np.abs((w)))))
print("Minimum value:" + str(min((w))))


def degress_method(A, s, TYPE):
    E  = np.array([rand.random()] * s)
    E.resize(s, 1)
    x  = np.dot(A, E)
    x1 = float(x[0][0]) if TYPE == TYPE_MAX else 1/float(x[0][0])
    d  = abs(norm(x, ord = 2) - norm(E, ord = 2))
    n  = 1

    while d > PRECISION: 
        E  = np.array(x)
        x  = np.dot(A, E)
        x1 = float(x[0][0]) if TYPE == TYPE_MAX else 1/float(x[0][0])
        n += 1
        res = np.array(x)
        x /= x1
        d = abs(norm(E, ord=2) - norm(x, ord=2))
    return res, n


lambdamax, i = degress_method(m, gs, TYPE_MAX)
maxl = lambdamax[0, 0]
print(lambdamax)
print("Maximum by degress method: ", maxl, " by ", i,
        " iterations")

print(m)
print(np.identity(len(m)))
m -= np.identity(len(m)) * maxl
print(m)
#ainv = inv(m) 
lambdamin, i = degress_method(m, gs, TYPE_MAX)
minl = lambdamin[0, 0] + maxl
print(lambdamin)
print("Minimum by degress method: ", minl, " by ", i,
        " iterations")

