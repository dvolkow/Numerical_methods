#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import random as rand
from numpy.linalg import norm
from numpy.linalg import inv
import math
import functools
import fractions

PRECISION = 1e-11
print('Precision is set as: ', PRECISION)

gs = int(input('Enter s: '))
ga = float(input('Enter a: '))
gb = float(input('Enter b: '))

TYPE_MAX  = 0
TYPE_MIN  = 1


'''
Set value for my matrix
'''
def matrix_init(s, a, b):
    m = [[0]*s for _ in range(s)]
    for i in range(s):
        for k in range(i, s):
            if i == k:
                m[i][i] = i * (s + 1 - i) * (4 / (s + 1) ** 2) * a
            elif abs(i - k) > 1: 
                m[i][k] = 0.1 / ((i - k) ** 2 + 1)
                m[k][i] = m[i][k]
            else:
                m[i][k] = b 
                m[k][i] = b 
    return np.mat(m)




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


'''
Here a_pp != a_qq 
'''
def sign_generic(expr):
    return -1 if expr < 0 else 1

def sign(a_pp, a_qq, a_pq):
    return sign_generic(a_pq * (a_pp - a_qq)) 

def cs_calc(a_pp, a_qq, a_pq):
    d = math.sqrt((a_pp - a_qq) ** 2 + 4. * a_pq ** 2)
    c = math.sqrt(0.5 * (1 + math.fabs(a_pp - a_qq) / d))
    s = sign(a_pp, a_qq, a_pq) * math.sqrt(0.5 * (1 - math.fabs(a_pp - a_qq) / d)) 
    return c, s, d


'''
Rotate m_a with angle for 
'''
def jacoby_rotate(m_a, p, q):
    '''
    Probably you can use next condition
    for stop calculation:

    if math.fabs(m_a[p][q]) < PRECISION:
        return m_a
    '''

    c, s, d = cs_calc(m_a[p, p], m_a[q, q], m_a[p, q])
    m_c     = np.array(m_a)
    m_c[p, p] = (m_a[p, p] + m_a[q, q]) / 2. + sign_generic(m_a[p, p] - m_a[q, q]) * d / 2.
    m_c[q, q] = (m_a[p, p] + m_a[q, q]) / 2. - sign_generic(m_a[p, p] - m_a[q, q]) * d / 2.
    m_c[p, q] = 0.
    m_c[q, p] = 0.
    for i in range(len(m_a)):
        if i != p and i != q:
            m_c[i, p] = c * m_a[i, p] + s * m_a[i, q]
            m_c[p, i] = m_c[i, p]
            m_c[i, q] = -s * m_a[i, p] + c * m_a[i, q]
            m_c[q, i] = m_c[i, q]
    return m_c

def matrix_half_norm(matrix):
    hn = 0
    for i in range(len(matrix) - 1):
        for j in range(i + 1, len(matrix)):
            hn += matrix[i, j] ** 2
    return math.sqrt(2 * hn)
   

def trace(matrix):
    tr = 0
    for i in range(len(matrix)):
        tr += matrix[i, i]
    return tr

'''
'''
def jacoby(matrix):
    while(matrix_half_norm(matrix) > PRECISION):
    #    print(matrix_half_norm(matrix))
        for i in range(len(matrix) - 1):
            for j in range(i + 1, len(matrix)):
                matrix = jacoby_rotate(matrix, i, j)
    return matrix

'''
Get <<true>> values
'''
m = matrix_init(gs, ga, gb)
w, v = np.linalg.eig(m)
print("'True' Maximum abs value: " + str(max(np.abs((w)))))
print("'True' Minimum value: " + str(min((w))))


print("\nMatrix before jacoby procedure:")
print(matrix_init(gs, ga, gb))
jm = jacoby(matrix_init(gs, ga, gb))
print("\nMatrix after jacoby procedure:")
jm = jacoby(matrix_init(gs, ga, gb))
print("Trace:               ", trace(matrix_init(gs, ga, gb)))
print("Sum of self numbers: \n", trace(jm))
print(jm)
 

print("Self numbers by Jacoby:")
for i in range(len(jm)):
    print(jm[i, i])


print("\nMatrix before degress procedure:")
print(matrix_init(gs, ga, gb))
lambdamax, i = degress_method(matrix_init(gs, ga, gb), gs, TYPE_MAX)
maxl = lambdamax[0, 0]
print("Maximum by degress method: ", maxl, " by ", i,
        " iterations")

m -= np.identity(len(m)) * maxl
lambdamin, i = degress_method(m, gs, TYPE_MAX)
minl = lambdamin[0, 0] + maxl
print("Minimum by degress method: ", minl, " by ", i,
        " iterations")

