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
def sign(a_pp, a_qq, a_pq):
    return -1 if a_pq * (a_pp - a_qq) < 0 else 1

def sign_generic(expr):
    return -1 if expr < 0 else 1

def cs_calc(a_pp, a_qq, a_pq):
    d = math.sqrt((a_pp - a_qq) ** 2 + 4 * a_pq ** 2)
    c = math.sqrt(0.5 * (1 + math.fabs(a_pp - a_qq) / d))
    s = sign(a_pp, a_qq, a_pq) * math.sqrt(0.5 * (1 - math.fabs(a_pp - a_qq) / d)) 
    return c, s, d


def jacoby_rotate(m_a, p, q):
    '''
    if math.fabs(m_a[p][q]) < PRECISION:
        return m_a
    '''

    c, s, d = cs_calc(m_a[p, p], m_a[q, q], m_a[p, q])
    m_c     = np.array(m_a)
    m_c[p, p] = (m_a[p, p] + m_a[q, q]) / 2 + sign_generic(m_a[p, p] - m_a[q, q]) * d / 2
    m_c[q, q] = (m_a[p, p] + m_a[q, q]) / 2 - sign_generic(m_a[p, p] - m_a[q, q]) * d / 2
    m_c[p, q] = 0
    m_c[q, p] = 0
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
    return hn
    

'''
'''
def jacoby(matrix):
    while(matrix_half_norm(matrix) > PRECISION):
        for i in range(len(matrix) - 1):
            for j in range(i + 1, len(matrix)):
                matrix = jacoby_rotate(matrix, i, j)
    return matrix

'''
Get <<true>> values
'''
m = matrix_init(gs, ga, gb)
print(m)
w, v = np.linalg.eig(m)
print("Maximum abs value:" + str(max(np.abs((w)))))
print("Minimum abs value:" + str(min(np.abs((w)))))
print("Minimum value:" + str(min((w))))


print("Matrix after jacoby procedure:")
print(jacoby(matrix_init(gs, ga, gb)))


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

