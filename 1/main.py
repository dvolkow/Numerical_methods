#! /usr/bin/env python3
import math 
#import matplotlib.pyplot as plt

p = lambda x: (x + x * x) / (2.0 - x)
q = lambda x: math.log(x + (1 + x) ** 0.5, math.e)
f = lambda x: 1.0 + math.cos((math.pi / 2.0) * x)

# -- Edge conditions of III kind:
alpha = [2.0, 0.0]
beta = [-1.0, 1.0]
A, B = 0.0, 0.0

a, b = 0.0, 1.0

N = int(input("Enter N: "))


def thomas_algorithm(N, a = a, b = b, A = A, B = B ):
	X = [a + float(j) * (b - a) / N for j in range(N + 1)]
	Y = [0.0 for i in range(N + 1)]

	m, k, c, d, u, v = [], [], [], [], [], []
	h = (b - a) / N
	kap, nu = [], []
	
	r = (4.0 - 2.0*(2.0 - q(a+h) * h ** 2) / (2.0 + p(a+h) * h)) * alpha[1]
	z = 2.0 * alpha[0] * h - 3.0 * alpha[1] - alpha[1] * (p(a+h) * h - 2.0) / (2.0 + p(a+h) *h)

	kap.append((-r / z))
	nu.append(((2.0 * h * A + 2.0 * alpha[1] * f(a+h) * (h**2) / (2.0 + p(a+h) *h)) / z))

	r = (-4.0 + 2.0 * (2.0 - q(b-h) * (h**2)) / (2.0 - p(b-h) * h)) * beta[1]
	z = 2.0 * beta[0] * h + 3.0 * beta[1] - beta[1] * (p(b-h) * h + 2.0) / (2.0 - p(b-h) *h)
	kap.append((-r / z))
	nu.append(((2.0 * h * B - 2.0 * beta[1] * f(b-h) * (h**2) / (2.0 - p(b-h) * h)) /z))

	u.append(kap[0])
	v.append(nu[0])

	grid = ((j[0], a + float(j[1]) * (b-a) / N) for j in enumerate(range(0, N)))
	
	for (i, x) in grid:
		m.append(1.0 + (h/2.0)*p(x))
		k.append(2.0 - (h**2)*q(x))
		c.append(1.0 - (h/2.0)*p(x))
		d.append(((h**2)*f(x)))

	for i in range(1, N):
		u.append((m[i])/(k[i] - c[i]*u[i-1]))
		v.append(((c[i]*v[i-1]-d[i])/(k[i]-c[i]*u[i-1])))
#		print(u[i], v[i])

	Y[N] = (-kap[1]*v[N-1]-nu[1])/(kap[1]*u[N-1]-1.0)
	
	for i in range(N-1, 0 ,-1):
		Y[i] = u[i]*Y[i+1] + v[i]

	Y[0] = u[0]*Y[1] + v[0]
	print("X[i],    Y[i]")        
	for i in range(N+1):
		print(X[i], "   ", Y[i])
	return X, Y


'''
plt.title('The graphs of the solution of equation:')
plt.ylabel('y')
plt.xlabel('x')
plt.grid(True)
plt.xlim(a,b)

def plot_result(N, a = a, b = b, A = A, B = B):
	x, y = thomas_algorithm(N, a, b, A, B)
	plt.plot(x, y)
	
plot_result(N)
plt.savefig('result.png')
plt.show()
'''



def main():
    thomas_algorithm(N)


if __name__ == '__main__':
	main()

