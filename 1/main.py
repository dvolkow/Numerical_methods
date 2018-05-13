#! /usr/bin/env python3
import math 


class Equation:
    # Edge conditions of III kind:
    alpha = [2.0, 0.0]
    beta = [-1.0, 1.0]
    A = 0.
    B = 0.
    a = 0.
    b = 1.

    def p(self, x):
        return (x + x ** 2) / (2.0 - x)

    def q(self, x):
        return math.log(x + (1 + x) ** 0.5, math.e)

    def f(self, x):
        return 1.0 + math.cos((math.pi / 2.0) * x)

    '''
    Return (a, b, c, d) values
    '''
    def toward_run(self, x, N):
        h = float(self.b - self.a) / N
        return 1 - (h / 2.) * self.p(x), 2 - (h ** 2) * self.q(x), 1 + (h / 2.) * self.p(x), (h ** 2) * self.f(x)




class ThomasCoeffs:
    def __init__(self, eq, N):
        h = float(eq.b - eq.a) / N
        # first:
        self.a = [0] 
        self.b = [1 + h / 2 * eq.p(eq.a) - h ** 2 / 2 * eq.q(eq.a) - eq.alpha[0] * h]
        self.c = [1 + h / 2 * eq.p(eq.a)]
        self.d = [eq.A * h + h ** 2 / eq.f(eq.a)]

        for i in range(1, N):
            a_, b_, c_, d_ = eq.toward_run(eq.a + i * h, N)
            self.a.append(a_)
            self.b.append(b_)
            self.c.append(c_)
            self.d.append(d_)

        # last:
        self.a.append(-1 + h / 2. * eq.p(eq.b))
        self.b.append(-1 + h / 2. * eq.p(eq.b) + h ** 2 / 2. * eq.q(eq.b) + eq.beta[0] * h)
        self.c.append(0)
        self.d.append(eq.B * h - h ** 2 / 2. * eq.f(eq.b))

    def get(self):
        return self.a, self.b, self.c, self.d



class Thomas:
    def __init__(self, N):
        self.e   = Equation()
        self.h   = float(self.e.b - self.e.a) / N
        self.X   = [self.e.a + self.h * j for j in range(N + 1)]
        self.Y   = [0.      for i in range(N + 1)]
        self.N   = N

    def solve(self):
        cs = ThomasCoeffs(self.e, self.N)
        a, b, c, d = cs.get()
        c_ = [c[0] / b[0]]
        d_ = [d[0] / b[0]]
        # forward sweep:
        for i in range(1, self.N):
            c_.append(c[i] / (b[i] - a[i] * c_[i - 1]))
            d_.append((d[i] - a[i] * d_[i - 1]) / (b[i] - a[i] * c_[i - 1]))

        # back substitution:
        self.Y[self.N] = (d[self.N] - a[self.N] * d_[self.N - 1]) / (b[self.N] - a[self.N] * c_[self.N - 1])
        for i in range(self.N - 1, 0 ,-1):
            self.Y[i] = d_[i] - c_[i] * self.Y[i + 1]

    def print(self):
        print("By Thomas:")
        print("X[i],    Y[i]")        
        for i in range(self.N + 1):
            print(self.X[i], "   ", self.Y[i])
        





def thomas_algorithm(N):
    Eq = Equation()
    h = float(Eq.b - Eq.a) / N

    X = [Eq.a + float(j) * (Eq.b - Eq.a) / N for j in range(N + 1)]
    Y = [0.0 for i in range(N + 1)]

#       a  b  c  d
    m, k, c, d, u, v = [], [], [], [], [], []
    kap, nu = [], []
	
    r = (4.0 - 2.0 * (2.0 - Eq.q(Eq.a + h) * h ** 2) / (2.0 + Eq.p(Eq.a + h) * h)) * Eq.alpha[1]
    z = 2.0 * Eq.alpha[0] * h - 3.0 * Eq.alpha[1] - Eq.alpha[1] * (Eq.p(Eq.a+h) * h - 2.0) / (2.0 + Eq.p(Eq.a+h) *h)

    kap.append((-r / z))
    nu.append(((2.0 * h * Eq.A + 2.0 * Eq.alpha[1] * Eq.f(Eq.a+h) * (h**2) / (2.0 + Eq.p(Eq.a+h) *h)) / z))

    r = (-4.0 + 2.0 * (2.0 - Eq.q(Eq.b-h) * (h**2)) / (2.0 - Eq.p(Eq.b-h) * h)) * Eq.beta[1]
    z = 2.0 * Eq.beta[0] * h + 3.0 * Eq.beta[1] - Eq.beta[1] * (Eq.p(Eq.b-h) * h + 2.0) / (2.0 - Eq.p(Eq.b-h) *h)
    kap.append((-r / z))
    nu.append(((2.0 * h * Eq.B - 2.0 * Eq.beta[1] * Eq.f(Eq.b-h) * (h**2) / (2.0 - Eq.p(Eq.b-h) * h)) /z))

    u.append(kap[0])
    v.append(nu[0])

    grid = ((j[0], Eq.a + float(j[1]) * (Eq.b - Eq.a) / N) for j in enumerate(range(0, N)))
	
    for (i, x) in grid:
        m.append(1.0 + (h / 2.0) * Eq.p(x))
        k.append(2.0 - (h ** 2) * Eq.q(x))
        c.append(1.0 - (h / 2.0) * Eq.p(x))
        d.append(((h ** 2) * Eq.f(x)))

    # toward sweep:
    for i in range(1, N):
        u.append((m[i]) / (k[i] - c[i] * u[i - 1]))
        v.append(((c[i] * v[i - 1] - d[i]) / (k[i] - c[i] * u[i - 1])))

    Y[N] = (- kap[1] * v[N - 1] - nu[1]) / (kap[1] * u[N - 1] - 1.0)
    # back substitution:
    for i in range(N - 1, 0 ,-1):
        Y[i] = u[i] * Y[i + 1] + v[i]

    Y[0] = u[0] * Y[1] + v[0]
    print("X[i],    Y[i]")        
    for i in range(N + 1):
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
    N = int(input("Enter N: "))
    thomas_algorithm(N)
    test = Thomas(N)
    test.solve()
    test.print()


if __name__ == '__main__':
	main()

