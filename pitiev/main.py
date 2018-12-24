#! /usr/bin/env python3

RES_F = "res.txt"
x_0 = 0.01 + 0.005 * 6
y_0 = 0.01 + 0.005 * 6
z_0 = 0 
b = 8.0/3
sigma = 10.0
r = 21 - 6
delta_t = 1e-4

def get_dot_x(x, y, z):
    return sigma * (-x + y)

def get_dot_y(x, y, z):
    return r * x - y - x * z

def get_dot_z(x, y, z):
    return -b * z + x * y


def get_new_point(x, y, z, delta_t):
    x += get_dot_x(x, y, z) * delta_t
    y += get_dot_y(x, y, z) * delta_t
    z += get_dot_z(x, y, z) * delta_t
    return (x, y, z)

def do_loop(delta_t):
    t = 0
    x, y, z = (x_0, y_0, z_0)
    with open(RES_F, 'w') as f:
        while t < 10: 
            x, y, z = get_new_point(x, y, z, delta_t)
            f.write('{0} {1} {2}\n'.format(x, y, z))
            t += delta_t


if __name__ == '__main__':
    do_loop(delta_t)
