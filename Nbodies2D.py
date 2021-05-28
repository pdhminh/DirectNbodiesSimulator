import random as rd
import math
import numpy as np
from quadtree import *
from vpython import *
from decimal import Decimal
# from pprint import pprint

scene2 = canvas(title='N bodies simulator', x=0, y=0, width=800, height=800,
                center=vector(2500e6, 2500e6, 0))  # background=vector(1, 1, 1
G = 6.67408e-11
collided = False
planets = []
root = Node()


def flat_tree(Tree):
    for i in range(0, 4, 1):
        if isinstance(Tree['children'][i], Node):
            Tree['children'][i] = vars(Tree['children'][i])
            flat_tree(Tree['children'][i])
    return Tree


class Planets(object):
    def __init__(self, xo, yo, vox, voy, mass, radius):
        self.velocity = [vox, voy]
        self.mass = mass
        self.shape = simple_sphere(pos=vector(
            xo, yo, 0), radius=radius, color=vector(0, 1, 0))
        self.pos = [self.shape.pos.x, self.shape.pos.y]

    def get_x(self):
        return self.shape.pos.x

    def get_y(self):
        return self.shape.pos.y


def create_planet(x, y, vx, vy, m, d):
    V = m / d
    r = np.cbrt(V * (3 / 4) / np.pi)
    target = Planets(xo=x, yo=y, vox=vx, voy=vy, mass=m, radius=r)
    planets.append(target)
    root.insert(target, 2500e6, 2500e6, 5000e6, 5000e6)
    root.cal_approxiamtions(root)


def update_tree(payload, dimension):
    xo, yo, w, h = dimension
    # box(length=(y_max - y_min), height=(x_max - x_min), pos=vector((x_max + x_min)/2, (y_max + y_min)/2, 0))
    root.insert(payload, xo, yo, w, h)
    root.cal_approxiamtions(root)
    # bruh = root
    # pprint(flat_tree(vars(bruh)))


def mergingStars(p1, p2):
    if p1.mass > p2.mass:
        p1.velocity[0] = (p1.mass * p1.velocity[0] + p2.mass *
                          p2.velocity[0]) / (p1.mass + p2.mass)
        p1.velocity[1] = (p1.mass * p1.velocity[1] + p2.mass *
                          p2.velocity[1]) / (p1.mass + p2.mass)
        p1.mass = p1.mass + p2.mass
        V = p1.mass / d
        p1.shape.radius = np.cbrt(V * (3 / 4) / np.pi)
        ball = p2.shape
        ball.visible = False
        del ball
        planets.remove(p2)


def descend(root, target, data, kind, w):
    theta = 0.5
    w = float(w) / 2
    if root != target and root is not None:
        pos = [root.get_x(), root.get_y()]
        if kind == 0:
            distance = math.sqrt((pos[0] - data[0]) **
                                 2 + (pos[1] - target.pos[1])**2)
        else:
            distance = math.sqrt(
                (pos[0] - target.pos[0])**2 + (pos[1] - data[1])**2)
        if isinstance(root, Planets) or w / distance <= theta:
            F = G * root.mass * target.mass / (distance**2)
            return F * (pos[kind] - data[kind]) / distance
        else:
            res = 0
            for children in root.children:
                res = res + descend(root, target, data, kind, w)
            return res


def calca(target, kind, data):
    res = 0
    xo, yo, w, h = find_dimension()
    res = descend(root, target, data, kind, w)
    return res / target.mass


def ODE(planet, func, index, dt):
    # RK4th method
    data = [0, 0]
    data[index] = planet.pos[index]
    k1 = func(planet, index, data)
    data[index] = planet.pos[index] + planet.velocity[index] * dt / 2 + 0.5 * k1 * (dt / 2)**2
    k2 = func(planet, index, data)
    data[index] = planet.pos[index] + planet.velocity[index] * dt / 2 + 0.5 * k2 * (dt / 2)**2
    k3 = func(planet, index, data)
    data[index] = planet.pos[index] + planet.velocity[index] * dt + 0.5 * k3 * dt**2
    k4 = func(planet, index, data)
    a = (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    planet.pos[index] = planet.pos[index] + planet.velocity[index] * dt + 0.5 * a * dt ** 2
    planet.velocity[index] = planet.velocity[index] + a * dt
    return planet.pos[index]


def find_dimension():
    pos_max = [0, 0]
    pos_min = [1e10000, 1e1000]
    for planet in planets:
        planet.pos = [planet.shape.pos.x, planet.shape.pos.y]
        for i in range(0, 2, 1):
            if Decimal(planet.pos[i]) >= pos_max[i]:
                pos_max[i] = Decimal(planet.pos[i])
            if Decimal(planet.pos[i]) <= pos_min[i]:
                pos_min[i] = Decimal(planet.pos[i])
    x_max, y_max = pos_max
    x_min, y_min = pos_min
    return [(x_max + x_min) / Decimal(2), (y_max + y_min) / Decimal(2), x_max - x_min, y_max - y_min]


N = int(input("Bodies: "))
for i in range(1, N + 1, 1):
    x = rd.uniform(0, 5000e6)
    y = rd.uniform(0, 5000e6)
    vx = rd.uniform(0, 170.55 * 2)
    vy = rd.uniform(0, 170.55 * 2)
    m = rd.uniform(1.0224e+22, 2.1469e+24)
    d = 600
    create_planet(x, y, vx, vy, m, d)

print('=======================')
dt = float(input('time interval: '))
t = 0

# pprint(vars(root))
while True:
    rate(60)
    for planet in planets:
        planet.shape.pos = vector( ODE(planet, calca, 0, dt), ODE(planet, calca, 1, dt), 0)
        dem = find_dimension()
        root = root.remove(planet)
        update_tree(planet, dem)
    t = t + dt
