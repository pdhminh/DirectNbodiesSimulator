import random as rd
import math
import numpy as np
from vpython import *

scene2 = canvas(title='N bodies simulator', x=0, y=0, width=800, height=800, center=vector(0, 0, 0)) # background=vector(1, 1, 1
G = 6.67408e-11
collided = False
planets = []


class Planets(object):
    def __init__(self, xo, yo, zo, vox, voy, voz, mass, radius):
        self.velocity = [vox, voy, voz]
        self.mass = mass
        self.shape = simple_sphere(pos=vector(xo, yo ,zo), radius=radius, color=vector(0, 1, 0))
        self.pos = [self.shape.pos.x, self.shape.pos.y, self.shape.pos.z]


def create_planet(x, y, z, vx, vy, vz, m, d):
    V = m/d
    r = np.cbrt(V*(3/4)/np.pi)
    planets.append(Planets(xo=x, yo=y, zo=z, vox=vx, voy=vy, voz=vz, mass=m, radius=r))


def mergingStars(p1, p2):
    if p1.mass > p2.mass:
        p1.velocity[0] = (p1.mass*p1.velocity[0] + p2.mass*p2.velocity[0])/(p1.mass + p2.mass)
        p1.velocity[1] = (p1.mass*p1.velocity[1] + p2.mass*p2.velocity[1])/(p1.mass + p2.mass)
        p1.velocity[2] = (p1.mass*p1.velocity[2] + p2.mass*p2.velocity[2])/(p1.mass + p2.mass)
        p1.mass = p1.mass + p2.mass
        V = p1.mass / d
        p1.shape.radius = np.cbrt(V*(3/4)/np.pi)
        ball = p2.shape
        ball.visible = False
        del ball
        planets.remove(p2)
        collided = False


def calca(target, kind, data):
    res = 0
    for planet in planets:
        if planet != target:
            if kind == 0:
                distance = math.sqrt((planet.pos[0] - data[0])**2 + (planet.pos[1] - target.pos[1])**2 + (planet.pos[2] - target.pos[2])**2)
            if kind == 1:
                distance = math.sqrt((planet.pos[0] - target.pos[0])**2 + (planet.pos[1] - data[1])**2 + (planet.pos[2] - target.pos[2])**2)
            if kind == 2:
                distance = math.sqrt((planet.pos[0] - target.pos[0])**2 + (planet.pos[1] - target.pos[1])**2 + (planet.pos[2] - data[2])**2)
            
            if distance < planet.shape.radius + target.shape.radius:
                collided = True
                mergingStars(target, planet)
                break
            F = G * planet.mass * target.mass / (distance**2)
            res = res + F * (planet.pos[kind] - data[kind]) / distance
    return res/target.mass


def ODE(planet, func, index, dt):
    # RK4th method
    data = [0, 0, 0]
    data[index] = planet.pos[index]     
    k1 = func(planet, index, data)
    data[index] = planet.pos[index] + planet.velocity[index] * dt/2 + 0.5 * k1 * (dt/2)**2
    k2 = func(planet, index, data)
    data[index] = planet.pos[index] + planet.velocity[index] * dt/2 + 0.5 * k2 * (dt/2)**2
    k3 = func(planet, index, data)
    data[index] = planet.pos[index] + planet.velocity[index] * dt + 0.5 * k3 * dt**2
    k4 = func(planet, index, data)
    a = (1/6)*(k1 + 2*k2 + 2*k3 + k4)
    planet.pos[index] = planet.pos[index] + planet.velocity[index] * dt + 0.5 * a * dt ** 2
    planet.velocity[index] = planet.velocity[index] + a * dt
    return planet.pos[index]



N = int(input("Bodies: "))
for i in range(1, N + 1, 1):
    x = rd.uniform(-1000e6, 1000e6)
    y = rd.uniform(-1000e6, 1000e6)
    z = rd.uniform(-1000e6, 1000e6)
    vx = rd.uniform(-170.55, 170.55)
    vy = rd.uniform(-170.55, 170.55)
    vz = rd.uniform(-170.55, 170.55)
    m = rd.uniform(1.0224e+22, 2.1469e+24)
    d = 600
    create_planet(x, y, z, vx, vy, vz, m, d)

print('=======================')
dt = float(input('time interval: '))
time = 10000
t = 0

while not collided:
    rate(60)
    for planet in planets:
        planet.shape.pos = vector(ODE(planet, calca, 0, dt), ODE(planet, calca, 1, dt), ODE(planet, calca, 2, dt))
    t = t + dt