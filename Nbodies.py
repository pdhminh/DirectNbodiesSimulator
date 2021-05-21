import random as rd
import math
import numpy as np
from vpython import *
scene2 = canvas(title='N bodies simulator', x=0, y=0, width=800, height=800, center=vector(0, 0, 0)) # background=vector(1, 1, 1
G = 6.67408e-11
collided = False
planets = []
# a = float(input('softing const: '))


class Planets(object):
    def __init__(self, xo, yo, zo, vox, voy, voz, mass, radius):
        self.velocity = [vox, voy, voz]
        self.mass = mass
        self.shape = simple_sphere(pos=vector(xo, yo, zo), radius=radius, color=vector(0, 1, 0))


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


def calca(target, kind, x, y, z):
    res = 0
    for planet in planets:
        if planet != target:
            distance = math.sqrt((planet.shape.pos.x - x)**2 + (planet.shape.pos.y - y)**2 + (planet.shape.pos.z - z)**2)
            if distance < planet.shape.radius + target.shape.radius:
                collided = True
                mergingStars(target, planet)
                break
                # return 0
            F = G * planet.mass * target.mass / (distance**2)
            if kind == 0:
                res = res + F * (planet.shape.pos.x - x) / distance
            elif kind == 1:
                res = res + F * (planet.shape.pos.y - y) / distance
            elif kind == 2:
                res = res + F * (planet.shape.pos.z - z) / distance
    return res/target.mass


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
        # RK4th method
        k1 = calca(planet, 0, x=planet.shape.pos.x, y=planet.shape.pos.y, z=planet.shape.pos.z)
        x2 = planet.shape.pos.x + planet.velocity[0] * dt/2 + 0.5 * k1 * (dt/2)**2
        k2 = calca(planet, 0, x=x2, y=planet.shape.pos.y, z=planet.shape.pos.z)
        x3 = planet.shape.pos.x + planet.velocity[0] * dt/2 + 0.5 * k2 * (dt/2)**2
        k3 = calca(planet, 0, x=x3, y=planet.shape.pos.y, z=planet.shape.pos.z)
        x4 = planet.shape.pos.x + planet.velocity[0] * dt + 0.5 * k3 * dt**2
        k4 = calca(planet, 0, x=x4, y=planet.shape.pos.y, z=planet.shape.pos.z)
        ax = (1/6)*(k1 + 2*k2 + 2*k3 + k4)
        planet.shape.pos.x = planet.shape.pos.x + planet.velocity[0]*dt + 0.5*ax*dt**2
        planet.velocity[0] = planet.velocity[0] + ax * dt

        k1 = calca(planet, 1, x=planet.shape.pos.x, y=planet.shape.pos.y, z=planet.shape.pos.z)
        y2 = planet.shape.pos.y + planet.velocity[1] * dt/2 + 0.5 * k1 * (dt/2)**2
        k2 = calca(planet, 1, x=planet.shape.pos.x, y=y2, z=planet.shape.pos.z)
        y3 = planet.shape.pos.y + planet.velocity[1] * dt/2 + 0.5 * k2 * (dt/2)**2
        k3 = calca(planet, 1, x=planet.shape.pos.x, y=y3, z=planet.shape.pos.z)
        y4 = planet.shape.pos.y + planet.velocity[1] * dt + 0.5 * k3 * dt**2
        k4 = calca(planet, 1, x=planet.shape.pos.x, y=y4, z=planet.shape.pos.z)
        ay = (1/6)*(k1 + 2*k2 + 2*k3 + k4)
        planet.shape.pos.y = planet.shape.pos.y + planet.velocity[1]*dt + 0.5*ay*dt**2
        planet.velocity[1] = planet.velocity[1] + ay * dt

        k1 = calca(planet, 2, x=planet.shape.pos.x, y=planet.shape.pos.y, z=planet.shape.pos.z)
        z2 = planet.shape.pos.z + planet.velocity[2] * dt / 2 + 0.5 * k1 * (dt / 2) ** 2
        k2 = calca(planet, 2, x=planet.shape.pos.x, y=planet.shape.pos.y, z=z2)
        z3 = planet.shape.pos.z + planet.velocity[2] * dt / 2 + 0.5 * k2 * (dt / 2) ** 2
        k3 = calca(planet, 2, x=planet.shape.pos.x, y=planet.shape.pos.y, z=z3)
        z4 = planet.shape.pos.z + planet.velocity[2] * dt + 0.5 * k3 * dt ** 2
        k4 = calca(planet, 2, x=planet.shape.pos.x, y=planet.shape.pos.y, z=z4)
        az = (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        planet.shape.pos.z = planet.shape.pos.z + planet.velocity[2] * dt + 0.5 * az * dt ** 2
        planet.velocity[2] = planet.velocity[2] + az * dt

    t = t + dt