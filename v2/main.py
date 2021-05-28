import math
import numpy as np
from vpython import *
import random as rd

planets = []
G = 6.67408e-11
dt = float(input('dt: '))
N = int(input('N: '))
d = 600
planets = []

class Planet():
    def __init__(self, x, y, vx, vy, m, r):
        self.pos = [x, y]
        self.vel = [vx, vy]
        self.mass = m
        self.radius = r 
        self.shape = simple_sphere(pos=vector(self.pos[0], self.pos[1], 0), color=vector(0,1,0), radius=self.radius)
    

    def ODE(self, index):
        data = [0, 0]
        data[index] = self.pos[index]
        k1 = self.calca(index, data)
        data[index] = self.pos[index] + self.vel[index] * dt / 2 + 0.5 * k1 * (dt / 2)**2
        k2 = self.calca(index, data)
        data[index] = self.pos[index] + self.vel[index] * dt / 2 + 0.5 * k2 * (dt / 2)**2
        k3 = self.calca(index, data)
        data[index] = self.pos[index] + self.vel[index] * dt + 0.5 * k3 * dt**2
        k4 = self.calca(index, data)
        a = (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        self.pos[index] = self.pos[index] + self.vel[index] * dt + 0.5 * a * dt ** 2
        self.vel[index] = self.vel[index] + a * dt


    def calca(self, index, data):
        res = 0
        for planet in planets:
            if planet != self:
                if index == 0:
                    distance = math.sqrt((planet.pos[0] - data[0])**2 + (planet.pos[1] - self.pos[1])**2)
                else:
                    distance = math.sqrt((planet.pos[0] - self.pos[0])**2 + (planet.pos[1] - data[1])**2)
                if distance < planet.shape.radius + self.shape.radius:
                    self.mergingStars(planet)
                F = G * self.mass * planet.mass / (distance**2)
                res = res + F * (planet.pos[index] - data[index]) / distance
        return res/self.mass


    def mergingStars(self, planet):
        global planets
        if self.mass > planet.mass:
            for i in range(0, 2, 1):
                self.vel[i] = (self.vel[i]*self.mass + planet.vel[i]*planet.mass)/(self.mass + planet.mass)
            self.mass = self.mass + planet.mass
            V = self.mass / d
            self.shape.radius = np.cbrt(V * (3/4) / np.pi)
            ball = planet.shape
            ball.visible = False
            del ball
            planets.remove(planet)

    
    def move(self):
        for i in range(0, 2, 1):
            self.ODE(i)
        self.shape.pos = vector(self.pos[0], self.pos[1], 0)


scene2 = canvas(title='N bodies simulator', x=0, y=0, width=800,
                height=800, center=vector(0, 0, 0))
for i in range(0, N, 1):
    # print(i)
    x = rd.uniform(-1000e6, 1000e6)
    y = rd.uniform(-1000e6, 1000e6)
    vx = rd.uniform(-170.55, 170.55)
    vy = rd.uniform(-170.55, 170.55)
    m = rd.uniform(1.0224e+22, 2.1469e+24)
    V = m / d
    r = np.cbrt(V * (3/4) / np.pi)
    planets.append(Planet(x, y, vx, vy, m, r))


while True:
    rate(60)
    for planet in planets:
        planet.move()