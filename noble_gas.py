import random as rnd
from scipy.constants import Boltzmann
import numpy as np
import matplotlib.pyplot as plt

size_of_box = 100


def two_d_distance(ball1, ball2):
    x_dist = abs(ball1.position[0] - ball2.position[0])
    if x_dist > size_of_box / 2:
        x_dist = size_of_box - x_dist

    y_dist = abs(ball1.position[1] - ball2.position[1])
    if y_dist > size_of_box / 2:
        y_dist = size_of_box - y_dist

    return (x_dist**2 + y_dist**2)**(1/2)


class Ball:
    def __init__(self, position, mass, velocity):
        self.position = position
        self.mass = mass
        self.velocity = velocity
        self.previous_velocity = velocity
        self.neighbours = []

    def add_neighbour(self, neighbour):
        self.neighbours.append(neighbour)
        neighbour.neighbours.append(self)


class System:
    def __init__(self, list_of_balls):
        self.list_of_balls = list_of_balls
        self.history = []


class Potential:
    def calc_force(self):
        pass


class Harmonic(Potential):
    def __init__(self, k, x0):
        self.k = k
        self.x0 = x0

    def calc_force_2d(self, ball1, ball2):
        dist = two_d_distance(ball1, ball2)
        force_abs = -1 * self.k * (dist - self.x0)
        force_x = np.sign(ball1.position[0] - ball2.position[0]) * force_abs\
                  * ball1.position[0] / dist
        force_y = np.sign(ball1.position[1] - ball2.position[1]) * force_abs\
                  * ball1.position[1] / dist

        return np.array([force_x, force_y])


class Langevin(Potential):
    def __init__(self, temperature):
        self.temperature = temperature

    def calc_force(self, gamma):
        return rnd.gauss(0, 1) * (2 * gamma * Boltzmann * self.temperature) ** (1/2)

    def calc_force_2d(self, gamma):
        force_x = self.calc_force(gamma)
        force_y = self.calc_force(gamma)
        return np.array([force_x, force_y])


class LennardJones(Potential):
    def __init__(self, epsilon, r0):
        self.epsilon = epsilon
        self.r0 = r0

    def calc_force_2d(self, ball1, ball2):
        dist = two_d_distance(ball1, ball2)
        force_abs = 12 * self.epsilon * (self.r0**12 / dist**13 - self.r0**6 / dist**7)
        force_x = np.sign(ball1.position[0] - ball2.position[0]) * force_abs \
                  * ball1.position[0] / dist
        if abs(ball1.position[0] - ball2.position[0]) > size_of_box / 2:
            force_x = -1 * force_x

        force_y = np.sign(ball1.position[1] - ball2.position[1]) * force_abs \
                  * ball1.position[1] / dist
        if abs(ball1.position[1] - ball2.position[1]) > size_of_box / 2:
            force_y = -1 * force_y
        return np.array([force_x, force_y])


class Algorithm:
    def update_pos(self, ball):
        pass

    def update_vel(self, ball):
        pass

    def resistance(self, ball, gamma):
        return -1 * (ball.previous_velocity * gamma)


class LeapFrog(Algorithm):
    def __init__(self, dt, system, potential, gamma):
        self.dt = dt
        self.system = system
        self.potential = potential
        self.gamma = gamma

    def calc_all_harmonic(self, ball, pot):
        force = np.array([0.0, 0.0])
        for other_ball in ball.neighbours:
            force += pot.calc_force_2d(ball, other_ball)
        return force

    def calc_all_lennardjones(self, ball, pot):
        force = 0.0
        for row in self.system.list_of_balls:
            for other_ball in row:
                if other_ball != ball:
                    force += pot.calc_force_2d(ball, other_ball)
        return force

    def calc_total_force(self, ball):
        force = np.array([0.0, 0.0])
        for pot in self.potential:
            if isinstance(pot, Harmonic):
                force += self.calc_all_harmonic(ball, pot)
            if isinstance(pot, Langevin):
                force += pot.calc_force_2d(self.gamma)
            if isinstance(pot, LennardJones):
                force += self.calc_all_lennardjones(ball, pot)
            force += self.resistance(ball, self.gamma)
        return force

    def update_vel(self, ball):
        ai = self.calc_total_force(ball) / ball.mass
        vi_12 = ball.previous_velocity + self.dt * ai
        ball.velocity = vi_12

    def update_pos(self, ball):
        xi_1 = ball.position + ball.velocity * self.dt
        ball.position = np.remainder(xi_1, 100)


class Simulation:
    def __init__(self, system, dt, steps, potentials, gamma, algorithm, n):
        self.system = system
        self.dt = dt
        self.steps = steps
        self.potentials = potentials
        self.gamma = gamma
        self.algorithm = algorithm(self.dt, self.system, self.potentials, self.gamma)
        self.history = []
        self.n = n

    def run(self):
        print("Running simulation...\n")
        for x in range(self.steps):
            #print(f'\nStep {x + 1}: ')
            for row in self.system.list_of_balls:
                for ball in row:
                    self.algorithm.update_vel(ball)
            for row in self.system.list_of_balls:
                for ball in row:
                    ball.previous_velocity = ball.velocity
                    self.algorithm.update_pos(ball)
            #for i, ball in enumerate(self.system.list_of_balls):
            #    print(f'Position of ball {i + 1}: {ball.position}')
            self.system.history.append(two_d_distance(self.system.list_of_balls[0][0], self.system.list_of_balls[0][1]))


def square(n):
    network = []
    for x in range(n):
        network.append([Ball(np.array([y * 10.0, (n - x - 1) * 10.0]), 1, np.array([0.0, 0.0]))\
                        for y in range(n)])
    return network


def main():
    # n = int(input('Wpisz jakiej długości ma być krawędź sieci: '))
    #
    # list_of_balls = square(n)
    # list_of_balls[n-1][0].position = np.array([-1.0, -1.0])

    # for x in range(n):
    #     for y in range(n):
    #         if x + 1 != n:
    #             list_of_balls[x][y].add_neighbour(list_of_balls[x + 1][y])
    #         if y + 1 != n:
    #             list_of_balls[x][y].add_neighbour(list_of_balls[x][y + 1])
    n=2

    list_of_balls = [[Ball(np.array([91.0, 20.0]), 1, np.array([0.0, 0.0])),\
                     Ball(np.array([5.0, 20.0]), 1, np.array([0.0, 0.0]))]]
    S = System(list_of_balls)
    #Harmonic(k=10, x0=10), Langevin(temperature=10)
    sim = Simulation(system=S, dt=0.01, steps=2500, gamma=0.3, potentials=[LennardJones(10, 10)],\
                     algorithm=LeapFrog, n=n)
    S.history = []
    sim.run()
    plt.plot(S.history)
    plt.savefig('plot.png')


if __name__ == '__main__':
    main()
