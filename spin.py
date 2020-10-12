import numpy as np
import random
import copy
from math import exp
import matplotlib.pyplot as plt


def create_lattice():
    siec = np.zeros([m, m])
    for x in range(m):
        for y in range(m):
            siec[x][y] = random.choice([-1, 1])
    return siec


def calc_cell_energy(siec, x, y):
    e = 0
    cell = siec[x][y]
    e += cell * siec[(x + 1) % m][y]
    e += cell * siec[x][(y + 1) % m]
    e += cell * siec[x - 1][y]
    e += cell * siec[x][y - 1]
    return e


def calc_energy(siec):
    e = 0
    for x in range(m):
        for y in range(m):
            if (x % 2 == 0 and y % 2 == 0) or (x % 2 == 1 and y % 2 == 1):
                e += calc_cell_energy(siec, x, y)
    return e


def calc_prob(energy_old, energy_new, temp):
    delta = energy_new - energy_old
    p = exp(-1 * delta / temp)
    if random.random() > p:
        return False
    else:
        return True


def change_spin(siec, temp):
    siec_new = copy.deepcopy(siec)
    x = random.randint(0, m - 1)
    y = random.randint(0, m - 1)
    siec_new[x, y] = siec[x, y] * -1

    e_old = calc_cell_energy(siec, x, y)
    e_new = calc_cell_energy(siec_new, x, y)

    if e_new <= e_old:
        siec = siec_new
    else:
        change = calc_prob(e_old, e_new, temp)
        if change:
            siec = siec_new
        else:
            e_new = e_old

    return siec, e_new - e_old


def monte_carlo(siec, temp, steps):
    energy_all = np.zeros(steps)
    energy = calc_energy(siec)

    for x in range(steps):
        siec_new, e = change_spin(siec, temp)
        energy += e
        siec = siec_new
        energy_all[x] = energy
    return np.mean(energy_all[(steps // 2):]), np.std(energy_all[(steps // 2):])


def main():
    m = 200 # size of lattice, must be even
    steps = 100000
    temp_range = np.linspace(0.1, 20, 200)
    energy = np.zeros(len(temp_range))
    std = np.zeros(len(temp_range))

    for x in range(len(temp_range)):
        print('calculating for temperature ', temp_range[x])
        siec = create_lattice()
        e, sd = monte_carlo(siec, temp_range[x], steps)
        energy[x] = e
        std[x] = sd

    plt.plot(energy, label='energy')
    plt.plot(std, label='sd')
    plt.savefig('spin.png')


if __name__ == "__main__":
    main()
