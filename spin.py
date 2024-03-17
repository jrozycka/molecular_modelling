import numpy as np
import random
import copy
from math import exp
import matplotlib.pyplot as plt


def create_lattice(size):
    grid = np.zeros([size, size])
    for x in range(size):
        for y in range(size):
            grid[x][y] = random.choice([-1, 1])
    return grid


def calc_cell_energy(grid, x, y):
    e = 0
    size = grid.shape[0]
    cell = grid[x][y]
    e += cell * grid[(x + 1) % size][y]
    e += cell * grid[x][(y + 1) % size]
    e += cell * grid[x - 1][y]
    e += cell * grid[x][y - 1]
    return e


def calc_energy(grid):
    e = 0
    size = grid.shape[0]
    for x in range(size):
        for y in range(size):
            if (x % 2 == 0 and y % 2 == 0) or (x % 2 == 1 and y % 2 == 1):
                e += calc_cell_energy(grid, x, y)
    return e


def calc_prob(energy_old, energy_new, temp):
    delta = energy_new - energy_old
    p = exp(-delta / temp)
    p = min(p, 1.0)  # Ensure probability doesn't exceed 1
    if random.random() > p:
        return False
    else:
        return True


def change_spin(grid, temp):
    new_grid = copy.deepcopy(grid)
    size = new_grid.shape[0]
    x = random.randint(0, size - 1)
    y = random.randint(0, size - 1)
    new_grid[x, y] = grid[x, y] * -1

    e_old = calc_cell_energy(grid, x, y)
    e_new = calc_cell_energy(new_grid, x, y)

    if e_new <= e_old:
        grid = new_grid
    else:
        change = calc_prob(e_old, e_new, temp)
        if change:
            grid = new_grid
        else:
            e_new = e_old

    return grid, e_new - e_old


def monte_carlo(grid, temp, steps):
    energy_all = np.zeros(steps)
    energy = calc_energy(grid)

    for x in range(steps):
        new_grid, e = change_spin(grid, temp)
        energy += e
        grid = new_grid
        energy_all[x] = energy
    return np.mean(energy_all[(steps // 2) :]), np.std(energy_all[(steps // 2) :])


def main():
    m = 200  # size of lattice, must be even
    steps = 10000
    temp_range = np.linspace(0.1, 20, 20)
    energy = np.zeros(len(temp_range))
    std = np.zeros(len(temp_range))

    for x in range(len(temp_range)):
        print("calculating for temperature ", temp_range[x])
        grid = create_lattice(m)
        e, sd = monte_carlo(grid, temp_range[x], steps)
        energy[x] = e
        std[x] = sd

    plt.plot(energy, label="energy")
    plt.plot(std, label="sd")
    plt.savefig("size.png")


if __name__ == "__main__":
    main()
