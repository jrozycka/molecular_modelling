import numpy as np
import random
import matplotlib.pyplot as plt
import copy
from math import exp


def random_dir():
    directions = [[1, 0], [0, 1], [-1, 0], [0, -1]]
    return np.array(random.choice(directions))


def four_dirs():
    directions = [[1, 0], [0, 1], [-1, 0], [0, -1]]
    return [np.array(dir) for dir in directions]


class Mer:
    def __init__(self, position, id):
        self.position = position
        self.neighbours = [None, None]
        self.id = id
        self.spin = random.choice([-1, 1])

    def add_neighbour(self, neighbour):
        self.neighbours[1] = neighbour
        neighbour.neighbours[0] = self


class System:
    def __init__(self, n, m, temp):
        self.n = n #rozmiar sieci
        self.m = m #długość polimera
        self.temp = temp
        self.lattice = np.zeros([n, n])
        self.polymer = []

    def init_polymer(self):
        while True:
            self.lattice = np.zeros([self.n, self.n])
            self.polymer = []

            position = np.array([self.n // 2, self.n // 2])
            self.polymer.append(Mer(position, 1))
            self.lattice[tuple(position)] = 1

            for i in range(self.m - 1):
                good_structure = False

                while not good_structure:
                    dirs = four_dirs()
                    random.shuffle(dirs)
                    for dir in dirs:
                        new_position = position + dir

                        if self.lattice[tuple(new_position)] == 0: #można zrobić kolejny krok
                            self.polymer.append(Mer(new_position, i + 2))
                            self.lattice[tuple(new_position)] = i + 2
                            old_mer = self.polymer[i]
                            neighbour = self.polymer[i + 1]
                            old_mer.neighbours[1] = neighbour
                            neighbour.neighbours[0] = old_mer
                            position = new_position
                            good_structure = True
                            break

                    if not good_structure:
                        break
            break


class Move:
    def find_move(self, lattice, polimer):
        pass

    def update_position(self, lattice, polimer):
        pass


class MoveEnd(Move):
    def __init__(self):
        pass

    def move_one_end(self, mer, lattice):
        dirs = four_dirs()
        random.shuffle(dirs)

        for dir in dirs:
            move = mer.position + dir
            if lattice[tuple(move)] == 0:
                return move

        return None

    def find_move(self, lattice, polimer):
        ends = [polimer[0], polimer[-1]]
        end1 = ends.pop(random.randint(0, 1))
        end1_neighbour = next(item for item in end1.neighbours if item is not None)
        move = self.move_one_end(end1_neighbour, lattice)

        if move is not None:
            return end1, move
        else:
            end2 = ends[0]
            end2_neighbour = next(item for item in end2.neighbours if item is not None)
            move = self.move_one_end(end2_neighbour, lattice)

        if move is not None:
            return end2, move
        else:
            return None


class MoveCorner(Move):
    def __init__(self):
        pass

    def check_mer(self, mer_id, polimer, lattice):
        left_neighbour_pos = polimer[mer_id - 2].position
        right_neighbour_pos = polimer[mer_id].position
        mer_pos = polimer[mer_id - 1].position
        corner = right_neighbour_pos - left_neighbour_pos

        if abs(corner[0]) == 1 and abs(corner[1]) == 1:
            move = np.array([0, 0])
            left_coord = np.nonzero(mer_pos - left_neighbour_pos)
            move[left_coord] = left_neighbour_pos[left_coord]

            right_coord = np.nonzero(mer_pos - right_neighbour_pos)
            move[right_coord] = right_neighbour_pos[right_coord]

            if lattice[tuple(move)] == 0:
                return move
        return None

    def find_move(self, lattice, polimer):
        available_mers = [x for x in range(1, len(polimer) - 1)]
        random.shuffle(available_mers)

        for mer in available_mers:
            move = self.check_mer(mer + 1, polimer, lattice)
            if move is not None:
                return polimer[mer], move
        return None


class MoveLoop(Move):
    def __init__(self):
        pass

    def check_mer(self, mer_id, polimer, lattice):
        if polimer[mer_id - 1].position[0] == polimer[mer_id + 2].position[0]\
                and abs(polimer[mer_id - 1].position[1] - polimer[mer_id + 2].position[1]) == 1:
            #ta sama x współrzędna
            new_poz_x1 = polimer[mer_id - 1].position + np.array([1, 0])
            new_poz_x2 = polimer[mer_id + 2].position + np.array([1, 0])
            if lattice[tuple(new_poz_x1)] == 0 and lattice[tuple(new_poz_x2)] == 0:
                return new_poz_x1, new_poz_x2

            new_poz_x1 = polimer[mer_id - 1].position + np.array([-1, 0])
            new_poz_x2 = polimer[mer_id + 2].position + np.array([-1, 0])
            if lattice[tuple(new_poz_x1)] == 0 and lattice[tuple(new_poz_x2)] == 0:
                return new_poz_x1, new_poz_x2

        elif polimer[mer_id - 1].position[1] == polimer[mer_id + 2].position[1]\
                and abs(polimer[mer_id - 1].position[0] - polimer[mer_id + 2].position[0]) == 1:
            #ta sama y współrzędna
            new_poz_x1 = polimer[mer_id - 1].position + np.array([0, 1])
            new_poz_x2 = polimer[mer_id + 2].position + np.array([0, 1])
            if lattice[tuple(new_poz_x1)] == 0 and lattice[tuple(new_poz_x2)] == 0:
                return new_poz_x1, new_poz_x2

            new_poz_x1 = polimer[mer_id - 1].position + np.array([0, -1])
            new_poz_x2 = polimer[mer_id + 2].position + np.array([0, -1])
            if lattice[tuple(new_poz_x1)] == 0 and lattice[tuple(new_poz_x2)] == 0:
                return new_poz_x1, new_poz_x2

        else:
            return None

    def find_move(self, lattice, polimer):
        available_mers = [x for x in range(1, len(polimer) - 2)]
        random.shuffle(available_mers)

        for mer in available_mers:
            move = self.check_mer(mer, polimer, lattice)
            if move is not None:
                return (polimer[mer], polimer[mer + 1], move[0], move[1])

        return None


class MoveSnake(Move):
    def __init__(self):
        pass

    def check_end(self, end, lattice):
        dirs = four_dirs()
        random.shuffle(dirs)

        for dir in dirs:
            move = end.position + dir
            if lattice[tuple(move)] == 0:
                return end, move

        return None

    def find_move(self, lattice, polimer):
        ends = [polimer[0], polimer[-1]]
        random.shuffle(ends)

        for end in ends:
            move = self.check_end(end, lattice)
            if move is not None:
                return move

        return None


class Simulation:
    def __init__(self, lattice, polimer, temp, steps, mc_steps):
        self.lattice = lattice
        self.polimer = polimer
        self.temp = temp
        self.steps = steps
        self.mc_steps = mc_steps
        pass

    def calc_energy(self, polimer):
        energy = 0
        lattice_of_spins = np.zeros([len(self.lattice), len(self.lattice)])
        for mer in polimer:
            lattice_of_spins[tuple(mer.position)] = mer.spin
        for mer in polimer:
            e = 0
            for dir in four_dirs():
                if lattice_of_spins[tuple(mer.position + dir)] != 0:
                    e += (lattice_of_spins[tuple(mer.position + dir)] * mer.spin)
            energy += e
        return energy / 2

    def move_spin(self, polimer):
        polimer_copy = copy.deepcopy(polimer)
        mer = random.randrange(0, len(polimer))
        polimer_copy[mer].spin = polimer[mer].spin * (-1)
        return polimer_copy

    def monte_carlo(self, new_polimer, steps):
        energies = np.zeros(steps)
        for x in range(steps):
            energy_old = self.calc_energy(new_polimer)
            poli_copy = self.move_spin(new_polimer)
            energy_new = self.calc_energy(poli_copy)
            if self.calc_prob(energy_old, energy_new):
                new_polimer = poli_copy
                energy_old = energy_new
            energies[x] = energy_old
        return np.mean(energies[int(steps / 2):])

    def choose_move(self):
        move = random.choice([MoveEnd(), MoveCorner(), MoveLoop(), MoveSnake()])
        mer_move_pair = move.find_move(self.lattice, self.polimer)
        return mer_move_pair, move

    def calc_prob(self, energy_old, energy_new):
        delta = energy_new - energy_old
        p = exp(-1 * delta / self.temp)
        if random.random() > p:
            return False
        else:
            return True

    def step(self):
        mer_move_pair, move = self.choose_move()
        if not mer_move_pair:
            energy = self.calc_energy(self.polimer)
            return energy

        elif isinstance(move, MoveLoop):
            mer1, mer2, move1, move2 = mer_move_pair
            pos_old1 = mer1.position
            pos_old2 = mer2.position
            pos_new1 = move1
            pos_new2 = move2
            lattice_new = copy.deepcopy(self.lattice)
            lattice_new[tuple(pos_old1)] = 0
            lattice_new[tuple(pos_old2)] = 0
            lattice_new[tuple(pos_new1)] = mer1.id
            lattice_new[tuple(pos_new2)] = mer2.id
            polimer_new = copy.deepcopy(self.polimer)
            polimer_new[mer1.id - 1].position = pos_new1
            polimer_new[mer2.id - 1].position = pos_new2

        elif isinstance(move, MoveSnake):
            mer, move = mer_move_pair
            polimer_new = copy.deepcopy(self.polimer)
            lattice_new = copy.deepcopy(self.lattice)

            if mer.id == 1:
                lattice_new[tuple(self.polimer[-1].position)] = 0
                for x in range(len(self.polimer) - 1, 0, -1):
                    polimer_new[x].position = polimer_new[x - 1].position
                    lattice_new[tuple(polimer_new[x - 1].position)] = x + 1

                polimer_new[0].position = move
                lattice_new[tuple(move)] = 1

            else:
                lattice_new[tuple(self.polimer[0].position)] = 0
                for x in range(len(self.polimer) - 1):
                    polimer_new[x].position = polimer_new[x + 1].position
                    lattice_new[tuple(polimer_new[x + 1].position)] = x + 1

                polimer_new[-1].position = move
                lattice_new[tuple(move)] = mer.id

        else:
            mer, move = mer_move_pair

            pos_old = mer.position
            pos_new = move
            lattice_new = copy.deepcopy(self.lattice)
            lattice_new[tuple(pos_old)] = 0
            lattice_new[tuple(pos_new)] = mer.id
            polimer_new = copy.deepcopy(self.polimer)
            polimer_new[mer.id - 1].position = pos_new

        energy_old = self.calc_energy(self.polimer)
        energy_new = self.monte_carlo(polimer_new, self.mc_steps)

        if energy_new <= energy_old:
            self.lattice = lattice_new
            self.polimer = polimer_new
            energy_old = energy_new
        else:
            good_energy = self.calc_prob(energy_old, energy_new)
            if good_energy:
                self.lattice = lattice_new
                self.polimer = polimer_new
                energy_old = energy_new
        return energy_old

    def simulation(self):
        energy_list = np.empty(self.steps)
        for x in range(self.steps):
            energy_list[x] = self.step()

        return energy_list


def simulation_in_temperature(lattice, polymer_length, temp, steps, mc_steps):
    S = System(lattice, polymer_length, temp)
    S.init_polymer()
    sim = Simulation(S.lattice, S.polymer, S.temp, steps, mc_steps)
    energy_list = sim.simulation()
    mean = np.mean(energy_list[(steps // 2):])
    sd = np.std(energy_list[(steps // 2):])

    return mean, sd


def plot_mean_sd(temp_range, lattice, polimer_size, steps, monte_carlo_steps):
    mean_list = np.empty(len(temp_range))
    sd_list = np.empty(len(temp_range))

    for x in range(len(temp_range)):
        print('Building model for temperature ', temp_range[x])
        mean, sd = simulation_in_temperature(lattice, polimer_size, temp_range[x], steps, monte_carlo_steps)
        mean_list[x] = mean
        sd_list[x] = sd

    plt.plot(temp_range, mean_list, label='mean energy')
    plt.plot(temp_range, sd_list, label='energy sd')
    plt.xlabel('temperature')
    plt.legend()
    plt.savefig('energy.png')


def main():
    temp = np.linspace(0.1, 10, 20)
    plot_mean_sd(temp, 100, 20, 500, 100) #temp_range, lattice, polimer_size, steps, monte_carlo_steps


if __name__ == "__main__":
    main()
