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

    def add_neighbour(self, neighbour):
        self.neighbours[1] = neighbour
        neighbour.neighbours[0] = self


class System:
    def __init__(self, n, m, temp):
        self.n = n # size of lattice
        self.m = m # polimer length
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

                        if self.lattice[tuple(new_position)] == 0:
                            # next step is possible
                            self.polymer.append(Mer(new_position, i + 2))
                            self.lattice[tuple(new_position)] = i + 2
                            old_mer = self.polymer[i]
                            neighbour = self.polymer[i + 1]
                            old_mer.neighbours[1] = neighbour
                            neighbour.neighbours[0] = old_mer
                            position = new_position
                            good_structure = True
                            break
                        # if it is not, try again

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
        # id is location in polymer + 1
        right_neighbour_pos = polimer[mer_id].position
        mer_pos = polimer[mer_id - 1].position
        corner = right_neighbour_pos - left_neighbour_pos

        if abs(corner[0]) == 1 and abs(corner[1]) == 1:
            # corner
            move = np.array([0, 0])
            # which coordinate should be taken from the left neighbour
            left_coord = np.nonzero(mer_pos - left_neighbour_pos)
            # which one from the right one
            move[left_coord] = left_neighbour_pos[left_coord]

            right_coord = np.nonzero(mer_pos - right_neighbour_pos)
            move[right_coord] = right_neighbour_pos[right_coord]

            if lattice[tuple(move)] == 0:
                # empty space
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
            # the same x coordinate
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
            # same y coorinate
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
    def __init__(self, lattice, polimer, temp, steps):
        self.lattice = lattice
        self.polimer = polimer
        self.temp = temp
        self.steps = steps
        pass

    def calc_energy(self, lattice, polimer):
        energy = 0
        for mer in polimer:
            e = 0
            for dir in four_dirs():
                if lattice[tuple(mer.position + dir)] > 0:
                    e -= 1
            energy += e + 2
        return energy / 2 - 1 # because start and end only have 1 neighbour and each neighbouring is counted twice

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
            energy = self.calc_energy(self.lattice, self.polimer)
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

        energy_old = self.calc_energy(self.lattice, self.polimer)
        energy_new = self.calc_energy(lattice_new, polimer_new)

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


def simulation_in_temperature(lattice, polymer_length, temp, steps):
    S = System(lattice, polymer_length, temp)
    S.init_polymer()
    sim = Simulation(S.lattice, S.polymer, S.temp, steps)
    energy_list = sim.simulation()
    mean = np.mean(energy_list[(steps // 2):])
    sd = np.std(energy_list[(steps // 2):])
    return mean, sd


def plot_mean_sd(temp_range):
    mean_list = np.empty(len(temp_range))
    sd_list = np.empty(len(temp_range))

    for x in range(len(temp_range)):
        print('budowanie modelu dla temperatury ', temp_range[x])
        mean, sd = simulation_in_temperature(100, 20, temp_range[x], 500)
        mean_list[x] = mean
        sd_list[x] = sd

    plt.plot(temp_range, mean_list, label='średnia energia')
    plt.plot(temp_range, sd_list, label='sd energii')
    plt.xlabel('temperatura')
    plt.legend()
    plt.show()


temp = np.linspace(0.1, 10)
plot_mean_sd(temp)
