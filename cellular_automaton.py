import random
import numpy as np
import copy


def generate_cells(n):
    start = np.zeros(n)
    for x in range(n):
        start[x] = random.randint(0, 1)
    return start


def calc_score(cells):
    e = 0
    n = len(cells)
    for i in range(n):
        if cells[i] != cells[(i + 1) % n]:
            e += 1
    return e


def generate_rule():
    rule = bin(random.randint(0, 2**8 - 1))[2:].zfill(8)
    return np.array([int(x) for x in rule])


def generate_population(size):
    population = np.zeros([size, 8])
    for x in range(size):
        population[x] = generate_rule()
    return population


def calc_new_cell(rule, three):
    three_to_decimal = int(str(three[0])[0] + str(three[1])[0] + str(three[2])[0], 2)
    return rule[three_to_decimal]


def cells_printable(cell):
    printable = 'l: '
    for x in cell:
        if x == 0:
            printable += ' '
        else:
            printable += u"\u2588"
    print(printable)


def test_rule(rule, steps, n, printing=False):
    cell = generate_cells(n)
    for x in range(steps):
        new_cells = np.zeros(n)
        for i in range(n):
            new_cells[i] = calc_new_cell(rule, [cell[i-1], cell[i], cell[(i + 1) % n]])
        cell = new_cells
        if printing:
            cells_printable(cell)
    return calc_score(cell)


def point_mutation(cells):
    new_cells = copy.deepcopy(cells)
    where = random.randint(0, len(cells) - 1)
    if new_cells[where] == 0:
        new_cells[where] = 1
    else:
        new_cells[where] = 0
    return new_cells


def offspring(pair_of_cells):
    """input: list of two cells"""
    child = np.zeros(8)
    for x in range(8):
        child[x] = pair_of_cells[random.choice([0, 1])][x]
    return child


def genetic_algorithm(pop_size, cell_length, steps):
    population = generate_population(pop_size)
    scores = np.zeros(pop_size)

    for step in range(steps):
        for x in range(pop_size):
            scores[x] = test_rule(population[x], steps, cell_length)

        print(f'In step {step} mean score is {np.mean(scores)}')
        best_scores = scores > np.median(scores)
        # don't have to choose if many have the same score
        new_population = np.zeros([pop_size, 8])
        new_population[:sum(best_scores)] = population[best_scores]

        for x in range(sum(best_scores), pop_size):
            # supplying new populations
            if x < pop_size * (5/8):
                random_good = new_population[random.randrange(sum(best_scores))]
                mutant = point_mutation(random_good)
                new_population[x] = mutant

            elif x < pop_size * (6/8):
                random_pair = new_population[random.sample(range(sum(best_scores)), 2)]
                child = offspring(random_pair)
                new_population[x] = child

            else:
                new_population[x] = generate_rule()

        population = new_population

    best = np.argmax(scores)
    return population[best]


test_rule([0, 1, 1, 1, 0, 1, 1, 0], 200, 100, printing=True)

solution = genetic_algorithm(50, 30, 200)
if solution is not None:
    print(solution)
    test_rule(solution, 200, 100, printing=True)
