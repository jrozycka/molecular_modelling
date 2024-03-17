import random
from math import exp
import numpy as np


class Lattice:
    """
    Represents a lattice of cells with spins of 1 or -1.

    Attributes:
        size (int): The size of the lattice (size x size).
        temperature (float): The temperature parameter controlling the acceptance probability.
        lattice (numpy.ndarray): The lattice representing cells with spins.
        energies (numpy.ndarray): Energies of each cell in the lattice.
        total_energy (float): Total energy of the lattice.
    """

    def __init__(self, size, temperature):
        """
        Initializes a lattice of size x size with random spins.

        Args:
            size (int): The size of the lattice (size x size).
            temperature (float): The temperature parameter controlling the acceptance probability.
        """
        self.size = size
        self.temperature = temperature
        self.lattice = np.random.choice([-1, 1], size=(size, size))
        self.energies = np.zeros((size, size))
        self.total_energy = 0
        self.calculate_cell_energies()
        self.calculate_total_energy()

    def calculate_cell_energy(self, i, j):
        """
        Calculates the energy of a given cell.

        Args:
            i (int): Row index of the cell.
            j (int): Column index of the cell.

        Returns:
            float: Energy of the cell.
        """
        energy = 0
        spin = self.lattice[i, j]
        energy -= spin * self.lattice[(i + 1) % self.size, j]  # Bottom neighbor
        energy -= spin * self.lattice[(i - 1) % self.size, j]  # Top neighbor
        energy -= spin * self.lattice[i, (j + 1) % self.size]  # Right neighbor
        energy -= spin * self.lattice[i, (j - 1) % self.size]  # Left neighbor
        return energy

    def calculate_cell_energies(self):
        """
        Calculates energies of all cells in the lattice and stores them.
        """
        for i in range(self.size):
            for j in range(self.size):
                self.energies[i, j] = self.calculate_cell_energy(i, j)

    def calculate_total_energy(self):
        """
        Calculates the total energy of the lattice.
        """
        self.total_energy = np.sum(self.energies) / 2

    def metropolis_hastings_step(self):
        """
        Attempt a random spin flip using the Metropolis-Hastings algorithm.

        Returns:
            bool: True if the spin flip was accepted, False otherwise.
        """
        # Choose a random cell
        i, j = random.randint(0, self.size - 1), random.randint(0, self.size - 1)

        # Calculate the energy of the current configuration
        current_energy = self.calculate_cell_energy(i, j)

        # Flip the spin
        self.lattice[i, j] *= -1

        # Calculate the energy of the new configuration
        new_energy = self.calculate_cell_energy(i, j)

        # Calculate the energy difference
        delta_energy = new_energy - current_energy

        # Calculate acceptance probability
        if delta_energy <= 0:
            accept_prob = 1.0
        else:
            accept_prob = exp(-delta_energy / self.temperature)

        # Accept or reject the spin flip based on acceptance probability
        if random.random() < accept_prob:
            # Update energies
            self.update_energies_after_flip(i, j)
            return True
        else:
            return False

    def update_energies_after_flip(self, i, j):
        """
        Update energies after flipping the spin at cell (i, j).
        """
        # Update energy of the flipped cell
        self.energies[i, j] = self.calculate_cell_energy(i, j)

        # Update energy of neighboring cells
        for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            ni, nj = (i + di) % self.size, (j + dj) % self.size
            self.energies[ni, nj] = self.calculate_cell_energy(ni, nj)

        # Update total energy
        self.calculate_total_energy()


class Simulation:
    """
    Represents a simulation of a lattice with spins of 1 or -1.

    Attributes:
        temperature (float): The temperature parameter controlling the acceptance probability.
        size (int): The size of the lattice (size x size).
        steps (int): Number of steps to run the simulation.
        lattice (Lattice): The lattice object representing the system.
        energies (list): Energies of the lattice at each step.
    """

    def __init__(self, temperature, size, steps):
        """
        Initializes a simulation of a lattice with spins of 1 or -1.

        Args:
            temperature (float): The temperature parameter controlling the acceptance probability.
            size (int): The size of the lattice (size x size).
            steps (int): Number of steps to run the simulation.
        """
        self.temperature = temperature
        self.size = size
        self.steps = steps
        self.lattice = Lattice(size, temperature)
        self.energies = []

    def run_simulation(self):
        """
        Run the simulation for the specified number of steps.
        """
        for _ in range(self.steps):
            # Perform Metropolis-Hastings step
            self.lattice.metropolis_hastings_step()
            # Store total energy
            self.energies.append(self.lattice.total_energy)
