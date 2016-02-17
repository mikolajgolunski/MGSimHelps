import math
import itertools
import os.path
from collections import namedtuple, Counter
from enum import Enum
import numpy as np
from tqdm import tqdm

import MGReadFile
import MGSaveFile


class FileType(Enum):
    auto = 1
    lammpstrj = 2
    lammps_data = 3
    xyz = 4


class LammpsDataStyle(Enum):
    charge = 1


class ControlDict(Enum):
    lammps_data_style = 1
    lammpstrj_properties = 2
    lammpstrj_items = 3


class LammpstrjProp(Enum):
    id = 1
    element = 2
    type = 3
    name = 4
    x = 5
    y = 6
    z = 7
    charge = 8
    mass = 9
    vx = 10
    vy = 11
    vz = 12
    kinetic_energy = 13
    potential_energy = 14


class LammpstrjItem(Enum):
    timestep = 1
    time = 2
    number_atoms = 3
    bounds = 4
    atoms = 5


class StopError(Exception):
    def __init__(self, value="Stopped."):
        self.value = value

    def __str__(self):
        return repr(self.value)


class Universe:
    """Head class."""
    def __init__(self):
        self.systems = []
        self.reax = []

    def readFile(self, file_path, file_type=FileType.auto, control_dict=None):
        """Read file of specified type and save its content to the AtomsSystem.

        :param file_path: [string] - path to the file
        :param file_type: [string] - [default "auto"] type of file:
            "auto" - automatically find file_type using file's extension
            "lammpstrj" - LAMMPS trajectory file
            "lammps_data" - file read by LAMMPS's function read_data
            "xyz" - xyz file
        :param control_dict: [dict] - [Default: {}] dictionary consisting of additional controls that may be needed:
            "lammps_data_type":
                "charge"
        """
        print("Starting 'readFile' procedure.")
        if control_dict is None:
            control_dict = {}
        if file_type == FileType.auto:
            print("Detecting filetype automaticaly: ", end="")
            file_ext = os.path.splitext(file_path)[1][1:].strip().lower()
            if file_ext == "lammpstrj":
                file_type = FileType.lammpstrj
            elif file_ext == "xyz":
                file_type = FileType.xyz
            else:
                raise NotImplementedError(
                        "No implementation for the " + file_ext + " extension. Please set file_type explicitly."
                )
            print(file_type.name)

        if file_type == FileType.lammpstrj:
            system = (system for system in MGReadFile.readerLammpsFile(file_path))
        elif file_type == FileType.lammps_data:
            system = (system for system in MGReadFile.readerLammpsDataFile(file_path, control_dict))
        else:
            raise NotImplementedError("No implementation for the " + str(file_type.name) + " file_type.")
        for s in system:
            print("\nSuccessfully read system from " + file_path + " as a " + str(file_type.name) + " file.")
            print("----------")
            yield s

    # noinspection PyTypeChecker
    def readReax(self, file_path):
        print("Reading ReaxFF parameters file.")
        ReaxGeneral = namedtuple("ReaxGeneral", ["p_boc1", "p_boc2"])
        ReaxAtoms = namedtuple("ReaxAtoms", ["name", "ro_sigma", "Val", "mass", "ro_pi", "ro_pipi", "p_boc4", "p_boc3",
                                             "p_boc5", "Val_boc"])
        ReaxBonds = namedtuple("ReaxBonds", ["at1", "at2", "De_sigma", "De_pi", "De_pipi", "p_be1", "p_bo5", "p_bo6",
                                             "p_be2", "p_bo3", "p_bo4", "p_bo1", "p_bo2"])
        ReaxOffDiagonals = namedtuple("ReaxOffDiagonals", ["at1", "at2", "ro_sigma", "ro_pi", "ro_pipi"])
        with open(file_path, "r") as file:
            temp_list = []
            for i in range(2):
                file.readline()
            for i in range(2):
                temp_list.append(float(file.readline().strip().partition(" !")[0]))
            reax_general = ReaxGeneral._make(temp_list)
            for i in range(37):
                file.readline()
            nr = int(file.readline().strip().split()[0])
            for i in range(3):
                file.readline()
            reax_atoms = []
            temp_line = file.readline().strip().split()
            for at in range(nr):
                temp_list = [temp_line[0]]
                temp_list.extend([float(temp_line[num]) for num in (1, 2, 3, 7)])

                file.readline()

                temp_line = file.readline().strip().split()
                temp_list.extend([float(temp_line[num]) for num in (0, 3, 4, 5)])

                temp_line = file.readline().strip().split()
                temp_list.extend([float(temp_line[num]) for num in (3, )])

                reax_atoms.append(ReaxAtoms._make(temp_list))

                stop = False
                while True:
                    temp_line = file.readline().strip().split()
                    try:
                        float(temp_line[0])
                        float(temp_line[-1])
                    except ValueError:
                        stop = True
                        break
                if stop:
                    continue
            nr = int(temp_line[0])
            file.readline()
            reax_bonds = []
            for bond in range(nr):
                temp_line = file.readline().strip().split()
                temp_list = [int(temp_line[0]), int(temp_line[1])]
                temp_list.extend([float(temp_line[num]) for num in (2, 3, 4, 5, 6, 8)])

                temp_line = file.readline().strip().split()
                temp_list.extend([float(temp_line[num]) for num in (0, 1, 2, 4, 5)])

                reax_bonds.append(ReaxBonds._make(temp_list))
            nr = int(file.readline().strip().split()[0])
            reax_off_diagonals = []
            for off in range(nr):
                temp_line = file.readline().strip().split()
                temp_list = [int(temp_line[0]), int(temp_line[1])]
                temp_list.extend([float(temp_line[num]) for num in (5, 6, 7)])
                reax_off_diagonals.append(ReaxOffDiagonals._make(temp_list))
        self.reax.append([reax_general, reax_atoms, reax_bonds, reax_off_diagonals])
        print("----------")


class Bond:
    def __init__(self, atom):
        self.atom = atom
        self.reax = {"order": None, "params": None, "overcoord": None}


class Atom:
    """Object containing atoms.

    :param coords: [tuple] - tuple of xyz coordinates:
        [float] - x coordinate
        [float] - y coordinate
        [float] - z coordinate

    As a default contains:
        coords [tuple]:
            0.0 [float] - x coordinate
            0.0 [float] - y coordinate
            0.0 [float] - z coordinate
        id (incrementing) [int] - atom id
        charge 0.0 [float] - atom charge
        velocity [tuple]:
            0.0 [float] - x velocity
            0.0 [float] - y velocity
            0.0 [float] - z velocity
        type 1 [int] - atom type
        name "Xxx" [string] - atom name
        mass 1.0 [float] - atom mass
        energy [dict]:
            "kinetic": 0.0 [float] - kinetic energy
            "potential": 0.0 [float] - potential energy
        neighbours [] [list] - list of close neighbours
        bin None [list] - bin the atom is currently in:
            x [int] - x bin
            y [int] - y bin
            z [int] - z bin
        bonds [list] - atoms that the atom is bonded to
    """

    _new_id = next(itertools.count())

    def __init__(self, coords):
        """
        :param coords: [tuple] - tuple of xyz coordinates:
            [float] - x coordinate
            [float] - y coordinate
            [float] - z coordinate
        """
        self.coords = coords
        self.id = Atom._new_id
        self.charge = 0.0
        self.velocity = (0.0, 0.0, 0.0)
        self.type = 1
        self.name = "Xxx"
        self.mass = 1.0
        self.energy = {"kinetic": 0.0, "potential": 0.0}
        self.neighbours = []
        self.bin = None
        self.bonds = None

    def __str__(self):
        return "Atom: ID " + str(self.id) + "; type " + str(self.type) + "; name " + self.name + \
               "; coords " + str(self.coords)


class Molecule:
    """Molecule consisting of Atom objects.

    :param atom: [[Atom]] - list of Atom objects

    As a default contains:
        mass [float] - mass of the molecule
        centre_of_mass [list] - coords of centre of mass
    """
    def __init__(self, atom):
        self.atoms = [atom]
        self.mass = 0.0
        self.centre_of_mass = np.array([0.0, 0.0, 0.0])
        self.name = ""
        self.energy = {"kinetic": None}
        self.velocity = None

    def __str__(self):
        return "Molecule: number of atoms " + str(len(self.atoms)) + "; total mass " + str(self.mass)


class AtomsSystem:  # TODO: change tuples to lists
    """System of Atom objects.

    As a default contains:
        atoms [] [list]:
            [Atom] - list of Atom objects
        number 0 [int] - number of atoms in the AtomsSystem
        timestep 0 [int] - timestep
        bounds [dict]:
            "xlo" [tuple]: - lower x boundary
                "f" [string] - boundary type
                -1.0 [float] - boundary value
            "xhi" [tuple]: - higher x boundary
                "f" [string] - boundary type
                1.0 [float] - boundary value
            "ylo" ...
            "yhi" ...
            "zlo" ...
            "zhi" ...
        max_type 1 [int] - maximum type of atom (often equal to number of types)
        time 0.0 [float] - time of simulation in fs
        bins [] [list] - bins for atoms during binning process
        molecules [] [[Molecule]] - molecules in the system
        molecules_ejected [] [[Molecule]] - molecules that are regarded as ejected
        reax_params None [namedtuple] - reaxFF parameters
    """

    def __init__(self):
        self._atoms = []
        self.number = 0
        self._next_iter_id = next(itertools.count())
        self.timestep = 0
        self.bounds = {
            "xlo": ("f", -1.0),
            "xhi": ("f", 1.0),
            "ylo": ("f", -1.0),
            "yhi": ("f", 1.0),
            "zlo": ("f", -1.0),
            "zhi": ("f", 1.0)
        }
        self.max_type = 1
        self.time = 0.0
        self.bins = None
        self.molecules = None
        self.molecules_ejected = None
        self.reax_params = None

    def __iter__(self):
        return self

    def next(self):
        iter_id = self._next_iter_id
        if iter_id > self.number:
            self._next_iter_id = next(itertools.count())
            raise StopIteration
        return self._atoms[iter_id]

    def __str__(self):
        return "AtomsSystem: number of atoms " + str(self.number)

    def addAtom(self, atom):
        """Adds Atom object instance to the AtomsSystem.

        :param atom: [Atom] - Atom object instance"""
        self._atoms.append(atom)
        self.number += 1

    @property
    def atoms(self):
        """List of Atom objects."""
        return self._atoms

    @atoms.setter
    def atoms(self, atoms):
        self._atoms = atoms
        self.number = len(self._atoms)

    def saveSystem(self, file_path, file_type=FileType.auto, control_dict=None):
        """Save file of the set type using data in the AtomsSystem.

        :param file_path: [string] - path to the output file.
        :param file_type: [string] - [Default: "auto"] type of the file. If "auto" than type determined by extension
                of the file specified in file_path:
            "lammps_data" - file read by LAMMPS's function read_data
            "lammpstrj" - LAMMPS trajectory file
        :param control_dict: [dict] - [Default: {}] dictionary consisting of additional controls that may be needed:
            "lammps_data_type": - type of data (one of) for lammps_data
                "charge"
            "properties": - list of properties (any of) for lammpstrj
                "id"
                "element"
                "type"
                "name"
                "x"
                "y"
                "z"
                "charge"
                "mass"
                "vx"
                "vy"
                "vz"
                "kinetic_energy"
                "potential_energy"
        """
        print("Starting 'saveSystem' procedure.")
        if control_dict is None:
            control_dict = {}
        if file_type == FileType.auto:
            print(" Detecting filetype automaticaly: ", end="")
            file_ext = os.path.splitext(file_path)[1][1:].strip().lower()
            if file_ext == "lammpstrj":
                file_type = FileType.lammpstrj
            elif file_ext == "xyz":
                file_type = FileType.xyz
            else:
                raise NotImplementedError(
                        "No implementation for the " + file_ext + " extension. Please set file_type explicitly."
                )
            print(file_type)
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        with open(file_path, "w") as file:
            if file_type == FileType.lammps_data:
                MGSaveFile.saveLammpsDataFile(self, file, control_dict)
            elif file_type == FileType.lammpstrj:
                MGSaveFile.saveLammpsFile(self, file, control_dict)
            else:
                raise NotImplementedError("No implementation for the " + str(file_type.name) + " file_type.")
        print("Successfully saved system to " + file_path + " as a " + str(file_type.name) + " file.")
        print("----------")

    def recalculateBounds(self):
        """Recalculate boundaries of the system."""
        print("Starting recalculateBounds procedure.")
        max_coords = [0, 0, 0]
        min_coords = [0, 0, 0]
        for atom in self.atoms:
            if atom.coords[0] > max_coords[0]:
                max_coords[0] = atom.coords[0]
            if atom.coords[0] < min_coords[0]:
                min_coords[0] = atom.coords[0]
            if atom.coords[1] > max_coords[1]:
                max_coords[1] = atom.coords[1]
            if atom.coords[1] < min_coords[1]:
                min_coords[1] = atom.coords[1]
            if atom.coords[2] > max_coords[2]:
                max_coords[2] = atom.coords[2]
            if atom.coords[2] < min_coords[2]:
                min_coords[2] = atom.coords[2]
        self.bounds["xhi"] = (self.bounds["xhi"][0], max_coords[0])
        self.bounds["xlo"] = (self.bounds["xlo"][0], min_coords[0])
        self.bounds["yhi"] = (self.bounds["yhi"][0], max_coords[1])
        self.bounds["ylo"] = (self.bounds["ylo"][0], min_coords[1])
        self.bounds["zhi"] = (self.bounds["zhi"][0], max_coords[2])
        self.bounds["zlo"] = (self.bounds["zlo"][0], min_coords[2])
        print("----------")

    def stretchBounds(self, x, y, z):
        """Add or substract number to/from boundaries.

        :param x: [list] - x boundaries:
            [float] - x min
            [float] - x max
        :param y: [list] - y boundaries:
            ...
        :param z: [list] - z boundaries:
            ...
        """
        print("Starting 'stretchBounds' procedure.")
        self.bounds["xlo"] = (self.bounds["xlo"][0], self.bounds["xlo"][1] + x[0])
        self.bounds["xhi"] = (self.bounds["xhi"][0], self.bounds["xhi"][1] + x[1])
        self.bounds["ylo"] = (self.bounds["ylo"][0], self.bounds["ylo"][1] + y[0])
        self.bounds["yhi"] = (self.bounds["yhi"][0], self.bounds["yhi"][1] + y[1])
        self.bounds["zlo"] = (self.bounds["zlo"][0], self.bounds["zlo"][1] + z[0])
        self.bounds["zhi"] = (self.bounds["zhi"][0], self.bounds["zhi"][1] + z[1])
        print("----------")

    def recalculateIDs(self):
        """Recalculate IDs of atoms in the system"""
        print("Starting 'recalculateIDs' procedure.")
        for i, atom in enumerate(self.atoms):
            atom.id = i + 1
        print("----------")

    def recalculateTypes(self):
        """Recalculate types of atoms based on their names."""
        print("Starting 'recalculateTypes' procedure.")
        atoms_names = []
        for atom in self.atoms:
            if atom.name not in atoms_names:
                atoms_names.append(atom.name)
            if atom.type == 1:
                atom.type = atoms_names.index(atom.name) + 1
            if atom.type > self.max_type:
                self.max_type = atom.type
        print("----------")

    def _doBinning(self, new_bounds_q=False, method="few_bins"):
        """Binning atoms.

        :param new_bounds_q: [bool] - should boundaries be recalculated?
        :param method: [string] - method for binning:
            "few_bins"
            "long_bins"
        """
        print("Starting '_doBinning' procedure.")
        self.recalculateIDs()
        if new_bounds_q:
            print("Recalculating boundaries.")
            self.recalculateBounds()
            self.stretchBounds([-0.1, 0.1], [-0.1, 0.1], [-0.1, 0.1])
        deltax = self.bounds["xhi"][1] - self.bounds["xlo"][1]
        deltay = self.bounds["yhi"][1] - self.bounds["ylo"][1]
        deltaz = self.bounds["zhi"][1] - self.bounds["zlo"][1]
        dx = deltax / 1000
        dy = deltay / 1000
        dz = deltaz / 1000
        if dx > 100:
            dx = 100
        if dx < 3:
            dx = 3
        if dy > 100:
            dy = 100
        if dy < 3:
            dy = 3
        if dz > 100:
            dz = 100
        if dz < 3:
            dz = 3

        if method == "long_bins":
            self.bins = [[set([]) for x in range(math.ceil(deltax / dx))],
                         [set([]) for y in range(math.ceil(deltay / dy))],
                         [set([]) for z in range(math.ceil(deltaz / dz))]]
            print("Created",
                  len(self.bins[0]), "x bins,",
                  len(self.bins[1]), "y bins and",
                  len(self.bins[2]), "z bins.")
            for counter_atom, atom in tqdm(enumerate(self.atoms), total=self.number, unit="atom"):
                binx = math.floor((atom.coords[0] - self.bounds["xlo"][1]) / dx)
                biny = math.floor((atom.coords[1] - self.bounds["ylo"][1]) / dy)
                binz = math.floor((atom.coords[2] - self.bounds["zlo"][1]) / dz)
                atom.bin = [binx, biny, binz]
                self.bins[0][binx].add(counter_atom + 1)
                self.bins[1][biny].add(counter_atom + 1)
                self.bins[2][binz].add(counter_atom + 1)
        elif method == "few_bins":
            from operator import attrgetter
            dx, dy, dz = 3, 3, 3
            for atom in self.atoms:
                binx = math.floor((atom.coords[0] - self.bounds["xlo"][1]) / dx)
                biny = math.floor((atom.coords[1] - self.bounds["ylo"][1]) / dy)
                binz = math.floor((atom.coords[2] - self.bounds["zlo"][1]) / dz)
                atom.bin = [binx, biny, binz]
            atoms_sorted = sorted(self.atoms, key=attrgetter("bin"))
            temp_bins = []
            last_bin = [-1, -1, -1]
            bin_counter = [0, 0, 0]
            for counter_atom, atom in tqdm(enumerate(atoms_sorted), total=self.number, unit="atom"):
                if last_bin == atom.bin:
                    temp_bins[-1][-1][-1].append(atom.id)
                elif last_bin[0] != atom.bin[0]:
                    temp_bins.append([[[atom.id]]])
                    bin_counter[0] += 1
                elif last_bin[1] != atom.bin[1]:
                    temp_bins[-1].append([[atom.id]])
                    bin_counter[1] += 1
                elif last_bin[2] != atom.bin[2]:
                    temp_bins[-1][-1].append([atom.id])
                    bin_counter[2] += 1
                else:
                    raise NotImplementedError
                last_bin = atom.bin
            self.bins = temp_bins
            bin_total = bin_counter[2]
            bin_counter[2] = round(bin_counter[2] / bin_counter[1], 2)
            bin_counter[1] = round(bin_counter[1] / bin_counter[0], 2)
            print("\nCreated total", bin_total, "bins (means:", bin_counter[0], "x bins,", bin_counter[1], "y bins and",
                  bin_counter[2], "z).")
        else:
            print("Using horrible but easiest binning method.")
            self.bins = [  # TODO: Better list initialization (possible?)
                [
                    [np.array([], dtype=np.int_)
                        for z in range(math.ceil(deltaz / dz))]
                    for y in range(math.ceil(deltay / dy))]
                for x in range(math.ceil(deltax / dx))]

            print("Created " + str(math.ceil(deltaz / dz)) + " x " + str(math.ceil(deltay / dy)) + " x " +
                  str(math.ceil(deltay / dy)) + " bins.")
            for counter_atom, atom in tqdm(enumerate(self.atoms), total=self.number, unit="atom"):
                binx = math.floor((atom.coords[0] - self.bounds["xlo"][1]) / dx)
                biny = math.floor((atom.coords[1] - self.bounds["ylo"][1]) / dy)
                binz = math.floor((atom.coords[2] - self.bounds["zlo"][1]) / dz)
                atom.bin = [binx, biny, binz]
                self.bins[binx][biny][binz] = np.append(self.bins[binx][biny][binz], np.array([counter_atom + 1]))
        print("\n----------")

    def doCloseNeighbours(self, method="few_bins"):
        """Create list of closest neighbours for each atom in the system.

        :param method: [string] - method of binning:
            "few_bins"
            "long_bins"
        """
        if self.bins is None:
            print("Binning procedure execution not found. Starting '_doBinning' procedure.")
            self._doBinning(method=method)
        print("Starting 'doCloseNeighbours' procedure.")

        if method == "long_bins":
            for atom_counter, atom in tqdm(enumerate(self.atoms), total=self.number, unit="atom"):
                x = atom.bin[0]
                y = atom.bin[1]
                z = atom.bin[2]
                big_bin = set()
                for i in [x-1, x, x+1]:
                    for j in [y-1, y, y+1]:
                        for k in [z-1, z, z+1]:
                            try:
                                big_bin = big_bin.union(self.bins[0][i].intersection(self.bins[1][j], self.bins[2][k]))
                            except IndexError:
                                continue
                atom.neighbours = np.array([near for near in big_bin if near != atom_counter + 1])
        elif method == "few_bins":
            pbar = tqdm(total=self.number, unit="atom")
            for counter_x, bin_x in enumerate(self.bins):
                for bin_y in bin_x:
                    for bin_z in bin_y:
                        big_bin = []
                        bin_now = self.atoms[bin_z[0] - 1].bin
                        for x in [-1, 0, 1]:
                            try:
                                bin_check = self.atoms[self.bins[counter_x + x][0][0][0] - 1].bin
                            except IndexError:
                                continue
                            if bin_check[0] > bin_now[0] + 1 or counter_x + x < 0:
                                continue
                            counter_y = 0
                            while True:
                                try:
                                    bin_check = self.atoms[self.bins[counter_x + x][counter_y][0][0] - 1].bin
                                except IndexError:
                                    break
                                if bin_check[1] > bin_now[1] + 1:
                                    break
                                for y in [-1, 0, 1]:
                                    if bin_check[0:2] == [bin_now[0] + x, bin_now[1] + y]:
                                        counter_z = 0
                                        while True:
                                            try:
                                                bin_check = self.atoms[self.bins[counter_x + x][counter_y][counter_z]
                                                                       [0] - 1].bin
                                            except IndexError:
                                                break
                                            if bin_check[2] > bin_now[2] + 1:
                                                break
                                            for z in [-1, 0, 1]:
                                                if bin_check == [bin_now[0] + x, bin_now[1] + y, bin_now[2] + z]:
                                                    big_bin.extend(self.bins[counter_x + x][counter_y][counter_z])
                                            counter_z += 1
                                counter_y += 1
                        for atom_nr in bin_z:
                            self.atoms[atom_nr - 1].neighbours = np.array([near for near in big_bin if near !=
                                                                           atom_nr], dtype=np.int_)
                            pbar.update()
            pbar.close()
        else:
            print("Using old, horrible but easy method.")
            xmax = len(self.bins) - 1
            ymax = len(self.bins[0]) - 1
            zmax = len(self.bins[0][0]) - 1
            for atom_counter, atom in tqdm(enumerate(self.atoms), total=self.number, unit="atom"):
                x = atom.bin[0]
                y = atom.bin[1]
                z = atom.bin[2]
                big_bin = [neigh for neigh in self.bins[x][y][z] if neigh != atom_counter + 1]
                if 0 < x < xmax:
                    big_bin.extend(self.bins[x-1][y][z])
                    big_bin.extend(self.bins[x+1][y][z])
                    if 0 < y < ymax:
                        big_bin.extend(self.bins[x][y-1][z])
                        big_bin.extend(self.bins[x][y+1][z])
                        big_bin.extend(self.bins[x-1][y-1][z])
                        big_bin.extend(self.bins[x+1][y-1][z])
                        big_bin.extend(self.bins[x-1][y+1][z])
                        big_bin.extend(self.bins[x+1][y+1][z])
                        if 0 < z < zmax:
                            big_bin.extend(self.bins[x][y][z-1])
                            big_bin.extend(self.bins[x][y][z+1])
                            big_bin.extend(self.bins[x-1][y][z-1])
                            big_bin.extend(self.bins[x+1][y][z-1])
                            big_bin.extend(self.bins[x-1][y][z+1])
                            big_bin.extend(self.bins[x+1][y][z+1])
                            big_bin.extend(self.bins[x][y-1][z-1])
                            big_bin.extend(self.bins[x][y+1][z-1])
                            big_bin.extend(self.bins[x-1][y-1][z-1])
                            big_bin.extend(self.bins[x+1][y-1][z-1])
                            big_bin.extend(self.bins[x-1][y+1][z-1])
                            big_bin.extend(self.bins[x+1][y+1][z-1])
                            big_bin.extend(self.bins[x][y-1][z+1])
                            big_bin.extend(self.bins[x][y+1][z+1])
                            big_bin.extend(self.bins[x-1][y-1][z+1])
                            big_bin.extend(self.bins[x+1][y-1][z+1])
                            big_bin.extend(self.bins[x-1][y+1][z+1])
                            big_bin.extend(self.bins[x+1][y+1][z+1])
                        else:
                            if z > 0:
                                big_bin.extend(self.bins[x][y][z-1])
                                big_bin.extend(self.bins[x-1][y][z-1])
                                big_bin.extend(self.bins[x+1][y][z-1])
                                big_bin.extend(self.bins[x][y-1][z-1])
                                big_bin.extend(self.bins[x][y+1][z-1])
                                big_bin.extend(self.bins[x-1][y-1][z-1])
                                big_bin.extend(self.bins[x+1][y-1][z-1])
                                big_bin.extend(self.bins[x-1][y+1][z-1])
                                big_bin.extend(self.bins[x+1][y+1][z-1])
                            if z < zmax:
                                big_bin.extend(self.bins[x][y][z+1])
                                big_bin.extend(self.bins[x-1][y][z+1])
                                big_bin.extend(self.bins[x+1][y][z+1])
                                big_bin.extend(self.bins[x][y-1][z+1])
                                big_bin.extend(self.bins[x][y+1][z+1])
                                big_bin.extend(self.bins[x-1][y-1][z+1])
                                big_bin.extend(self.bins[x+1][y-1][z+1])
                                big_bin.extend(self.bins[x-1][y+1][z+1])
                                big_bin.extend(self.bins[x+1][y+1][z+1])
                    else:
                        if y > 0:
                            big_bin.extend(self.bins[x][y-1][z])
                            big_bin.extend(self.bins[x-1][y-1][z])
                            big_bin.extend(self.bins[x+1][y-1][z])
                            if 0 < z < zmax:
                                big_bin.extend(self.bins[x][y][z-1])
                                big_bin.extend(self.bins[x][y][z+1])
                                big_bin.extend(self.bins[x-1][y][z-1])
                                big_bin.extend(self.bins[x+1][y][z-1])
                                big_bin.extend(self.bins[x-1][y][z+1])
                                big_bin.extend(self.bins[x+1][y][z+1])
                                big_bin.extend(self.bins[x][y-1][z-1])
                                big_bin.extend(self.bins[x-1][y-1][z-1])
                                big_bin.extend(self.bins[x+1][y-1][z-1])
                                big_bin.extend(self.bins[x][y-1][z+1])
                                big_bin.extend(self.bins[x-1][y-1][z+1])
                                big_bin.extend(self.bins[x+1][y-1][z+1])
                            else:
                                if z > 0:
                                    big_bin.extend(self.bins[x][y][z-1])
                                    big_bin.extend(self.bins[x-1][y][z-1])
                                    big_bin.extend(self.bins[x+1][y][z-1])
                                    big_bin.extend(self.bins[x][y-1][z-1])
                                    big_bin.extend(self.bins[x-1][y-1][z-1])
                                    big_bin.extend(self.bins[x+1][y-1][z-1])
                                if z < zmax:
                                    big_bin.extend(self.bins[x][y][z+1])
                                    big_bin.extend(self.bins[x-1][y][z+1])
                                    big_bin.extend(self.bins[x+1][y][z+1])
                                    big_bin.extend(self.bins[x][y-1][z+1])
                                    big_bin.extend(self.bins[x-1][y-1][z+1])
                                    big_bin.extend(self.bins[x+1][y-1][z+1])
                        if y < ymax:
                            big_bin.extend(self.bins[x][y+1][z])
                            big_bin.extend(self.bins[x-1][y+1][z])
                            big_bin.extend(self.bins[x+1][y+1][z])
                            if 0 < z < zmax:
                                big_bin.extend(self.bins[x][y][z-1])
                                big_bin.extend(self.bins[x][y][z+1])
                                big_bin.extend(self.bins[x-1][y][z-1])
                                big_bin.extend(self.bins[x+1][y][z-1])
                                big_bin.extend(self.bins[x-1][y][z+1])
                                big_bin.extend(self.bins[x+1][y][z+1])
                                big_bin.extend(self.bins[x][y+1][z-1])
                                big_bin.extend(self.bins[x-1][y+1][z-1])
                                big_bin.extend(self.bins[x+1][y+1][z-1])
                                big_bin.extend(self.bins[x][y+1][z+1])
                                big_bin.extend(self.bins[x-1][y+1][z+1])
                                big_bin.extend(self.bins[x+1][y+1][z+1])
                            else:
                                if z > 0:
                                    big_bin.extend(self.bins[x][y][z-1])
                                    big_bin.extend(self.bins[x-1][y][z-1])
                                    big_bin.extend(self.bins[x+1][y][z-1])
                                    big_bin.extend(self.bins[x][y+1][z-1])
                                    big_bin.extend(self.bins[x-1][y+1][z-1])
                                    big_bin.extend(self.bins[x+1][y+1][z-1])
                                if z < zmax:
                                    big_bin.extend(self.bins[x][y][z+1])
                                    big_bin.extend(self.bins[x-1][y][z+1])
                                    big_bin.extend(self.bins[x+1][y][z+1])
                                    big_bin.extend(self.bins[x][y+1][z+1])
                                    big_bin.extend(self.bins[x-1][y+1][z+1])
                                    big_bin.extend(self.bins[x+1][y+1][z+1])
                else:
                    if x > 0:
                        big_bin.extend(self.bins[x-1][y][z])
                        if 0 < y < ymax:
                            big_bin.extend(self.bins[x][y-1][z])
                            big_bin.extend(self.bins[x][y+1][z])
                            big_bin.extend(self.bins[x-1][y-1][z])
                            big_bin.extend(self.bins[x-1][y+1][z])
                            if 0 < z < zmax:
                                big_bin.extend(self.bins[x][y][z-1])
                                big_bin.extend(self.bins[x][y][z+1])
                                big_bin.extend(self.bins[x-1][y][z-1])
                                big_bin.extend(self.bins[x-1][y][z+1])
                                big_bin.extend(self.bins[x][y-1][z-1])
                                big_bin.extend(self.bins[x][y+1][z-1])
                                big_bin.extend(self.bins[x-1][y-1][z-1])
                                big_bin.extend(self.bins[x-1][y+1][z-1])
                                big_bin.extend(self.bins[x][y-1][z+1])
                                big_bin.extend(self.bins[x][y+1][z+1])
                                big_bin.extend(self.bins[x-1][y-1][z+1])
                                big_bin.extend(self.bins[x-1][y+1][z+1])
                            else:
                                if z > 0:
                                    big_bin.extend(self.bins[x][y][z-1])
                                    big_bin.extend(self.bins[x-1][y][z-1])
                                    big_bin.extend(self.bins[x][y-1][z-1])
                                    big_bin.extend(self.bins[x][y+1][z-1])
                                    big_bin.extend(self.bins[x-1][y-1][z-1])
                                    big_bin.extend(self.bins[x-1][y+1][z-1])
                                if z < zmax:
                                    big_bin.extend(self.bins[x][y][z+1])
                                    big_bin.extend(self.bins[x-1][y][z+1])
                                    big_bin.extend(self.bins[x][y-1][z+1])
                                    big_bin.extend(self.bins[x][y+1][z+1])
                                    big_bin.extend(self.bins[x-1][y-1][z+1])
                                    big_bin.extend(self.bins[x-1][y+1][z+1])
                        else:
                            if y > 0:
                                big_bin.extend(self.bins[x][y-1][z])
                                big_bin.extend(self.bins[x-1][y-1][z])
                                if 0 < z < zmax:
                                    big_bin.extend(self.bins[x][y][z-1])
                                    big_bin.extend(self.bins[x][y][z+1])
                                    big_bin.extend(self.bins[x-1][y][z-1])
                                    big_bin.extend(self.bins[x-1][y][z+1])
                                    big_bin.extend(self.bins[x][y-1][z-1])
                                    big_bin.extend(self.bins[x-1][y-1][z-1])
                                    big_bin.extend(self.bins[x][y-1][z+1])
                                    big_bin.extend(self.bins[x-1][y-1][z+1])
                                else:
                                    if z > 0:
                                        big_bin.extend(self.bins[x][y][z-1])
                                        big_bin.extend(self.bins[x-1][y][z-1])
                                        big_bin.extend(self.bins[x][y-1][z-1])
                                        big_bin.extend(self.bins[x-1][y-1][z-1])
                                    if z < zmax:
                                        big_bin.extend(self.bins[x][y][z+1])
                                        big_bin.extend(self.bins[x-1][y][z+1])
                                        big_bin.extend(self.bins[x][y-1][z+1])
                                        big_bin.extend(self.bins[x-1][y-1][z+1])
                            if y < ymax:
                                big_bin.extend(self.bins[x][y+1][z])
                                big_bin.extend(self.bins[x-1][y+1][z])
                                if 0 < z < zmax:
                                    big_bin.extend(self.bins[x][y][z-1])
                                    big_bin.extend(self.bins[x][y][z+1])
                                    big_bin.extend(self.bins[x-1][y][z-1])
                                    big_bin.extend(self.bins[x-1][y][z+1])
                                    big_bin.extend(self.bins[x][y+1][z-1])
                                    big_bin.extend(self.bins[x-1][y+1][z-1])
                                    big_bin.extend(self.bins[x][y+1][z+1])
                                    big_bin.extend(self.bins[x-1][y+1][z+1])
                                else:
                                    if z > 0:
                                        big_bin.extend(self.bins[x][y][z-1])
                                        big_bin.extend(self.bins[x-1][y][z-1])
                                        big_bin.extend(self.bins[x][y+1][z-1])
                                        big_bin.extend(self.bins[x-1][y+1][z-1])
                                    if z < zmax:
                                        big_bin.extend(self.bins[x][y][z+1])
                                        big_bin.extend(self.bins[x-1][y][z+1])
                                        big_bin.extend(self.bins[x][y+1][z+1])
                                        big_bin.extend(self.bins[x-1][y+1][z+1])
                    if x < xmax:
                        big_bin.extend(self.bins[x+1][y][z])
                        if 0 < y < ymax:
                            big_bin.extend(self.bins[x][y-1][z])
                            big_bin.extend(self.bins[x][y+1][z])
                            big_bin.extend(self.bins[x+1][y-1][z])
                            big_bin.extend(self.bins[x+1][y+1][z])
                            if 0 < z < zmax:
                                big_bin.extend(self.bins[x][y][z-1])
                                big_bin.extend(self.bins[x][y][z+1])
                                big_bin.extend(self.bins[x+1][y][z-1])
                                big_bin.extend(self.bins[x+1][y][z+1])
                                big_bin.extend(self.bins[x][y-1][z-1])
                                big_bin.extend(self.bins[x][y+1][z-1])
                                big_bin.extend(self.bins[x+1][y-1][z-1])
                                big_bin.extend(self.bins[x+1][y+1][z-1])
                                big_bin.extend(self.bins[x][y-1][z+1])
                                big_bin.extend(self.bins[x][y+1][z+1])
                                big_bin.extend(self.bins[x+1][y-1][z+1])
                                big_bin.extend(self.bins[x+1][y+1][z+1])
                            else:
                                if z > 0:
                                    big_bin.extend(self.bins[x][y][z-1])
                                    big_bin.extend(self.bins[x+1][y][z-1])
                                    big_bin.extend(self.bins[x][y-1][z-1])
                                    big_bin.extend(self.bins[x][y+1][z-1])
                                    big_bin.extend(self.bins[x+1][y-1][z-1])
                                    big_bin.extend(self.bins[x+1][y+1][z-1])
                                if z < zmax:
                                    big_bin.extend(self.bins[x][y][z+1])
                                    big_bin.extend(self.bins[x+1][y][z+1])
                                    big_bin.extend(self.bins[x][y-1][z+1])
                                    big_bin.extend(self.bins[x][y+1][z+1])
                                    big_bin.extend(self.bins[x+1][y-1][z+1])
                                    big_bin.extend(self.bins[x+1][y+1][z+1])
                        else:
                            if y > 0:
                                big_bin.extend(self.bins[x][y-1][z])
                                big_bin.extend(self.bins[x+1][y-1][z])
                                if 0 < z < zmax:
                                    big_bin.extend(self.bins[x][y][z-1])
                                    big_bin.extend(self.bins[x][y][z+1])
                                    big_bin.extend(self.bins[x+1][y][z-1])
                                    big_bin.extend(self.bins[x+1][y][z+1])
                                    big_bin.extend(self.bins[x][y-1][z-1])
                                    big_bin.extend(self.bins[x+1][y-1][z-1])
                                    big_bin.extend(self.bins[x][y-1][z+1])
                                    big_bin.extend(self.bins[x+1][y-1][z+1])
                                else:
                                    if z > 0:
                                        big_bin.extend(self.bins[x][y][z-1])
                                        big_bin.extend(self.bins[x+1][y][z-1])
                                        big_bin.extend(self.bins[x][y-1][z-1])
                                        big_bin.extend(self.bins[x+1][y-1][z-1])
                                    if z < zmax:
                                        big_bin.extend(self.bins[x][y][z+1])
                                        big_bin.extend(self.bins[x+1][y][z+1])
                                        big_bin.extend(self.bins[x][y-1][z+1])
                                        big_bin.extend(self.bins[x+1][y-1][z+1])
                            if y < ymax:
                                big_bin.extend(self.bins[x][y+1][z])
                                big_bin.extend(self.bins[x+1][y+1][z])
                                if 0 < z < zmax:
                                    big_bin.extend(self.bins[x][y][z-1])
                                    big_bin.extend(self.bins[x][y][z+1])
                                    big_bin.extend(self.bins[x+1][y][z-1])
                                    big_bin.extend(self.bins[x+1][y][z+1])
                                    big_bin.extend(self.bins[x][y+1][z-1])
                                    big_bin.extend(self.bins[x+1][y+1][z-1])
                                    big_bin.extend(self.bins[x][y+1][z+1])
                                    big_bin.extend(self.bins[x+1][y+1][z+1])
                                else:
                                    if z > 0:
                                        big_bin.extend(self.bins[x][y][z-1])
                                        big_bin.extend(self.bins[x+1][y][z-1])
                                        big_bin.extend(self.bins[x][y+1][z-1])
                                        big_bin.extend(self.bins[x+1][y+1][z-1])
                                    if z < zmax:
                                        big_bin.extend(self.bins[x][y][z+1])
                                        big_bin.extend(self.bins[x+1][y][z+1])
                                        big_bin.extend(self.bins[x][y+1][z+1])
                                        big_bin.extend(self.bins[x+1][y+1][z+1])
                atom.neighbours = np.array(big_bin, dtype=np.int_)
        print("\n----------")

    def findBonds(self, method=None):  # TODO: Add additional methods of calculating bonds
        """Find bonds between atoms in the system.

        :param method: [string] - method used
            None
            "reaxFF"
        """
        print("Starting 'findBonds' procedure.")
        if method is None:
            print("Using default method: atoms closer to each other than 2 A are regarded as bonded.")
            for atom in tqdm(self.atoms, total=self.number, unit="atom"):
                atom.bonds = []
                for neighbour in atom.neighbours:
                    distance = sum([(at_coord - neigh_coord)**2 for at_coord, neigh_coord in
                                    zip(atom.coords, self.atoms[neighbour - 1].coords)])
                    if distance < 2**2:
                        atom.bonds.append(Bond(neighbour))
        elif method is "reaxFF":
            print("Using reax method.")
            Reax = namedtuple("Reax", ["p_boc1", "p_boc2", "p_boc3", "p_boc4", "p_boc5", "p_bo1", "p_bo2", "p_bo3",
                                       "p_bo4", "p_bo5", "p_bo6", "Val_at", "Val_boc_at", "Val_neigh", "Val_boc_neigh",
                                       "ro_sigma", "ro_pi", "ro_pipi"])
            r_cutoff = 3

            def BO_sigma_prim(reax, r_ij):
                if reax.ro_sigma <= 0:
                    return 0
                else:
                    return math.exp(reax.p_bo1 * (r_ij / reax.ro_sigma)**reax.p_bo2)

            def BO_pi_prim(reax, r_ij):
                if reax.ro_pi <= 0:
                    return 0
                else:
                    return math.exp(reax.p_bo3 * (r_ij / reax.ro_pi)**reax.p_bo4)

            def BO_pipi_prim(reax, r_ij):
                if reax.ro_pipi <= 0:
                    return 0
                else:
                    return math.exp(reax.p_bo5 * (r_ij / reax.ro_pipi)**reax.p_bo6)

            def BO_sigma(reax):
                return reax["order"][1][0] * f1(reax) * f4(reax) * f5(reax)

            def BO_pi(reax):
                return reax["order"][1][1] * f1(reax)**2 * f4(reax) * f5(reax)

            def BO_pipi(reax):
                return reax["order"][1][2] * f1(reax)**2 * f4(reax) * f5(reax)

            def f1(reax):
                temp_f2 = f2(reax)
                temp_f23 = temp_f2 + f3(reax)
                return 1/2*(((reax["params"].Val_at + temp_f2) / (reax["params"].Val_at + temp_f23)) +
                            ((reax["params"].Val_neigh + temp_f2) / (reax["params"].Val_neigh + temp_f23)))

            def f2(reax):
                return math.exp(-reax["params"].p_boc1 * reax["overcoord"][0][0]) + \
                       math.exp(-reax["params"].p_boc1 * reax["overcoord"][1][0])

            def f3(reax):
                return -1/reax["params"].p_boc2 * \
                       math.log(1/2 * (math.exp(-reax["params"].p_boc2 * reax["overcoord"][0][0]) +
                                       math.exp(-reax["params"].p_boc2 * reax["overcoord"][1][0])))

            def f4(reax):
                return 1 / (1 + math.exp(-reax["params"].p_boc3 *
                                         (reax["params"].p_boc4 * reax["order"][0]**2 -
                                          reax["overcoord"][0][1]) + reax["params"].p_boc5))

            def f5(reax):
                return 1 / (1 + math.exp(-reax["params"].p_boc3 *
                                         (reax["params"].p_boc4 * reax["order"][0]**2 -
                                          reax["overcoord"][1][1]) + reax["params"].p_boc5))

            for atom in tqdm(self.atoms, total=self.number, unit="atoms"):
                atom.bonds = []
                for neighbour in atom.neighbours:
                    neigh_data = [self.atoms[neighbour - 1].name, self.atoms[neighbour - 1].coords]
                    r = sum([(at_coord - neigh_coord)**2 for at_coord, neigh_coord in zip(atom.coords, neigh_data[1])])
                    if r > r_cutoff**2:
                        continue
                    atom.bonds.append(Bond(neighbour))

                    r = math.sqrt(r)
                    names = [param[0] for param in self.reax_params[1]]
                    idxs = [names.index(atom.name), names.index(neigh_data[0])]
                    atoms = [self.reax_params[1][idxs[i]] for i in (0, 1)]
                    bonds_idxs = [[param[0] - 1, param[1] - 1] for param in self.reax_params[2]]
                    bond = [self.reax_params[2][i] for i, idx in enumerate(bonds_idxs) if idx == idxs or idx ==
                            [idxs[1], idxs[0]]][0]
                    if atom.name != neigh_data[0]:
                        diagonal_idxs = [[param[0] - 1, param[1] - 1] for param in self.reax_params[3]]
                        diagonal = [self.reax_params[3][i] for i, idx in enumerate(diagonal_idxs) if idx == idxs or
                                    idx == [idxs[1], idxs[0]]][0]
                        diagonal = [diagonal.ro_sigma, diagonal.ro_pi, diagonal.ro_pipi]
                    else:
                        diagonal = [atoms[0].ro_sigma, atoms[0].ro_pi, atoms[0].ro_pipi]
                    atom.bonds[-1].reax["params"] = Reax(self.reax_params[0].p_boc1, self.reax_params[0].p_boc2,
                                                         math.sqrt(atoms[0].p_boc3 * atoms[1].p_boc3),
                                                         math.sqrt(atoms[0].p_boc4 * atoms[1].p_boc4),
                                                         math.sqrt(atoms[0].p_boc5 * atoms[1].p_boc5),
                                                         bond.p_bo1, bond.p_bo2, bond.p_bo3, bond.p_bo4, bond.p_bo5,
                                                         bond.p_bo6, atoms[0].Val, atoms[0].Val_boc, atoms[1].Val,
                                                         atoms[1].Val_boc, diagonal[0], diagonal[1], diagonal[2])
                    atom.bonds[-1].reax["order"] = [None, [BO_sigma_prim(atom.bonds[-1].reax["params"], r),
                                                           BO_pi_prim(atom.bonds[-1].reax["params"], r),
                                                           BO_pipi_prim(atom.bonds[-1].reax["params"], r)]]
                    atom.bonds[-1].reax["order"][0] = sum(atom.bonds[-1].reax["order"][1])
            for atom in tqdm(self.atoms, unit="atoms"):
                BO_sum = sum([bond.reax["order"][0] for bond in atom.bonds])
                for bond in atom.bonds:
                    Delta = -bond.reax["params"].Val_at + BO_sum
                    Delta_boc = -bond.reax["params"].Val_boc_at + BO_sum
                    bond.reax["overcoord"] = [[Delta, Delta_boc]]
            for atom in tqdm(self.atoms, unit="atoms"):
                for bond in atom.bonds:
                    bond.reax["overcoord"].append(self.atoms[bond.atom - 1].bonds[0].reax["overcoord"][0])
            for atom in tqdm(self.atoms, total=self.number, unit="atoms"):
                new_bonds = []
                for bond in atom.bonds:
                    bond.reax["order"][1] = [BO_sigma(bond.reax), BO_pi(bond.reax), BO_pipi(bond.reax)]
                    bond.reax["order"][0] = sum(bond.reax["order"][1])
                    if bond.reax["order"][0] > 0.3:
                        new_bonds.append(bond)
                atom.bonds = new_bonds
        print("\n----------")

    def findMolecules(self):
        """Find molecules in the system."""
        print("Starting 'findMolecules' procedure.")
        atoms_in_molecules = set()
        self.molecules = []
        for atom_counter, atom in tqdm(enumerate(self.atoms), total=self.number, unit="atom"):
            if (atom_counter + 1) not in atoms_in_molecules:
                self.molecules.append(Molecule(atom))
                atoms_in_molecules.add(atom_counter + 1)
                for mol_atom in self.molecules[-1].atoms:
                    for bond in mol_atom.bonds:
                        if bond.atom not in atoms_in_molecules:
                            self.molecules[-1].atoms.append(self.atoms[bond.atom - 1])
                            atoms_in_molecules.add(bond.atom)
        for molecule in self.molecules:
            elements = []
            for atom in molecule.atoms:
                molecule.mass += atom.mass
                elements.append(atom.name)
            elements = Counter(elements)
            name_list = []
            for letter in ("C", "N", "O"):
                if letter in elements.keys():
                    if elements[letter] != 1:
                        name_list.append(letter + str(elements[letter]))
                    else:
                        name_list.append(letter)
            if "H" in elements.keys():
                if elements["H"] != 1:
                    name_list.insert(1, "H" + str(elements["H"]))
                else:
                    name_list.insert(1, "H")
            molecule.name = "".join(name_list)
        print("\nFound " + str(len(self.molecules)) + " molecules.")
        print("----------")

    def findMoleculesCOM(self):
        """Find molecules' centres of mass."""
        if self.molecules is None:
            print("No molecules' lookup found. Executing 'findMolecules' procedure.")
            self.findMolecules()
        print("Starting 'findMoleculesCOM' procedure.")
        for molecule in tqdm(self.molecules, total=len(self.molecules), unit="mol"):
            mass_sum = 0
            coords_sum = np.array([0.0, 0.0, 0.0])
            for atom in molecule.atoms:
                mass_sum += atom.mass
                atom_coords = np.array(atom.coords)
                coords_sum += atom_coords * atom.mass
            molecule.centre_of_mass = coords_sum / mass_sum
        print("\n----------")

    def findMoleculesKineticEnergy(self):
        """Find molecules' kinetic energies."""
        if self.molecules is None:
            print("No molecules' lookup found. Executing 'findMolecules' procedure.")
            self.findMolecules()
            self.findMoleculesCOM()
        print("Starting 'findMoleculesKineticEnergy' procedure.")
        stop = False
        while True:
            for molecule in tqdm(self.molecules, unit="mol"):
                if molecule.centre_of_mass is None:
                    break
                v = np.array([0.0, 0.0, 0.0])
                for atom in molecule.atoms:
                    v += np.array(atom.velocity) * atom.mass
                v = v / molecule.mass
                molecule.velocity = [math.sqrt(sum(v**2)), v]
                molecule.energy["kinetic"] = [molecule.mass * molecule.velocity[0]**2 / 2, molecule.mass * v**2 / 2]
            else:
                stop = True

            if stop:
                break
            else:
                print("Molecules lack centre of mass. Executing 'findMoleculesCOM' procedure.")
                self.findMoleculesCOM()

    def findEjected(self, height):
        """Find ejected molecules.

        Find molecules that are higher than specified height (in z direction).
        :param height: [float] - height above which molecule is regarded as ejected
        """
        if self.molecules is None:
            print("No molecules' lookup found. Executing 'findMolecules' procedure.")
            self.findMolecules()
        print("Starting 'findEjected' procedure.")
        self.molecules_ejected = []
        for molecule in tqdm(self.molecules, total=len(self.molecules), unit="mol"):
            try:
                for atom in molecule.atoms:
                    if atom.coords[2] < height:
                        raise StopError
                self.molecules_ejected.append(molecule)
            except StopError:
                continue
        print("\nFound " + str(len(self.molecules_ejected)) + " molecules that are ejected.")
        print("----------")

    def saveMassSpectrum(self, file_path, rounding=0):
        """Save mass spectrum of the system to the file.

        Save mass spectrum of the ejected molecules to the file.
        :param file_path: [string] - path to the mass spectrum file
        :param rounding: [int] - number of digits after decimal point for masses
        """
        if self.molecules_ejected is None:
            print("No lookup of ejected molecules found in the system. Starting 'findEjected' procedure with height 0."
                  "If this is not what you want run 'findEjected' procedure before this one.")
            self.findEjected(0.0)
        print("Starting 'saveMassSpectrum' procedure.")
        if len(self.molecules_ejected) > 0:
            spectrum_numbers = []
            spectrum_masses = []
            spectrum_names = []
            for molecule in tqdm(self.molecules_ejected, total=len(self.molecules_ejected), unit="mol"):
                mass = round(molecule.mass, rounding)
                if mass not in spectrum_masses:
                    spectrum_masses.append(mass)
                    spectrum_names.append([molecule.name])
                    spectrum_numbers.append(1)
                else:
                    spectrum_numbers[spectrum_masses.index(mass)] += 1
                    if molecule.name not in spectrum_names[spectrum_masses.index(mass)]:
                        spectrum_names[spectrum_masses.index(mass)].append(molecule.name)
            spectrum = list(zip(spectrum_masses, spectrum_numbers, spectrum_names))
            spectrum.sort()
            print("\nFound", len(spectrum), "distinct masses, from", spectrum[0][0], "to", str(spectrum[-1][0]) +
                  ", consisting of", sum(map(len, spectrum_names)), "distinct compounds.")
            print("Saving data to the file.")
            os.makedirs(os.path.dirname(file_path), exist_ok=True)
            with open(file_path, "w") as f:
                for item in tqdm(spectrum):
                    f.write("\t".join([str(i) for i in item[:2]]) + "\t" + ", ".join([str(i) for i in item[2]]) + "\n")
            print("Saved data to file " + file_path)
        else:
            print("No ejected molecules in the system. Nothing to prepare mass spectrum from. Nothing was saved.")
        print("----------")

    def saveMasses(self, file_path):
        """Save different kinds of masses.

        :param file_path: [string] - path to the file
        """
        print("Starting 'saveMasses' procedure.")
        C60 = 0
        graph = 0
        phe = 0
        for molecule in tqdm(self.molecules_ejected, total=len(self.molecules_ejected), unit="mol"):
            for atom in molecule.atoms:
                if atom.type == 6:
                    C60 += atom.mass
                elif atom.type == 1:
                    graph += atom.mass
                else:
                    phe += atom.mass
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        with open(file_path, "w") as file:
            file.write("Total mass:\t" + str(C60 + graph + phe) + "\n")
            file.write("C60 mass:\t" + str(C60) + "\n")
            file.write("Graphene mass:\t" + str(graph) + "\n")
            file.write("Phenylalanine mass:\t" + str(phe) + "\n")
        print("----------")

    def saveEnergies(self, file_path):
        """Save different kinds of energies.

        :param file_path: [str] - path to the file
        """
        print("Starting 'saveEnergies' procedure.")
        COM = []
        internal = []
        for molecule in tqdm(self.molecules_ejected, unit="mol"):
            if molecule.name == "C9H11NO2":
                COM.append(molecule.energy["kinetic"][0])
                energy_total = 0
                for atom in molecule.atoms:
                    v = np.array(atom.velocity) - molecule.velocity[1]
                    energy_total += atom.mass * math.sqrt(sum(v**2))**2 / 2
                internal.append(energy_total)
        COM.sort()
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        with open(file_path, "w") as file:
            file.write("Phenylalanine COM kinetic energies:\t" + "\t".join([str(energy) for energy in COM]) + "\n")
            file.write("Phenylalanine internal kinetic energies:\t" + "\t".join([str(energy) for energy in internal]) +
                       "\n")
        print("----------")

    def saveSpotPlot(self, file_path, r=10000000):  # TODO: Write better docstring
        """Save data for spot plot.

        :param file_path: [string] - path to the file
        :param r: [float] - radius of the hemisphere
        """
        print("Starting 'saveSpotPlot' procedure.")
        coords_all = []
        masses_all = []
        for molecule in tqdm(self.molecules_ejected, unit="mol"):
            masses_all.append(molecule.mass)
            v = np.array(molecule.velocity[1])
            coords = np.array(molecule.centre_of_mass)
            a = sum(v**2)
            b = sum(2 * v * coords)
            c = -r**2+sum(coords**2)
            t = (-b + math.sqrt(b**2 - 4 * a * c)) / (2 * a)
            if t < 0:
                t2 = (-b - math.sqrt(b**2 - 4 * a * c)) / (2 * a)
                if t2 < 0:
                    if t2 > t:
                        t = t2
                else:
                    t = t2
            coords_new = coords + v * t
            theta = math.atan2(coords_new[1], coords_new[0])
            phi = math.acos(coords_new[2] / r)
            coords_all.append([phi, theta])
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        with open(file_path, "w") as file:
            file.write("phi [rad]\ttheta [rad]\tmass [g/mol]\n")
            for coord, mass in zip(coords_all, masses_all):
                file.write(str(coord[0]) + "\t" + str(coord[1]) + "\t" + str(mass) + "\n")
        print("----------")
