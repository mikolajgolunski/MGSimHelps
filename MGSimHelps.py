import math
import itertools
import os.path
from enum import Enum
import numpy as np

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
        self.bonds = set()

    def __str__(self):
        return "Atom ID " + str(self.id) + ", type " + str(self.type) + ", name " + self.name + \
               ", at " + str(self.coords)


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
        self.bins = []

    def __iter__(self):
        return self

    def next(self):
        iter_id = self._next_iter_id
        if iter_id > self.number:
            self._next_iter_id = next(itertools.count())
            raise StopIteration
        return self._atoms[iter_id]

    def __str__(self):
        return "AtomsSystem with " + str(self.number) + " atoms."

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
        print("Reading file.")
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

        with open(file_path, "r") as file:
            if file_type == FileType.lammpstrj:
                MGReadFile.readLammpsFile(self, file)
            elif file_type == FileType.lammps_data:
                MGReadFile.readLammpsDataFile(self, file, control_dict)
            else:
                raise NotImplementedError("No implementation for the " + str(file_type.name) + " file_type.")
        print("Successfully read system from " + file_path + " as a " + str(file_type.name) + " file.")

        self.recalculateTypes()

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
        print("Saving system.")
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
        with open(file_path, "w") as file:
            if file_type == FileType.lammps_data:
                MGSaveFile.saveLammpsDataFile(self, file, control_dict)
            elif file_type == FileType.lammpstrj:
                MGSaveFile.saveLammpsFile(self, file, control_dict)
            else:
                raise NotImplementedError("No implementation for the " + str(file_type.name) + " file_type.")
        print("Successfully saved system to " + file_path + " as a " + str(file_type.name) + " file.")

    def recalculateBounds(self):
        """Recalculate boundaries of the system."""
        print("Recalculating boundaries of the system.")
        offset = 0.1
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
        self.bounds["xhi"] = (self.bounds["xhi"][0], max_coords[0] + offset)
        self.bounds["xlo"] = (self.bounds["xlo"][0], min_coords[0] - offset)
        self.bounds["yhi"] = (self.bounds["yhi"][0], max_coords[1] + offset)
        self.bounds["ylo"] = (self.bounds["ylo"][0], min_coords[1] - offset)
        self.bounds["zhi"] = (self.bounds["zhi"][0], max_coords[2] + offset)
        self.bounds["zlo"] = (self.bounds["zlo"][0], min_coords[2] - offset)

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
        self.bounds["xlo"] = (self.bounds["xlo"][0], self.bounds["xlo"][1] + x[0])
        self.bounds["xhi"] = (self.bounds["xhi"][0], self.bounds["xhi"][1] + x[1])
        self.bounds["ylo"] = (self.bounds["ylo"][0], self.bounds["ylo"][1] + y[0])
        self.bounds["yhi"] = (self.bounds["yhi"][0], self.bounds["yhi"][1] + y[1])
        self.bounds["zlo"] = (self.bounds["zlo"][0], self.bounds["zlo"][1] + z[0])
        self.bounds["zhi"] = (self.bounds["zhi"][0], self.bounds["zhi"][1] + z[1])

    def recalculateIDs(self):
        """Recalculate IDs of atoms in the system"""
        print("Recalculating IDs.")
        for i, atom in enumerate(self.atoms):
            atom.id = i

    def recalculateTypes(self):
        """Recalculate types of atoms."""
        print("Recalculating types of atoms.")
        atoms_names = []
        for atom in self.atoms:
            if atom.name not in atoms_names:
                atoms_names.append(atom.name)
            if atom.type == 1:
                atom.type = atoms_names.index(atom.name) + 1
            if atom.type > self.max_type:
                self.max_type = atom.type

    def doBinning(self):  # TODO: dynamic bins (changing size depending on nr of atoms in bin
        """Binning atoms."""
        self.recalculateBounds()
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

        self.bins = [[[np.array([], dtype=np.int_)]*math.ceil(deltaz / dz) for y in range(math.ceil(deltay / dy))]
                     for x in range(math.ceil(deltax / dx))]
        print("Created " + str(math.ceil(deltaz / dz)) + " x " + str(math.ceil(deltay / dy)) + " x " +
              str(math.ceil(deltay / dy)) + " bins.")
        modulo = round(self.number / 10)
        for i, atom in enumerate(self.atoms):
            binx = math.floor((atom.coords[0] - self.bounds["xlo"][1]) / dx)
            biny = math.floor((atom.coords[1] - self.bounds["ylo"][1]) / dy)
            binz = math.floor((atom.coords[2] - self.bounds["zlo"][1]) / dz)
            atom.bin = [binx, biny, binz]
            self.bins[binx][biny][binz] = np.append(self.bins[binx][biny][binz], i)
            if i % modulo == 0:
                print("Binned " + str(i) + " out of " + str(self.number) + " atoms.")

    def doCloseNeighbours(self):
        """Create list of closest neighbours for each atom in the system."""
        if len(self.bins) == 0:
            self.doBinning()
        modulo = round(self.number / 10)
        xmax = len(self.bins) - 1
        ymax = len(self.bins[0]) - 1
        zmax = len(self.bins[0][0]) - 1
        for atomCounter, atom in enumerate(self.atoms):
            x = atom.bin[0]
            y = atom.bin[1]
            z = atom.bin[2]
            big_bin = [neigh for neigh in self.bins[x][y][z] if neigh != atomCounter]
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
            if atomCounter % modulo == 0:
                print("Found neighbours for " + str(atomCounter) + " out of " + str(self.number) + " atoms.")

    def findBonds(self, method=None):  # TODO: Add additional methods of calculating bonds
        """Find bonds between atoms in the system.

        :param method: [string] - method used
        """
        modulo = round(self.number / 10)
        if method is None:
            for atomCounter, atom in enumerate(self.atoms):
                for neighbour in atom.neighbours:
                    distance = math.sqrt((atom.coords[0] - self.atoms[neighbour].coords[0])**2 +
                                         (atom.coords[1] - self.atoms[neighbour].coords[1])**2 +
                                         (atom.coords[2] - self.atoms[neighbour].coords[2])**2)
                    if distance < 2.0:
                        atom.bonds.add(neighbour)
                if atomCounter % modulo == 0:
                    print("Found bonds for " + str(atomCounter) + " out of " + str(self.number) + " atoms.")

    def findMolecules(self):
        """Find molecules in the system."""
        pass