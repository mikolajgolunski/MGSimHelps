import itertools
import os.path

import MGReadFile
import MGSaveFile


class Atom:
    """Object containing atoms.

    :param coords: [tuple] - tuple of xyz coordinates:
        [float] - x coordinate
        [float] - y coordinate
        [float] - z coordinate

    As a default contains:
        coords [tuple]:
            [float] - x coordinate
            [float] - y coordinate
            [float] - z coordinate
        id [int] - atom id
        charge [float] - atom charge
        velocity [tuple]:
            [float] - x velocity
            [float] - y velocity
            [float] - z velocity
        type [int] - atom type
        name [string] - atom name
    """

    _new_id = next(itertools.count())

    def __init__(self, coords):
        self.coords = coords
        self.id = Atom._new_id
        self.charge = 0.0
        self.velocity = (0.0, 0.0, 0.0)
        self.type = 1
        self.name = "Xxx"

    def __str__(self):
        return "Atom ID " + str(self.id) + ", type " + str(self.type) + ", name " + self.name + \
               ", at " + str(self.coords)


class AtomsSystem:
    """System of Atom objects.

    As a default contains:
        atoms [list]:
            [Atom] - list of Atom objects
        number [int] - number of atoms in the AtomsSystem
        timestep [int] - timestep
        bounds [dict]:
            "xlo" [tuple]: - lower x boundary
                [string] - boundary type
                [float] - boundary value
            "xhi" [tuple]: - higher x boundary
                [string] - boundary type
                [float] - boundary value
            "ylo" ...
            "yhi" ...
            "zlo" ...
            "zhi" ...
        max_type [int] - maximum type of atom (often equal to number of types)
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

    def readFile(self, file_path, file_type="auto", control_dict=None):
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
        if control_dict is None:
            control_dict = {}
        if file_type == "auto":
            file_ext = os.path.splitext(file_path)[1][1:].strip().lower()
            if file_ext == "lammpstrj":
                file_type = "lammpstrj"
            elif file_ext == "xyz":
                file_type = "xyz"
            else:
                raise NotImplementedError(
                        "No implementation for the " + file_ext + " extension. Please set file_type explicitly."
                )

        with open(file_path, "r") as file:
            if file_type == "lammpstrj":
                MGReadFile.readLammpsFile(self, file)
            elif file_type == "lammps_data":
                MGReadFile.readLammpsDataFile(self, file, control_dict)
            else:
                raise NotImplementedError("No implementation for the " + file_type + " file_type.")

        atoms_names = []
        for atom in self.atoms:
            if atom.name not in atoms_names:
                atoms_names.append(atom.name)
            if atom.type == 1:
                atom.type = atoms_names.index(atom.name) + 1
            if atom.type > self.max_type:
                self.max_type = atom.type

    def saveFile(self, file_path, file_type="auto", control_dict=None):
        """Save file of the set type using data in the AtomsSystem.

        :param file_path: [string] - path to the output file.
        :param file_type: [string] - [Default: "auto"] type of the file. If "auto" than type determined by extension
                of the file specified in file_path:
            "lammps_data" - file read by LAMMPS's function read_data
        :param control_dict: [dict] - [Default: {}] dictionary consisting of additional controls that may be needed:
            "lammps_data_type":
                "charge"
        """
        if control_dict is None:
            control_dict = {}
        if file_type == "auto":
            file_ext = os.path.splitext(file_path)[1][1:].strip().lower()
            if file_ext == "lammpstrj":
                file_type = "lammpstrj"
            elif file_ext == "xyz":
                file_type = "xyz"
            else:
                raise NotImplementedError(
                        "No implementation for the " + file_ext + " extension. Please set file_type explicitly."
                )
        with open(file_path, "w") as file:
            if file_type == "lammps_data":
                MGSaveFile.saveLammpsDataFile(self, file, control_dict)
            else:
                raise NotImplementedError("No implementation for the " + file_type + " file_type.")
