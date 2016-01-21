import MGSimHelps


def readLammpsFile(system, file):  # TODO Implement additional information that can be found in file
    """Read content of the LAMMPS Trajectory file (extension .lammpstrj).

    :param system: [AtomsSystem] - AtomsSystem object
    :param file: [file] - reference to the file object returned by standard open() method
    """
    axis = 0
    for line in file:
        keyvalue = line.strip().partition(": ")
        if keyvalue[0] == "ITEM":
            if keyvalue[2].isupper():
                if keyvalue[2] == "TIMESTEP":
                    command = "timestep"
                elif keyvalue[2] == "TIME":
                    command = "time"
                elif keyvalue[2] == "NUMBER OF ATOMS":
                    command = "nr_atoms"
                else:
                    raise NotImplementedError("Unknown data header in LAMMPS file: " + keyvalue[2])
            else:
                split = keyvalue[2].split()
                key_elements = []
                value_elements = []
                for element in split:
                    if element.isupper():
                        key_elements.append(element)
                    else:
                        value_elements.append(element)
                key = " ".join(key_elements)
                if key == "BOX BOUNDS":
                    command = "box_bounds"
                elif key == "ATOMS":
                    command = "atoms"
                else:
                    raise NotImplementedError("Unknown data header in LAMMPS file: " + keyvalue[2])
        else:
            if command == "timestep":
                system.timestep = int(keyvalue[0])
            elif command == "time":
                system.time = float(keyvalue[0])
            elif command == "nr_atoms":
                pass  # nr of atoms
            elif command == "box_bounds":
                axis += 1
                split = keyvalue[0].split()
                if axis == 1:
                    system.bounds["xlo"] = (value_elements[0][0], float(split[0]))
                    system.bounds["xhi"] = (value_elements[0][1], float(split[1]))
                elif axis == 2:
                    system.bounds["ylo"] = (value_elements[1][0], float(split[0]))
                    system.bounds["yhi"] = (value_elements[1][1], float(split[1]))
                elif axis == 3:
                    system.bounds["zlo"] = (value_elements[2][0], float(split[0]))
                    system.bounds["zhi"] = (value_elements[2][1], float(split[1]))
                else:
                    raise NotImplementedError("Too many numbers in BOX BOUNDS part of the LAMMPS file.")
            elif command == "atoms":
                split = keyvalue[0].split()
                temp_dict = {}
                for key, value in zip(value_elements, split):
                    temp_dict[key] = value
                temp_coords = [0.0, 0.0, 0.0]
                if "x" in temp_dict:
                    temp_coords[0] = float(temp_dict["x"])
                if "y" in temp_dict:
                    temp_coords[1] = float(temp_dict["y"])
                if "z" in temp_dict:
                    temp_coords[2] = float(temp_dict["z"])
                atom = MGSimHelps.Atom(tuple(temp_coords))
                if "id" in temp_dict:
                    atom.id = int(temp_dict["id"])
                if "type" in temp_dict:
                    atom.type = int(temp_dict["type"])
                if "name" in temp_dict:
                    atom.name = temp_dict["name"]
                if "element" in temp_dict:
                    atom.name = temp_dict["element"]
                if "q" in temp_dict:
                    atom.charge = float(temp_dict["q"])
                if "vx" in temp_dict or "vy" in temp_dict or "vz" in temp_dict:
                    temp_v = [0.0, 0.0, 0.0]
                    if "vx" in temp_dict:
                        temp_v[0] = float(temp_dict["vx"])
                    if "vy" in temp_dict:
                        temp_v[1] = float(temp_dict["vy"])
                    if "vz" in temp_dict:
                        temp_v[2] = float(temp_dict["vz"])
                    atom.velocity = tuple(temp_v)
                if "mass" in temp_dict:
                    atom.mass = float(temp_dict["mass"])
                if "c_sput_ke" in temp_dict:
                    atom.energy["kinetic"] = float(temp_dict["c_sput_ke"])
                if "c_sput_pe" in temp_dict:
                    atom.energy["potential"] = float(temp_dict["c_sput_pe"])
                system.addAtom(atom)
            else:
                raise NameError("Unknown command: " + command)


def readLammpsDataFile(system, file, control_dict):
    """Read content of the LAMMPS Data file (extension often .dat).

    :param system: [AtomsSystem] - AtomsSystem object
    :param file: [file] - reference to the file object returned by standard open() method
    :param control_dict: [dict]: - dictionary consisting of additional controls that may be needed
            "lammps_data_type":
                "charge"
    """
    atoms_flag = False
    if "lammps_data_type" not in control_dict:
        raise NameError("Please provide the data type for LAMMPS data file.")
    if control_dict["lammps_data_type"] == "charge":
        for line in file:
            content = line.strip().split()
            if "atoms" in content:
                pass  # nr of atoms
            elif "types" in content:
                system.max_type = int(content[0])
            elif "xlo" in content:
                system.bounds["xlo"] = (system.bounds["xlo"][0], float(content[0]))
                system.bounds["xhi"] = (system.bounds["xhi"][0], float(content[1]))
            elif "ylo" in content:
                system.bounds["ylo"] = (system.bounds["ylo"][0], float(content[0]))
                system.bounds["yhi"] = (system.bounds["yhi"][0], float(content[1]))
            elif "zlo" in content:
                system.bounds["zlo"] = (system.bounds["zlo"][0], float(content[0]))
                system.bounds["zhi"] = (system.bounds["zhi"][0], float(content[1]))
            elif "Atoms" in content and len(content) == 1:
                atoms_flag = True
            else:
                if atoms_flag and len(content) > 1:
                    atom = MGSimHelps.Atom((float(content[3]), float(content[4]), float(content[5])))
                    atom.id = int(content[0])
                    atom.type = int(content[1])
                    atom.charge = float(content[2])
                    system.addAtom(atom)
    else:
        raise NotImplementedError("Data type " + control_dict["lammps_data_type"] + " is not implemented.")
