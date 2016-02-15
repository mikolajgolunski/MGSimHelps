import MGSimHelps as MG
from tqdm import tqdm


def saveLammpsDataFile(system, file, control_dict):
    """Save the LAMMPS Data file (extension often .dat).

    :param system: [AtomsSystem] - AtomsSystem object
    :param file: [file] - reference to the file object returned by standard open() method
    :param control_dict: [dict]: - dictionary consisting of additional controls that may be needed
            "lammps_data_type": - type of data (one of)
                "charge"
    """
    if MG.ControlDict.lammps_data_style not in control_dict:
        raise NameError("Please provide the data style for LAMMPS data file.")
    print("Saving system into LAMMPS Data file using " + control_dict[MG.ControlDict.lammps_data_style].name +
          " style.")
    if control_dict[MG.ControlDict.lammps_data_style] == MG.LammpsDataStyle.charge:
        print("Writing header.")
        file.write(
            "LAMMPS data file. CGCMM style. atom_style charge. Converted by Mikolaj Golunski (MGSimHelp utility).\n"
            " " + str(system.number) + " atoms\n"
            " " + str(system.max_type) + " atom types\n"
            " " + str(system.bounds["xlo"][1]) + " " + str(system.bounds["xhi"][1]) + "  xlo xhi\n"
            " " + str(system.bounds["ylo"][1]) + " " + str(system.bounds["yhi"][1]) + "  ylo yhi\n"
            " " + str(system.bounds["zlo"][1]) + " " + str(system.bounds["zhi"][1]) + "  zlo zhi\n"
            "\n"
            " Atoms\n"
            "\n"
        )
        for i, atom in tqdm(enumerate(system.atoms), total=system.number, unit="atom"):
            coord_str = " ".join(str(coord) for coord in atom.coords)
            file.write(str(atom.id) + " " + str(atom.type) + " " + str(atom.charge) + " " + coord_str + "\n")
        print("\nFinished writing atoms.")
    else:
        raise NotImplementedError("Data type " + control_dict[MG.ControlDict.lammps_data_style] +
                                  " is not implemented.")


def saveLammpsFile(system, file, control_dict):
    """Save the LAMMMPS trajectory file (.lammpstrj).

    :param system: [AtomsSystem] - AtomsSyste object
    :param file: [file] - reference to the file object returned by standard open() method
    :param control_dict: [dict]: - dictionary consisting of additional controls that may be needed
        "properties": - list of properties (any of)
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
        "items": - list of items to save (any of)
            "timestep"
            "time"
            "number_atoms"
            "bounds"
            "atoms"
    """
    print("Saving LAMPPS trajectory file with following items: ", end="")
    if MG.ControlDict.lammpstrj_items not in control_dict:
        control_dict[MG.ControlDict.lammpstrj_items] = [
            MG.LammpstrjItem.timestep, MG.LammpstrjItem.time, MG.LammpstrjItem.number_atoms, MG.LammpstrjItem.bounds,
            MG.LammpstrjItem.atoms
        ]
        print("[default properties set] ", end="")
    print(" ".join([str(item.name) for item in control_dict[MG.ControlDict.lammpstrj_items]]))

    print("Writing header.")
    if MG.LammpstrjItem.timestep in control_dict[MG.ControlDict.lammpstrj_items]:
        file.write(
            "ITEM: TIMESTEP\n"
            "" + str(system.timestep) + "\n"
        )
    if MG.LammpstrjItem.time in control_dict[MG.ControlDict.lammpstrj_items]:
        file.write(
            "ITEM: TIME\n"
            "" + str(system.time) + "\n"
        )
    if MG.LammpstrjItem.number_atoms in control_dict[MG.ControlDict.lammpstrj_items]:
        file.write(
            "ITEM: NUMBER OF ATOMS\n"
            "" + str(system.number) + "\n"
        )
    if MG.LammpstrjItem.bounds in control_dict[MG.ControlDict.lammpstrj_items]:
        file.write(
            "ITEM: BOX BOUNDS " + str(system.bounds["xlo"][0]) + str(system.bounds["xhi"][0]) + " " +
                str(system.bounds["ylo"][0]) + str(system.bounds["yhi"][0]) + " " + str(system.bounds["zlo"][0]) +
                str(system.bounds["zhi"][0]) + "\n"
            "" + str(system.bounds["xlo"][1]) + " " + str(system.bounds["xhi"][1]) + "\n"
            "" + str(system.bounds["ylo"][1]) + " " + str(system.bounds["yhi"][1]) + "\n"
            "" + str(system.bounds["zlo"][1]) + " " + str(system.bounds["zhi"][1]) + "\n"
        )
    if MG.LammpstrjItem.atoms in control_dict[MG.ControlDict.lammpstrj_items]:
        print("Writing atoms using following properties: ", end="")
        if MG.ControlDict.lammpstrj_properties not in control_dict:
            control_dict[MG.ControlDict.lammpstrj_properties] = [
                MG.LammpstrjProp.id, MG.LammpstrjProp.element, MG.LammpstrjProp.type, MG.LammpstrjProp.x,
                MG.LammpstrjProp.y, MG.LammpstrjProp.z, MG.LammpstrjProp.charge, MG.LammpstrjProp.mass,
                MG.LammpstrjProp.vx, MG.LammpstrjProp.vy, MG.LammpstrjProp.vz, MG.LammpstrjProp.kinetic_energy,
                MG.LammpstrjProp.potential_energy
            ]
            print("[default properties set] ", end="")

        properties_keys = [str(key.name) for key in control_dict[MG.ControlDict.lammpstrj_properties]]
        properties_str = [str(key.name) for key in control_dict[MG.ControlDict.lammpstrj_properties]]

        if str(MG.LammpstrjProp.charge.name) in properties_str:
            temp_index = properties_str.index(str(MG.LammpstrjProp.charge.name))
            properties_str[temp_index] = "q"
        if str(MG.LammpstrjProp.kinetic_energy.name) in properties_str:
            temp_index = properties_str.index(str(MG.LammpstrjProp.kinetic_energy.name))
            properties_str[temp_index] = "c_sput_ke"
        if str(MG.LammpstrjProp.potential_energy.name) in properties_str:
            temp_index = properties_str.index(str(MG.LammpstrjProp.potential_energy.name))
            properties_str[temp_index] = "c_sput_pe"
        properties_str = " ".join(properties_str)
        print(properties_str)
        file.write(
            "ITEM: ATOMS "
        )
        file.write(properties_str + "\n")
        for i, atom in tqdm(enumerate(system.atoms), total=system.number, unit="atom"):
            temp_properties = []
            for key in properties_keys:
                if key == MG.LammpstrjProp.id:
                    property_value = atom.id
                elif key == MG.LammpstrjProp.element:
                    property_value = atom.name
                elif key == MG.LammpstrjProp.name:
                    property_value = atom.name
                elif key == MG.LammpstrjProp.type:
                    property_value = atom.type
                elif key == MG.LammpstrjProp.x:
                    property_value = atom.coords[0]
                elif key == MG.LammpstrjProp.y:
                    property_value = atom.coords[1]
                elif key == MG.LammpstrjProp.z:
                    property_value = atom.coords[2]
                elif key == MG.LammpstrjProp.q:
                    property_value = atom.charge
                elif key == MG.LammpstrjProp.mass:
                    property_value = atom.mass
                elif key == MG.LammpstrjProp.vx:
                    property_value = atom.velocity[0]
                elif key == MG.LammpstrjProp.vy:
                    property_value = atom.velocity[1]
                elif key == MG.LammpstrjProp.vz:
                    property_value = atom.velocity[2]
                elif key == MG.LammpstrjProp.kinetic_energy:
                    property_value = atom.energy["kinetic"]
                elif key == MG.LammpstrjProp.potential_energy:
                    property_value = atom.energy["potential"]
                else:
                    raise NotImplementedError("Property " + str(key) + "not implemented in file save.")
                temp_properties.append(str(property_value))
            file.write(" ".join(temp_properties) + "\n")
        print("\nFinished writing atoms.")
