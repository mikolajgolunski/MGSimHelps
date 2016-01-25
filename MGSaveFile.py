def saveLammpsDataFile(system, file, control_dict):
    """Save the LAMMPS Data file (extension often .dat).

    :param system: [AtomsSystem] - AtomsSystem object
    :param file: [file] - reference to the file object returned by standard open() method
    :param control_dict: [dict]: - dictionary consisting of additional controls that may be needed
            "lammps_data_type": - type of data (one of)
                "charge"
    """
    if "lammps_data_type" not in control_dict:
        raise NameError("Please provide the data type for LAMMPS data file.")
    print("Saving system into LAMMPS Data file using " + control_dict["lammps_data_type"] + " style.")
    modulo = round(system.number / 10)
    if control_dict["lammps_data_type"] == "charge":
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
        for i, atom in enumerate(system.atoms):
            coord_str = " ".join(str(coord) for coord in atom.coords)
            file.write(str(atom.id) + " " + str(atom.type) + " " + str(atom.charge) + " " + coord_str + "\n")
            if i % modulo == 0:
                print("Written " + str(i) + " out of " + str(system.number) + " atoms.")
        print("Finished writing atoms.")
    else:
        raise NotImplementedError("Data type " + control_dict["lammps_data_type"] + " is not implemented.")


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
    if "items" not in control_dict:
        control_dict["items"] = ["timestep", "time", "number_atoms", "bounds", "atoms"]
        print("[default properties set] ", end="")
    print(" ".join(control_dict["items"]))

    print("Writing header.")
    if "timestep" in control_dict["items"]:
        file.write(
            "ITEM: TIMESTEP\n"
            "" + str(system.timestep) + "\n"
        )
    if "time" in control_dict["items"]:
        file.write(
            "ITEM: TIME\n"
            "" + str(system.time) + "\n"
        )
    if "number_atoms" in control_dict["items"]:
        file.write(
            "ITEM: NUMBER OF ATOMS\n"
            "" + str(system.number) + "\n"
        )
    if "bounds" in control_dict["items"]:
        file.write(
            "ITEM: BOX BOUNDS " + str(system.bounds["xlo"][0]) + str(system.bounds["xhi"][0]) + " " +
                str(system.bounds["ylo"][0]) + str(system.bounds["yhi"][0]) + " " + str(system.bounds["zlo"][0]) +
                str(system.bounds["zhi"][0]) + "\n"
            "" + str(system.bounds["xlo"][1]) + " " + str(system.bounds["xhi"][1]) + "\n"
            "" + str(system.bounds["ylo"][1]) + " " + str(system.bounds["yhi"][1]) + "\n"
            "" + str(system.bounds["zlo"][1]) + " " + str(system.bounds["zhi"][1]) + "\n"
        )
    if "atoms" in control_dict["items"]:
        print("Writing atoms using following properties: ", end="")
        if "properties" not in control_dict:
            control_dict["properties"] = ["id", "element", "type", "x", "y", "z", "charge", "mass", "vx", "vy", "vz",
                                          "kinetic_energy", "potential_energy"]
            print("[default properties set] ", end="")
        else:
            properties_keys = control_dict["properties"]
            properties_str = control_dict["properties"]

        if "charge" in properties_str:
            temp_index = properties_str.index("charge")
            properties_str[temp_index] = "q"
        if "kinetic_energy" in properties_str:
            temp_index = properties_str.index("kinetic_energy")
            properties_str[temp_index] = "c_sput_ke"
        if "potential_energy" in properties_str:
            temp_index = properties_str.index("potential_energy")
            properties_str[temp_index] = "c_sput_pe"
        properties_str = " ".join(properties_str)
        print(properties_str)
        file.write(
            "ITEM: ATOMS "
        )
        file.write(properties_str + "\n")
        modulo = round(system.number / 10)
        for i, atom in enumerate(system.atoms):
            temp_properties = []
            for key in properties_keys:
                if key == "id":
                    property_value = atom.id
                elif key == "element":
                    property_value = atom.name
                elif key == "name":
                    property_value = atom.name
                elif key == "type":
                    property_value = atom.type
                elif key == "x":
                    property_value = atom.coords[0]
                elif key == "y":
                    property_value = atom.coords[1]
                elif key == "z":
                    property_value = atom.coords[2]
                elif key == "q":
                    property_value = atom.charge
                elif key == "mass":
                    property_value = atom.mass
                elif key == "vx":
                    property_value = atom.velocity[0]
                elif key == "vy":
                    property_value = atom.velocity[1]
                elif key == "vz":
                    property_value = atom.velocity[2]
                elif key == "kinetic_energy":
                    property_value = atom.energy["kinetic"]
                elif key == 'potential_energy':
                    property_value = atom.energy["potential"]
                else:
                    raise NotImplementedError("Property " + key + "not implemented in file save.")
                temp_properties.append(str(property_value))
            file.write(" ".join(temp_properties) + "\n")
            if i % modulo == 0:
                print("Written " + str(i) + " out of " + str(system.number) + " atoms.")
        print("Finished writing atoms.")
