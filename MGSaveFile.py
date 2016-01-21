def saveLammpsDataFile(system, file, control_dict):
    """Save the LAMMPS Data file (extension often .dat).

    :param system: [AtomsSystem] - AtomsSystem object
    :param file: [file] - reference to the file object returned by standard open() method
    :param control_dict: [dict]: - dictionary consisting of additional controls that may be needed
            "lammps_data_type":
                "charge"
    """
    if "lammps_data_type" not in control_dict:
        raise NameError("Please provide the data type for LAMMPS data file.")
    if control_dict["lammps_data_type"] == "charge":
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
        for atom in system.atoms:
            coord_str = " ".join(str(coord) for coord in atom.coords)
            file.write(str(atom.id) + " " + str(atom.type) + " " + str(atom.charge) + " " + coord_str + "\n")
    else:
        raise NotImplementedError("Data type " + control_dict["lammps_data_type"] + " is not implemented.")
