import pickle
import os.path
import gzip

from tqdm import tqdm

import MGSimHelps
import MGTiming
# import cProfile

testing = True
read_archive = False

# path = "C:/Users/mikol/Documents/simulations/Phe/C60/phe1L/50keV"
path = "C:/Users/mikol/Dropbox/ZPGROUP/Wyniki/Mikolaj/50keVC60/lammpstrj"
path_save = "C:/Users/mikol/Documents/simulations/Phe/C60/phe1L/50keV/many_hits"

# system.readFile("test/water_lammps.dat", "lammps_data", {"lammps_data_type": "charge"})
# system.readFile(
#     path + "5keV.lammpstrj",
#     MGSimHelps.FileType.lammpstrj
# )

universe = MGSimHelps.Universe()
universe.readReax("potentials/ffield.reax.lg")
import csv
if read_archive:
    path_archive = os.path.join(path, "archive")
    files = [os.path.join(path_archive, file) for file_index, file in enumerate(os.listdir(path_archive))
             if os.path.isfile(os.path.join(path_archive, file))]
    for counter_files, file_name in enumerate(files):
        if counter_files > 145:
            with gzip.open(file_name) as file:
                system = pickle.load(file)
                name = "dist_" + '{0:03d}'.format(counter_files) + ".txt"
                system.save_angle_distribution(os.path.join(path, "calc/dist/", name))
    # with open(os.path.join(path, "C60_phe1L_500eV.lammps"), "w", newline="") as f:
    #     csvwriter = csv.writer(f, delimiter="\t", quoting=csv.QUOTE_NONE, quotechar="", escapechar="")
    #     for counter_system, path_file in tqdm(enumerate(files)):
    #         with gzip.open(path_file, "rb") as file:
    #             system = pickle.load(file)
    #         csvwriter.writerow(["ITEM:", "TIMESTEP"])
    #         csvwriter.writerow([system.timestep])
    #         csvwriter.writerow(["ITEM:", "NUMBER OF ATOMS"])
    #         csvwriter.writerow([system.number])
    #         csvwriter.writerow(["ITEM:", "BOX BOUNDS ff ff ff"])
    #         csvwriter.writerow([system.bounds["xlo"][1], system.bounds["xhi"][1]])
    #         csvwriter.writerow([system.bounds["ylo"][1], system.bounds["yhi"][1]])
    #         csvwriter.writerow([system.bounds["zlo"][1], system.bounds["zhi"][1]])
    # \\ZPGROUP\ZPGroup\Mikolaj\500eVC60_1L\Lammpstrj        csvwriter.writerow(["ITEM:", "ATOMS", "id", "element", "x", "y", "z", "mass", "vx", "vy", "vz"])
    #         for molecule in system.molecules:
    #             for atom in molecule.atoms:
    #                 csvwriter.writerow([atom.id, atom.name, atom.coords[0], atom.coords[1], atom.coords[2],
    #                                     atom.mass, molecule.centre_of_mass["v"][0], molecule.centre_of_mass["v"][1],
    #                                     molecule.centre_of_mass["v"][2]])
            # name = "helper_" + '{0:03d}'.format(counter_system) + "_"
            # system.saveHelper(os.path.join(path, "calc/helper/", name))
            # name = "mass_spectrum_" + '{0:03d}'.format(counter_system) + ".txt"
            # system.saveMassSpectrum(os.path.join(path, "calc/mass_spectrum/", name), rounding=3)
            # name = "mass_" + '{0:03d}'.format(counter_system) + ".txt"
            # system.saveMasses(os.path.join(path, "calc/masses/", name))
            # name = "energy_" + '{0:03d}'.format(counter_system) + ".txt"
            # system.saveEnergies(os.path.join(path, "calc/energies/", name))
            # name = "spot_plot_" + '{0:03d}'.format(counter_system) + ".txt"
            # system.saveSpotPlot(os.path.join(path, "calc/spot_plots/", name))
            del system
else:
    path_archive = os.path.join(path_save, "archive")
    files = [os.path.join(path, file) for file_index, file in enumerate(os.listdir(path))
             if os.path.isfile(os.path.join(path, file))]
    for counter_files, file_name in tqdm(enumerate(files)):
        for counter_system, system in enumerate(universe.readFile(file_name, control_dict={"start": -1, "end": -1},
                                                                  gz=True)):
            print("System nr", counter_files)
            print(system)
            if testing:
                # cProfile.run("system._doBinning()")
                system.recalculateBounds()
                system.stretchBounds(x=[-0.1, 0.1], y=[-0.1, 0.1], z=[-0.1, 0.1])
                # system.recalculateIDs()
                system.recalculateTypes()
                system.doCloseNeighbours(method="few_bins")
                system.reax_params = universe.reax[0]
                system.findBonds(method="reaxFF")
                system.findMolecules()
                system.findMoleculesCOM()
                system.findMoleculesEnergies()
                system.findEjected(height=25)  # TODO: better height
                name = "mass_spectrum_" + '{0:03d}'.format(counter_files) + ".txt"
                system.saveMassSpectrum(os.path.join(path_save, "calc/mass_spectrum/", name), rounding=3)
                name = "mass_" + '{0:03d}'.format(counter_files) + ".txt"
                system.saveMasses(os.path.join(path_save, "calc/masses/", name))
                name = "energy_" + '{0:03d}'.format(counter_files) + ".txt"
                system.saveEnergies(os.path.join(path_save, "calc/energies/", name))
                name = "spot_plot_" + '{0:03d}'.format(counter_files) + ".txt"
                system.saveSpotPlot(os.path.join(path_save, "calc/spot_plots/", name))
                name = "angle_dist_" + '{0:03d}'.format(counter_files) + ".txt"
                system.save_angle_distribution(os.path.join(path_save, "calc/angle_dist/", name))

                # system.archiveSystem(path_save)
            else:
                system.saveSystem(
                    os.path.join(path, "system.dat"),
                    MGSimHelps.FileType.lammps_data,
                    {MGSimHelps.ControlDict.lammps_data_style: MGSimHelps.LammpsDataStyle.charge}
                )

print("Finished.")
