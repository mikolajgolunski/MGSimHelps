import os
from time import asctime

from tqdm import tqdm

import MGSimHelps

final = True

path = "\\\\ZPGROUP\\ZPGroup\\Mikolaj\\50keVC60_1L\\Lammpstrj"
path_save = "C:\\Users\\mikol\\Desktop\\test\\50keV\\new"

universe = MGSimHelps.Universe()
universe.readReax("potentials/ffield.reax.lg")

if final:

    universe_working = MGSimHelps.Universe()

    files = [os.path.join(path, file) for file_index, file in enumerate(os.listdir(path)) if os.path.isfile(os.path.join(path, file))]
    # files = [os.path.join(path, "50keV_C60_phe1L_0_0.lammpstrj.gz")]

    for counter_files, file_name in tqdm(enumerate(files)):
        for counter_system, system in enumerate(universe.readFile(file_name, control_dict={"start": -1, "end": -1}, gz=True)):
            print(asctime())
            print("System nr", counter_system)
            print(system)
            if system.number != 0:
                # system.correct_coords()
                system.recalculateBounds()
                system.stretchBounds(x=[-0.1, 0.1], y=[-0.1, 0.1], z=[-0.1, 0.1])
                system.recalculateTypes()
                system.doCloseNeighbours(method="few_bins")
                system.reax_params = universe.reax[0]
                system.findBonds(method="reaxFF")
                system.findMolecules()
                system.findMoleculesCOM()
                system.findMoleculesEnergies()
                system.findEjected(height=25)  # TODO: better height
                name = "mass_spectrum_" + '{0:03d}'.format(counter_system) + ".txt"
                system.saveMassSpectrum(os.path.join(path_save, "calc/mass_spectrum/", name), rounding=3)
                name = "mass_" + '{0:03d}'.format(counter_system) + ".txt"
                system.saveMasses(os.path.join(path_save, "calc/masses/", name))
                name = "energy_" + '{0:03d}'.format(counter_system) + ".txt"
                system.saveEnergies(os.path.join(path_save, "calc/energies/", name))
                name = "spot_plot_" + '{0:03d}'.format(counter_system) + ".txt"
                system.saveSpotPlot(os.path.join(path_save, "calc/spot_plots/", name))
                name = "angle_dist_" + '{0:03d}'.format(counter_system) + ".txt"
                system.save_angle_distribution(os.path.join(path_save, "calc/angle_dist/", name))

                # universe_working.systems = [None, system]
                # for system_previous in universe.readFile(file_name, control_dict={"start": 1, "end": 1}, gz=True):
                #     system_previous.correct_coords()
                #     system_previous.recalculateBounds()
                #     system_previous.stretchBounds(x=[-0.1, 0.1], y=[-0.1, 0.1], z=[-0.1, 0.1])
                #     system_previous.recalculateTypes()
                #     system_previous.doCloseNeighbours(method="few_bins")
                #     system_previous.reax_params = universe.reax[0]
                #     system_previous.findBonds(method="reaxFF")
                #     system_previous.findMolecules()
                #     system_previous.findMoleculesCOM()
                #     system_previous.findMoleculesEnergies()
                #     system_previous.findEjected(height=25)
                #
                #     universe_working.systems[0] = system_previous
                #     break
                #
                # name = "emitted_coords_full_" + '{0:03d}'.format(counter_files) + ".txt"
                # path_full = os.path.join(path_save, "calc/emitted_coords/", name)
                # name = "emitted_coords_partial_" + '{0:03d}'.format(counter_files) + ".txt"
                # path_partial = os.path.join(path_save, "calc/emitted_coords/", name)
                # name = "energy_coords_" + '{0:03d}'.format(counter_files) + ".txt"
                # path_energies = os.path.join(path_save, "calc/emitted_coords/energies/", name)
                # universe_working.saveEmittedCoords(path_full, path_partial, path_energies)
else:
    file_name = os.path.join(path, "0_5keV_C60_phe1L_0_5.lammpstrj.gz")
    for counter_system, system in enumerate(
            universe.readFile(file_name, control_dict={"start": -1, "end": -1}, gz=True)):
        print("System nr", counter_system)
        print(system)
        system.recalculateBounds()
        system.stretchBounds(x=[-0.1, 0.1], y=[-0.1, 0.1], z=[-0.1, 0.1])
        system.recalculateTypes()
        system.doCloseNeighbours(method="few_bins")
        system.reax_params = universe.reax[0]
        system.findBonds(method="reaxFF")
        system.findMolecules()
        system.findMoleculesCOM()
        system.findMoleculesEnergies()
        system.findEjected(height=25)  # TODO: better height
        name = "mass_spectrum_" + '{0:03d}'.format(counter_system) + ".txt"
        system.saveMassSpectrum(os.path.join(path_save, "calc/mass_spectrum/", name), rounding=3)
        name = "mass_" + '{0:03d}'.format(counter_system) + ".txt"
        system.saveMasses(os.path.join(path_save, "calc/masses/", name))
        name = "energy_" + '{0:03d}'.format(counter_system) + ".txt"
        system.saveEnergies(os.path.join(path_save, "calc/energies/", name))
        name = "spot_plot_" + '{0:03d}'.format(counter_system) + ".txt"
        system.saveSpotPlot(os.path.join(path_save, "calc/spot_plots/", name))
        name = "angle_dist_" + '{0:03d}'.format(counter_system) + ".txt"
        system.save_angle_distribution(os.path.join(path_save, "calc/angle_dist/", name))
