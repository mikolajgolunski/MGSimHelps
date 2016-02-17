import MGSimHelps
import MGTiming
import cProfile

testing = True

path = "test/"

# system.readFile("test/water_lammps.dat", "lammps_data", {"lammps_data_type": "charge"})
# system.readFile(
#     path + "5keV.lammpstrj",
#     MGSimHelps.FileType.lammpstrj
# )

universe = MGSimHelps.Universe()
universe.readReax("test/ffield.reax.lg")
for counter_system, system in enumerate(universe.readFile("test/test.lammpstrj")):
    print("System nr", counter_system)
    print(system)
    if testing:
        # cProfile.run("system._doBinning()")
        system.recalculateBounds()
        system.stretchBounds(x=[-0.1, 0.1], y=[-0.1, 0.1], z=[-0.1, 0.1])
        system.recalculateIDs()
        system.recalculateTypes()
        system.doCloseNeighbours(method="few_bins")
        system.reax_params = universe.reax[0]
        system.findBonds(method="reaxFF")
        system.findMolecules()
        system.findMoleculesCOM()
        system.findMoleculesKineticEnergy()
        system.findEjected(height=15)
        system.saveMassSpectrum("test/output/mass_spectra/mass_spectrum_" + str(counter_system) + ".txt", rounding=0)
        system.saveMasses("test/output/masses/mass_" + str(counter_system) + ".txt")
        system.saveEnergies("test/output/energies/energy_" + str(counter_system) + ".txt")
        system.saveSpotPlot("test/output/spot_plots/spot_plot_" + str(counter_system) + ".txt")
    else:
        system.saveSystem(
            path + "system.dat",
            MGSimHelps.FileType.lammps_data,
            {MGSimHelps.ControlDict.lammps_data_style: MGSimHelps.LammpsDataStyle.charge}
        )

print("Finished.")
