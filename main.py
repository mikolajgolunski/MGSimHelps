import MGSimHelps
import MGTiming
import cProfile

testing = True

path = "test/"

universe = MGSimHelps.Universe()
universe.readFile("test/test.lammpstrj")
# system.readFile("test/water_lammps.dat", "lammps_data", {"lammps_data_type": "charge"})
# system.readFile(
#     path + "5keV.lammpstrj",
#     MGSimHelps.FileType.lammpstrj
# )
for system_name in universe.systems:
    print(system_name)

if testing:
    for counter_system, system in enumerate(universe.systems):
        print("System nr", counter_system)
        # cProfile.run("system.doBinning()")
        system.recalculateBounds()
        system.stretchBounds(x=[-0.1, 0.1], y=[-0.1, 0.1], z=[-0.1, 0.1])
        system.recalculateIDs()
        system.recalculateTypes()
        system.doBinning(new_bounds_q=False)
        system.doCloseNeighbours()
        system.findBonds()
        system.findMolecules()
        system.findMoleculesCOM()
        system.findEjected(height=15)
        system.saveMassSpectrum(file_path="test/mass_spectrum_" + str(counter_system) + ".txt", rounding=0)
else:
    universe.systems[0].saveSystem(
        path + "system.dat",
        MGSimHelps.FileType.lammps_data,
        {MGSimHelps.ControlDict.lammps_data_style: MGSimHelps.LammpsDataStyle.charge}
    )

print("Finished.")
