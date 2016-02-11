import MGSimHelps
import cProfile

testing = True

path = "test/"

system = MGSimHelps.AtomsSystem()
# system.readFile("test/water_lammps.dat", "lammps_data", {"lammps_data_type": "charge"})
system.readFile(
    path + "5keV.lammpstrj",
    MGSimHelps.FileType.lammpstrj
)
print(system)

if testing:
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
    system.saveMassSpectrum(file_path="test/mass_spectrum.txt", rounding=0)
else:
    system.saveSystem(
        path + "system.dat",
        MGSimHelps.FileType.lammps_data,
        {MGSimHelps.ControlDict.lammps_data_style: MGSimHelps.LammpsDataStyle.charge}
    )

print("Finished.")
