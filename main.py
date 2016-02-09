import MGSimHelps

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
    system.recalculateBounds()

    # system.recalculateIDs()

    system.doBinning()

    system.doCloseNeighbours("test/neighbours.dat")

    # system.findMolecules()
else:
    system.saveSystem(
            path + "system.dat",
            MGSimHelps.FileType.lammps_data,
            {MGSimHelps.ControlDict.lammps_data_style: MGSimHelps.LammpsDataStyle.charge}
    )

print("Finished.")
