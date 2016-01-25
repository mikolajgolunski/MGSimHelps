import MGSimHelps

system = MGSimHelps.AtomsSystem()
# system.readFile("test/water_lammps.dat", "lammps_data", {"lammps_data_type": "charge"})
system.readFile(
        "test/5keV.lammpstrj",
        MGSimHelps.FileType.auto
)
print(system)

system.recalculateBounds()

system.saveFile(
        "test/system.dat",
        MGSimHelps.FileType.lammps_data,
        {MGSimHelps.ControlDict.lammps_data_style: MGSimHelps.LammpsDataStyle.charge}
)

print("Finished.")
