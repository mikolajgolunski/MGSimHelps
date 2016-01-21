import MGSimHelps

system = MGSimHelps.AtomsSystem()
print("Reading file.")
# system.readFile("test/water_lammps.dat", "lammps_data", {"lammps_data_type": "charge"})
system.readFile("test/500eV.lammpstrj", "lammpstrj")
print(system)
system.saveFile("test/out.dat", "lammps_data", {"lammps_data_type": "charge"})
print("Finished.")
