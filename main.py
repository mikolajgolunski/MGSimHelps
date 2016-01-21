import MGSimHelps

system = MGSimHelps.AtomsSystem()
print("Reading file.")
# system.readFile("test/water_lammps.dat", "lammps_data", {"lammps_data_type": "charge"})
system.readFile("test/last_frame.lammpstrj", "lammpstrj")
print(system)
system.saveFile("test/out.dat", "lammps_data", {"lammps_data_type": "charge"})
print("Finished.")