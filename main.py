from enum import Enum

import MGSimHelps

system = MGSimHelps.AtomsSystem()
# system.readFile("test/water_lammps.dat", "lammps_data", {"lammps_data_type": "charge"})
system.readFile("test/5keV.lammpstrj", "lammpstrj")
print(system)

system.recalculateBounds()

system.saveFile("test/system.dat", "lammps_data", {"lammps_data_type": "charge"})

print("Finished.")
