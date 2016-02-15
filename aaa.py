import MGSimHelps as MG

u = MG.Universe()
s_gen = u.readFile("test/water.dat", MG.FileType.lammps_data, {MG.ControlDict.lammps_data_style: MG.LammpsDataStyle.charge})
system = next(s_gen)
