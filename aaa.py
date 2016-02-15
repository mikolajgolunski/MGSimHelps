import MGSimHelps as MG
#import MGTiming

u = MG.Universe()
# u.readFile("test/water.dat", file_type=MG.FileType.lammps_data, control_dict={MG.ControlDict.lammps_data_style: MG.LammpsDataStyle.charge})
u.readFile("test/500eV.lammpstrj")
# u.readFile("test/500eV.lammpstrj")
# s_no = u.systems[0]
# s_no.doCloseNeighbours(method="no_bins")
# s_no.findBonds()
# s_no.findMolecules()
s_few = u.systems[0]
s_few.doCloseNeighbours(method="few_bins")
# s_few.findBonds()
# s_few.findMolecules()