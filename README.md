# MGSimHelps
Help scripts for my simulations.

What it does:
- reads data from CSV or LAMMPS positions files (or pickle files)
- creates system of atoms based on the data read from files
- loads specific force field fiels (mainly ReaxFF) and calculates bonds between atoms
- based on the data read and calculated allows various analysis such as finding sputtering yields, angular distributions, changes in bonds, molecules' breaking
- writing analysis results in CSV
- writing whole system into pickle files to save all of the calculations
