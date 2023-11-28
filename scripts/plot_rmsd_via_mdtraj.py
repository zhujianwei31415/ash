from ash import *
import mdtraj

#Plotting RMSD of trajectory using mdtraj

#The name of the PDB-file to use for topology
pdbfile="frag-minimized.pdb"
#The name of the trajectory file to load
trajfile="relaxbox_NPT.dcd"

#Loading using mdtraj
system = mdtraj.load(pdbfile)
traj = mdtraj.load(trajfile, top=system)

#Calculating full RMSD
rmsd_all= mdtraj.rmsd(traj, system, 0)
#Defining heavy atoms
heavy_atoms = [atom.index for atom in system.topology.atoms if atom.element.symbol != 'H']
#RMSD for heavy atoms
rmsd_heavy = mdtraj.rmsd(traj, system, 0, atom_indices=heavy_atoms)


#Plotting using ASH-plot (matplotlib)
x_label="Time (ps)"
y_label="RMSD (Å)"
filelabel="RMSD"
eplot = ASH_plot(filelabel, num_subplots=1, x_axislabel=x_label, y_axislabel=y_label)
eplot.addseries(0, x_list=traj.time, y_list=rmsd_heavy, label=y_label, color='blue', line=True, scatter=True)
eplot.savefig(filelabel)
