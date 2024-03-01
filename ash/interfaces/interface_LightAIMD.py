import numpy as np
import shutil
import os
import subprocess as sp
import glob
import time

import ash.modules.module_coords
from ash.functions.functions_general import ashexit, BC,print_time_rel


#LightAIMD Theory object.
class LightAIMDTheory:
    def __init__(self, lightaimddir=None, settings=None,
                 filename='lightaimd_', label=None, numcores=1, printlevel=2):
        # label
        self.theorynamelabel = 'LightAIMD'
        self.theorytype = 'QM'
        self.analytic_hessian = False

        #Printlevel
        self.printlevel = printlevel
        self.numcores = numcores
        self.label = label
        self.filename = filename

        #Determining binary path
        print("Defining LightAIMD object")
        if lightaimddir is not None:
            print("Path to LightAIMD provided:", lightaimddir)
            self.lightaimdpath = os.path.join(lightaimddir, 'lightaimd')
        else:
            self.lightaimdpath = shutil.which('lightaimd')
            if self.lightaimdpath == None:
                print("Found no lightaimd in path. Add LightAIMD to Shell environment or provide lightaimddir variable")
                ashexit()
            else:
                print("Found lightaimd in path:", self.lightaimdpath)

        #Settings dict
        if settings == None:
            print("settings dict not set. Exiting")
            ashexit()
        self.settings = settings

    #Set numcores
    def set_numcores(self, numcores):
        self.numcores = numcores

    #Cleanup after run.
    def cleanup(self):
        print("Cleaning up old LightAIMD files")
        try:
            os.remove('timer.dat')
            os.remove(self.filename+'.out')
        except:
            pass

    #Run function. Takes coords, elems etc. arguments and computes E or E+G.
    def run(self, current_coords=None, current_MM_coords=None, MMcharges=None,
            qm_elems=None, label=None, mm_elems=None, elems=None,
            Grad=False, PC=False, numcores=None, restart=False,
            charge=None, mult=None):

        module_init_time=time.time()

        if numcores is None:
            numcores = self.numcores

        print(BC.OKBLUE, BC.BOLD, "------------RUNNING LIGHTAIMD INTERFACE-------------", BC.END)

        #Checking if charge and mult has been provided
        if charge == None or mult == None:
            print(BC.FAIL, "Error. charge and mult has not been defined for LightAIMDTheory.run method", BC.END)
            ashexit()

        #Coords provided to run
        if current_coords is not None:
            pass
        else:
            print("no current_coords")
            ashexit()

        #What elemlist to use. If qm_elems provided then QM/MM job, otherwise use elems list
        if qm_elems is None:
            if elems is None:
                print("No elems provided")
                ashexit()
            else:
                qm_elems = elems

        print("LightAIMD Running")
        print("Current directory:", os.getcwd())
        print("LightAIMD settings:", self.settings)
       
        # prepare command line parameters
        cmds = []

        # mpi setting
        if 'mpicmd' in self.settings:
            cmds.extend( self.settings['mpicmd'].split() )

        # binarypath
        cmds.extend( [self.lightaimdpath] )

        #Write input corrds
        with open(self.filename+'.xyz', 'w') as fin:
            fin.write(f'{len(current_coords)}\n')
            fin.write(f'charge {charge} multiplicity {mult} unit angstrom\n')
            for el, c in zip(qm_elems, current_coords):
                fin.write(f'{el} {c[0]} {c[1]} {c[2]}\n')
        cmds.extend(['--mol', f'{self.filename}.xyz'])

        # write MM charges as pointcharges if PC=True
        if PC == True:
            with open(self.filename+'.pc', 'w') as fin:
                fin.write(f'{len(MMcharges)}\n')
                # Mmcoords in Angstrom
                for mmcharge, mmcoord in zip(MMcharges, current_MM_coords):
                    fin.write(f'{mmcharge} {mmcoord[0]} {mmcoord[1]} {mmcoord[2]}\n')
            cmds.extend(['--mm-file', f'{self.filename}.pc'])

        for k, v in self.settings.items():
            if k == 'mpicmd':
                continue
            if v is not None:
                cmds.extend([f'--{k}', f'{v}'])
            else:
                cmds.extend([f'--{k}'])

        #Running lightaimd
        with open(self.filename + '.out', 'w') as ofile:
            print(" ".join(cmds))
            p = sp.run(cmds, check=True, stdout=ofile, stderr=ofile, universal_newlines=True)
        
        #Grab energy and possibly gradient
        self.energy, self.gradient, self.pcgradient = grabLightAIMDEandG(
            self.filename+'.out', len(qm_elems), Grad, len(MMcharges), PC)

        print(BC.OKBLUE, BC.BOLD, "------------ENDING LIGHTAIMD INTERFACE-------------", BC.END)

        print("Single-point energy:", self.energy)
        print_time_rel(module_init_time, modulename='LightAIMD run', moduleindex=2)

        if Grad and PC:
            return self.energy, self.gradient, self.pcgradient
        elif Grad:
            return self.energy, self.gradient
        else:
            return self.energy


#grab energy from output
def grabLightAIMDEandG(outfile, numatoms, Grad, numcharges, PC):
    energy, gradient, pcgradient = None, [], []
    energygrab, gradgrab, pcgradgrab = False, False, False
    with open(outfile, 'r') as fin:
        for line in fin:
            if '] SCF converged: steps' in line:
                energygrab = True
            elif energygrab and line.startswith('[20'):
                energygrab = False
            elif energygrab:
                if line.startswith('Total Energy'):
                    energy = float(line.split()[-1])

            if Grad and '] Total gradient:' in line:
                gradgrab = True
            elif gradgrab and line.startswith('[20'):
                gradgrab = False
            elif gradgrab:
                cols = line.split()
                if len(cols) == 5:
                    gradient.append([float(_) for _ in cols[-3:]])

            if Grad and PC and '] Total external gradient:' in line:
                pcgradgrab = True
            elif pcgradgrab and line.startswith('[20'):
                pcgradgrab = False
            elif pcgradgrab:
                cols = line.split()
                if len(cols) == 4:
                    pcgradient.append([float(_) for _ in cols[-3:]])
    if energy == None:
        print("Found no energy in LightAIMD outputfile:", outfile)
        ashexit()
    if Grad and len(gradient) != numatoms:
        print("Wrong number of gradient in file", outfile)
        ashexit()
    if PC and len(pcgradient) != numcharges:
        print("Wrong number of pcgradient in file", outfile)
        ashexit()
    return energy, np.asarray(gradient), np.asarray(pcgradient)
