import time

from ash.functions.functions_general import ashexit, BC,print_time_rel, print_line_with_mainheader,listdiff
import ash.modules.module_coords
from ash.modules.module_results import ASH_Results
import numpy as np
import os
import shutil

##########################
#MLatom Theory interface
##########################
#NOTE: we only intend to support ML functionality in MLatom (not basic QM methods)

#MLatom has 3 types of models:
#1) methods: these are pre-trained or at least standalone electronic structure methods
#2) ml_model: these are models that require training
#3) model_tree_node: some composite models (unclear)

#methods examples: AIQM1, AIQM1@DFT, AIQM1@DFT*, ANI-1ccx, ANI-1x, ANI-1x-D4, ANI-2x, ANI-2x-D4
#ml_model examples: ani, dpmd, gap, physnet, sgdml

#MLatom may use/require interfaces to : TorchANI, DeepMD-kit, GAP/QUIP, Physnet, sGDML

#NOTES on interface:
# AIQMx methods require either MNDO or Sparrow as QM-program. Also dftd4.
   #Sparrow lacks AIQMx gradient (for ODM2 part), only energies available.



class MLatomTheory:
    def __init__(self, printsetting=False, printlevel=2, numcores=1, label="mlatom", filename="sdf",
                 method=None, ml_model=None, qm_program=None):
        module_init_time=time.time()
        self.theorynamelabel="MLatom"
        self.theorytype="QM"
        self.analytic_hessian=False

        print_line_with_mainheader(f"{self.theorynamelabel} initialization")

        try:
            #
            import mlatom as ml
        except ModuleNotFoundError:
            print("MLatom  requires installation of mlatom")
            print("See: http://mlatom.com/docs/installation.html")
            print("Try: pip install mlatom")
            print("You probably also have to do: pip install scipy torch torchani tqdm matplotlib statsmodels h5py pyh5md")
            ashexit()

        #EARLY EXITS

        if method is None and ml_model is None:
            print("Neither a method or ml_model was selected for MLatomTheory interface. Exiting.")
            ashexit()

        #Store optional properties of run in dict
        self.properties ={}

        #Printlevel
        self.printlevel=printlevel
        self.label=label
        self.filename=filename


        #METHODS: pre-trained models
        #Note: useful method keywords in MLAatom below 
        # AIQMx models require either MNDO or Sparrow. Also dftd4 (except AIQM1@DFT* but that one is bad anyway)
        #'AIQM1', 'AIQM1@DFT', 'AIQM1@DFT*', 
        #ANI models: will download parameters automatically
        #'ANI-1ccx', 'ANI-1x', 'ANI-1x-D4', 'ANI-2x', 'ANI-2x-D4'
        self.method=method
        self.qm_program=qm_program
        print("Checking if method or ml_model was selected")
        print("Method:", self.method)
        print("QM program:", self.qm_program)

        if self.method is not None:
            if 'AIQM' in self.method :    
                print("An AIQMx method was selected")
                print("Warning: this requires setting qm_program keyword as either mndo or sparrow.")
                print("Also dftd4 D4-dispersion program")
                if self.qm_program == 'mndo':
                    print("QM program is mndo")
                    print("Make sure executable mndo is in your environment")
                    print("See https://mndo.kofo.mpg.de about MNDO licenses")
                    try:
                        mndodir = os.path.dirname(shutil.which('mndo2020'))
                        os.environ['mndobin'] = mndodir+"/mndo2020"
                        print("Found mndo2020 executable in:", mndodir)
                    except TypeError:
                        print("Found no mndo2020 executable in your environment. Exiting.")
                        ashexit()
                    try:
                        dftd4dir = os.path.dirname(shutil.which('dftd4'))
                        os.environ['dftd4bin'] = dftd4dir+"/dftd4"
                        print("Found dftd4 executable in:", dftd4dir)
                    except TypeError:
                        print("Found no dftd4 executable in your environment. Exiting.")
                        ashexit()
                elif self.qm_program == 'sparrow':
                    print("QM program is sparrow")
                    print("Make sure executable sparrow is in your environment")
                    print("See https://github.com/qcscine/sparrow. Possible installation  via: conda install scine-sparrow-python")
                    print("Also make sure dftd4 (https://github.com/dftd4/dftd4) is in your environment. Possible installation  via:  conda install dftd4")
                    print("Warning: sparrow lacks AIQMx gradient (for ODM2 part), only energies available.")
                    try:
                        sparrowdir = os.path.dirname(shutil.which('sparrow'))
                        os.environ['sparrowbin'] = sparrowdir+"/sparrow"
                        print("Found sparrow executable in:", sparrowdir)
                    except TypeError:
                        print("Found no sparrow executable in your environment. Exiting.")
                        ashexit()
                    try:
                        dftd4dir = os.path.dirname(shutil.which('dftd4'))
                        os.environ['dftd4bin'] = dftd4dir+"/dftd4"
                        print("Found dftd4 executable in:", dftd4dir)
                    except TypeError:
                        print("Found no dftd4 executable in your environment. Exiting.")
                        ashexit()
                    
                else:
                    print("QM program keyword is neither mndo or sparrow. Not allowed, exiting.")
                    ashexit()
            elif 'ANI' in self.method:
                print("An ANI type method was selected")
                print("This requires TorchANI and pytorch")
                print("Note: ANI parameters will be downloaded automatically if needed")
            elif 'ODM' in self.method:
                print("A ODMx type semi-empirical method was selected. This requires MNDO")
            elif 'OM' in self.method:
                print("A OMx type semi-empirical method was selected. This requires MNDO")
            else:
                print(f"Either an invalid  {self.method} or unknown method (to MLatomTheory interface) was selected. Exiting.")
                ashexit()

        print_time_rel(module_init_time, modulename='MLatom creation', moduleindex=2)

    #General run function
    def run(self, current_coords=None, current_MM_coords=None, MMcharges=None, qm_elems=None, mm_elems=None,
            elems=None, Grad=False, PC=False, numcores=None, restart=False, label=None,
            charge=None, mult=None):
        module_init_time=time.time()
        import mlatom as ml

        print(BC.OKBLUE,BC.BOLD, f"------------RUNNING {self.theorynamelabel} INTERFACE-------------", BC.END)

        #Prepare for run
        #mlatom.data
        molecule = ml.data.molecule(charge, mult)
        molecule.read_from_numpy(coordinates=current_coords, species=np.array(elems))

        #mlatom.models
        #Comp chem models, 3 types: methods (used as is), ml_model (requires training), model_tree_node (composite)
        if self.method != None:
            print("A method was selected: ", self.method)
            print("QM program:", self.qm_program)
            print("Creating model")
            model = ml.models.methods(method=self.method, qm_program=self.qm_program) 
        else:
            print("No method was defined yet.")
            ashexit()        

        #Create dftd4.json file before running if required
        if 'AIQM' in self.method:
            print("An AIQMx method was selected")
            #NOTE: dftd4 interface of MLatom has a bug for current dftd4 release
            if 'AIQM1@DFT*' in self.method:
                print("AIQM1@DFT* method was selected, no disp. correction needed")
            elif 'AIQM1@DFT' in self.method:
                print("AIQM1@DFT method was selected. Disp. correction needed")
            elif 'AIQM1' in self.method:
                print("AIQM1 method was selected. Disp. correction needed")

        #Run
        if PC is True:
            print("PC is not yet supported by MLatomTheory interface")
            #Note: MNDO should support PCs, not sure about sparrow
            ashexit()
        else:
            if Grad is True:
                model.predict(molecule=molecule,calculate_energy=True,
                            calculate_energy_gradients=True,
                            calculate_hessian=False)
                self.energy = molecule.energy
                self.gradient = molecule.get_energy_gradients()
                print("Single-point MLatom energy:", self.energy)

                print_time_rel(module_init_time, modulename='MLatom run', moduleindex=2)
                return self.energy,self.gradient

            else:
                model.predict(molecule=molecule, calculate_energy=True)
                self.energy = molecule.energy
                print("Single-point MLatom energy:", self.energy)
                #std = molecule.aiqm1_nn.energy_standard_deviation
                #print("std:", std)

                print_time_rel(module_init_time, modulename='MLatom run', moduleindex=2)
                return self.energy

