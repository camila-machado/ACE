"""
##------------------------------------------------------------------------------
# CNPEM - National Brazilian Center for Research in Energy and Materials
# Sirius - Research Group EMA
#
# Code project: TC_caculus
# 
# Objectif:
# Automate the use of Quantum Espresso (QE) for the calculation of 
# Superconductivity Critical Temperature (Tc) of diferent molecules and 
# cell structures
# 
##------------------------------------------------------------------------------

This module provide class interfaces for quantum espresso programs
"""
import os
import subprocess
import sys
import pickledb as pk

import tools 
from classtools import AttrDisplay
from models_program import CellStructure

class Pseudo(AttrDisplay):
    db_dict = 'pseudo'

    def __init__(self, folder = None):
        
        if folder:
            folder = self._checkfolder(path= folder)
            self.folder = folder
            self._setFiles()
    
    def _checkfolder(self, path):
        folder = os.path.join(os.getcwd(), path)
        folder = os.path.normpath(folder)
        assert os.path.isdir(folder), 'This is not an existing directory'
     
        return folder

    def _setFiles(self):     
        pseudo_files = os.listdir(self.folder)
        self.files = {}
        for file_name in pseudo_files:
            element = file_name.split('.')[0]
            self.files.update({element:file_name})
    
    def selectFiles (self, elements):
        pseudo_dict = {}
        elements = list( dict.fromkeys(elements) )#remove repeated
        for element in elements:
            pseudo_file = self.files.get(element)
            pseudo_dict.update({element:pseudo_file})
        return pseudo_dict

    def load (self, database):
        db = pk.load(database, False)
        folder =  db.dget(self.db_dict,'pseudo_folder')
        folder = self._checkfolder(path= folder)
        self.folder = folder
        self._setFiles()
        db.dump()
    
    def dump(self, database):
        db = pk.load(database, False)
        db.dadd(self.db_dict,('pseudo_folder', self.folder))
        db.dump()

class Grid(AttrDisplay):
    db_dict = 'grid'

    def __init__(self, name, div = (0,0,0), off = (0,0,0)):
        self.name = name
        self.div = div
        self.off = off

    def load(self, database):
        db = pk.load( database, False)
        self.div = tuple(db.dget(self.db_dict, self.name+'_div')) 
        self.off = tuple(db.dget(self.db_dict, self.name+'_off'))

    def dump(self, database):
        db = pk.load(database, False)
        db.dadd(self.db_dict,(self.name+'_div', self.div))
        db.dadd(self.db_dict,(self.name+'_off', self.off))
        db.dump()

class Mpi(AttrDisplay):
    """

    This class gather information to run quantum espresso programs  
    with parallel version
    """
    db_dict = 'mpi'

    def __init__(self, np = 1, nk = 1): 
        self.np = str(np)
        self.nk = str(nk)

    def load(self, database):
        db = pk.load(database, False)
        self.np = str(db.dget(self.db_dict,'np'))  
        self.nk = str(db.dget(self.db_dict,'nk'))
    
    def dump(self, database):
        db = pk.load(database, False)
        db.dadd(self.db_dict,('np', self.np) )
        db.dadd(self.db_dict,('nk', self.nk) )
        db.dump()        

class IQuantumEspresso(AttrDisplay):

    write_strategy = tools.IWriteStrategy()
    program = None
    db_dct = None
    variables = None
    section = 'input'

    def __init__(self, **keyargs):  
        assert self.variables, 'File parameters order must be defined!'
        assert self.program, 'Program name must be defined!'
        assert self.db_dict, 'Database dict string name must be defined!'

        if len(keyargs) > 0:
            self.parameters(**keyargs)

    def parameters(self, prefix, name = None):

        if name:
            self.name = name
        else:
            self.name = self.program.replace('.x','') 

        self.input =  prefix + '.' + self.name + '.in'
        self.output =  prefix + '.' + self.name + '.out'   
        self.parameters = {}
        self._prefixParameters(prefix)
    
    def _prefixParameters(self, prefix):
        assert False, '"setParameters must be defined!"'

    def set_Dir (self, dir):
        self.input = os.path.join(dir.input, self.input)
        self.output = os.path.join(dir.output, self.output)
        self.calc = dir.calc

    def _command(self, mpi):
        command = 'mpiexec'+' -np '+mpi.np+' '+self.program+' -nk '+ mpi.nk+' -in '+self.input
        return command
    
    def load(self, database):
        db = pk.load(database, False)

        par = dict(db.dgetall(self.db_dict).items())
        keys = par.keys()
        
        for key in keys:
            self.parameters.update({key: par.get(key)})
    
    def run(self, mpi):
        command = self._command(mpi= mpi)

        text = '\nrunning {program_txt}...'
        print(text.format(program_txt = self.program))
        try:
            process = subprocess.run(command, shell = True,
                                            stdout = subprocess.PIPE, 
                                            universal_newlines = True, 
                                            stderr=subprocess.PIPE,
                                            check = True)
        except subprocess.CalledProcessError as err:
            print("Execution failed:", err, file=sys.stderr)
            raise
        except OSError as err:
            print("Execution failed:", err, file=sys.stderr)
            raise
        else:
            with open(self.output, 'w') as f:
                f.write(str(process.stdout))
            text = '{program_txt} fineshed! :D'
            print(text.format(program_txt = self.name))

    def write(self):
        self.write_strategy.write(section = self.section, 
                                 filename = self.input, 
                                 file_order = self.variables, 
                                 parameters = self.parameters)

class Pwscf (IQuantumEspresso):
    write_strategy = tools.WriteScf()
    program = 'pw.x'
    db_dict = 'pw_par'
    variables = ['prefix','restart_mode','pseudo_dir','outdir','occupations',
                 'smearing','degauss','ecutwfc','ecutrho','conv_thr', 'la2F']
    
    def parameters(self, prefix, pseudo, cell, grid, name = None):
        IQuantumEspresso.parameters(self, prefix= prefix, name= name)
        self.grid = grid
        self.cell = cell
        if cell:
            self.pseudo = pseudo.selectFiles(cell.elements())
            self.parameters.update({'pseudo_dir': pseudo.folder})
            self._setEcut()

    def _prefixParameters(self, prefix):
        self.parameters.update({'prefix': prefix})
        self.parameters.update({'outdir': './out'})

    def _setEcut(self):
        pseudo_files =[]
        for name in self.pseudo.values():
            file_path = os.path.join(self.parameters.get('pseudo_dir'), name)
            with open(file_path, 'r') as f:
                file_content = f.readlines()
            pseudo_files.append(file_content)
        
        ecutwfc_values = []
        ecutrho_values = []
        for content in pseudo_files:
            for line in content:
                if line.find('Suggested minimum cutoff for wavefunctions') > -1:
                    aux = line.split(':')[1]
                    aux2 = aux.strip(' .Ry\n')
                    ecutwfc_values.append(int(aux2))
                elif line.find('Suggested minimum cutoff for charge density') > -1:
                    aux = line.split(':')[1]
                    aux2 = aux.strip(' .Ry\n')
                    ecutrho_values.append(int(aux2))
                    
        ecutwfc = max(ecutwfc_values)
        ecutrho = max(ecutrho_values)
        
        self.parameters.update({'ecutwfc': ecutwfc})
        self.parameters.update({'ecutrho': ecutrho})

    def write(self):
        self.write_strategy.write(filename = self.input,
                                images = self.cell.read(),
                                input_data = self.parameters,
                                pseudopotentials = self.pseudo,
                                kpts = self.grid.div,
                                koffset = self.grid.off)

class Pwscf1(Pwscf):

    def write(self):
        Pwscf.write(self)
        self.write_strategy.addParameter(key= 'la2F', value = '.true.', 
                                         section = 'system', filename = self.input)

class Phonon (IQuantumEspresso):

    write_strategy = tools.WriteSection()
    
    program = 'ph.x'
    db_dict = 'ph_par'

    variables = ['prefix', 'outdir', 'tr2_ph', 'fildyn','ldisp', 'nq1', 'nq2', 'nq3', 
                  'electron_phonon','fildvscf','recover', 'el_ph_sigma', 'el_ph_nsigma']
                
    section = 'inputph'

    def parameters(self, prefix, name = None, grid = None, vector= None):
        IQuantumEspresso.parameters(self, prefix= prefix, name = name)
        self.vector = vector
        if grid:
            self.parameters.update({'nq1': grid.div[0]})
            self.parameters.update({'nq2': grid.div[1]})
            self.parameters.update({'nq3': grid.div[2]})

    def _prefixParameters(self, prefix):
        self.parameters.update({'prefix': prefix})
        self.parameters.update({'outdir': './out'})
        self.parameters.update({'fildyn': prefix+'.dyn'})
        self.parameters.update({'fildvscf': prefix+'.dv'})

    def write(self):
        IQuantumEspresso.write(self)
        if self.vector:
            self.write_strategy.addLine(line= self.vector, 
                                       position= '/\n', 
                                       filename= self.input)

class Q2r (IQuantumEspresso):

    write_strategy = tools.WriteSection()

    program = 'q2r.x'
    db_dict = 'q2r_par'

    variables = ['zasr', 'fildyn', 'flfrc', 'la2F']

    def _prefixParameters(self, prefix):
        self.parameters.update({'fildyn': prefix+'.dyn'})
        self.parameters.update({'flfrc': prefix+'.frc'})
    
    def _command(self, mpi):
        command = self.program + ' < ' + self.input
        return command

class Matdyn (IQuantumEspresso):

    write_strategy = tools.WriteSection()

    program = 'matdyn.x'
    db_dict = 'matdyn_par'

    variables = ['asr', 'flfrc', 'flfrq', 'la2F', 'dos', 'fldos', 'nk1',
                  'nk2','nk3','ndos']

    def _prefixParameters(self, prefix):
        self.parameters.update({'flfrc': prefix+'.frc'})
        self.parameters.update({'flfrq': prefix+'.freq'})
        self.parameters.update({'fldos': prefix+'.phonon.dos'})

    def _command(self, mpi):
        command = 'mpiexec'+' -np '+ mpi.np +' '+ self.program +' -in '+ self.input
        return command

class Lambda (IQuantumEspresso):

    write_strategy = tools.WriteLambda()

    program = 'lambda.x'
    db_dict = 'lambda_par'

    variables = ['zasr', 'fildyn', 'flfrc', 'la2F']

    def parameters(self, prefix, name = None, dyndir= ''):
        IQuantumEspresso.parameters(self, prefix= prefix, name = name)
        self.dyndir = dyndir

    def _prefixParameters(self, prefix):
        self.parameters.update({'fildyn': prefix+'.dyn'})

    def _command(self, mpi):
        command = self.program + ' < ' + self.input
        return command

    def write(self):
        self.write_strategy.write(filename= self.input, 
                                  dyndir= self.dyndir,
                                  sigma_omega= self.parameters.get('sigma_omega'), 
                                  mu= self.parameters.get('mu'), 
                                  fildyn= self.parameters.get('fildyn'))    

################################################################################
##----------------------------------------------------------------------------##
################################################################################

if __name__ == '__main__':    
        
    #Unit Tests
    database = 'Nb_config.db'
    
    #Test class Pseudopotentials
    pseudo = Pseudo()
    print('\nTEST PSEUDO\n')
    print('\nTEST 1: criate empty pseudo obj\n', pseudo)

    pseudo = Pseudo(folder= 'pseudo')
    print('\nTEST 2: criate pseudo obj from folder of pseudopotential files:\n', pseudo)
    
    selection = pseudo.selectFiles(['Nb','H'])
    print('\nTEST 3: select pseudopotentials form list of elements\n', selection)

    pseudo.dump(database= database)
    pseudo2 = Pseudo()
    pseudo2.load(database= database)
    print('\nTEST 4: dump and load from database\n', pseudo2)
   
    #Test class Grid
    print('\nTEST GRID\n')

    print('\nTEST 1: criate grid obj and dump at database\n')

    kgrid = Grid('dense', div=(18,18,18))
    cgrid = Grid('coarse', div=(9,9,9))
    qgrid = Grid('qpoints', div=(3,3,3))

    for obj in [kgrid, cgrid, qgrid]:
        print(obj.name, ':', obj)
        obj.dump(database= database)

    print('\nTEST 2: criate empty grid obj and load from database\n')

    kgrid2 = Grid('dense')
    cgrid2 = Grid('coarse')
    qgrid2 = Grid('qpoints')

    for obj in [kgrid2, cgrid2, qgrid2]:
        print(obj.name, ':', obj)
        obj.load(database= database)
        print(obj.name, ':', obj)

    
    #Test class Mpi
    print('\nTEST MPI\n')

    print('\nTEST 1: criate mpi obj and dump at database\n')
 
    mpi = Mpi(nk= 3,np= 3)
    print(mpi)
    mpi.dump(database= database)

    print('\nTEST 2: criate empty mpi obj and load from database\n')

    mpi = Mpi()
    print(mpi)
    mpi.load(database= database)
    print(mpi)   
     
    #Test classes inheritated from IQuantumEspresso
    print('\nTEST IQE PROGRAMS\n')

    #load input information for porgrams
    database = 'Nb_config.db'
    prefix = 'Nb'
    cell = CellStructure('Nb.cif')

    mpi = Mpi()
    pseudo = Pseudo()
    grid1  = Grid('dense')
    grid2  = Grid('coarse')
    qgrid = Grid('qpoints')

    for obj in [pseudo, grid1, grid2, qgrid, mpi]:
        obj.load(database= database)

    print("\nTEST 1: criate program's objs and load from database\n")

    pw1    = Pwscf1(prefix= prefix, name= 'scf1', grid= grid1, cell= cell, pseudo= pseudo)
    pw2    = Pwscf(prefix= prefix, name= 'scf2', grid= grid2, cell= cell, pseudo= pseudo)
    ph     = Phonon(prefix= prefix, grid= qgrid)
    q2r    = Q2r(prefix= prefix)
    matdyn = Matdyn(prefix= prefix)
    lamb   = Lambda(prefix=prefix)

    for obj in [pw1, pw2, ph, q2r, matdyn, lamb]:
        obj.load(database= database)
        print(obj, '\n')
    
    print("\nTEST 2: write program's input files\n")
    for program in [pw1, pw2, ph, q2r, matdyn]:
        program.write()

    print("\nTEST 3: run program's input files\n")
    for program in [pw1, pw2, ph, q2r, matdyn]:
        program.run(mpi=mpi)

    print("\nTEST 4: write and run lambda input files\n")
    lamb.write()
    lamb.run(mpi=mpi)
    