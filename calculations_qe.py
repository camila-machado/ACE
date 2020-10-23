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

This module provides class interfaces to calculate the superconductivity 
critical temperature using quantum espresso programs and qeporgrams module.
"""
import numpy as np
import pickledb as pk
import os

from classtools import AttrDisplay
import models_qe as qe
import models_program as prg
import tools

class TcGridsWrapper(AttrDisplay):
    db_dict = 'grid'

    def __init__(self, cell , kdistance = 0, qdistance = 0):
        self.cell = cell
        self.kdistance = kdistance
        self.qdistance = qdistance
        self.coarse = qe.Grid('coarse')
        self.dense = qe.Grid('dense')
        self.qpoints = qe.Grid('qpoints')

    def calculate(self):

        print('\n-> Grids')

        if self.calc:
            self.reciprocal  = self._calcReciprocalModule(cell= self.cell)
            self.minimalgrid = self._calcMinimalGrid(reciprocal= self.reciprocal)
            self.c_factor = self._calcGridFactor(distance= self.kdistance)
            self.q_factor = self._calcGridFactor(distance= self.qdistance)
            self.d_factor = self._calcDenseFactor(c_factor= self.c_factor,
                                                q_factor= self.q_factor)
            
            self._setGrids()

            print('\nCalculated from kdistance and qdistance')

        else:

            print('\nUsed input grids')

    def _calcReciprocalModule(self, cell):
        a = np.zeros((3,3)) #direct vectors
        b = np.zeros((3,3)) #reciprocal vectors
        nb = np.zeros(3)    #reciprocal vectors' modules

        structure = cell.read()
        a = structure.cell
        
        pi = np.pi
        b[0] = 2*pi*np.cross(a[1],a[2])/np.dot(a[0],np.cross(a[1],a[2]))
        b[1] = 2*pi*np.cross(a[2],a[0])/np.dot(a[1],np.cross(a[2],a[0]))
        b[2] = 2*pi*np.cross(a[0],a[1])/np.dot(a[2],np.cross(a[0],a[1]))
        
        for i in range (0,3):
            nb[i] = np.sqrt(np.dot(b[i],b[i]))

        return nb

    def _calcMinimalGrid(self, reciprocal):
        div = np.zeros(3, dtype= np.int8)
        for i in range (0,3):
            div[i] = int(reciprocal[i] + 0.5)

        return div

    def _calcGridFactor(self, distance):
        factor = int(1/distance + 0.5)
        if factor == 0:
            factor = 1

        return factor 

    def _calcDenseFactor(self, c_factor, q_factor):
        factor = int(np.lcm(c_factor, q_factor))
        if factor == c_factor:
            factor = 2*c_factor

        return factor 

    def _correctType(self, factor, min_grid):

        div = [0,0,0]
        vector = factor * min_grid
        for i in range(0,3):
            div[i] = int(vector[i])

        return div

    def _setGrids(self):

        c_div = self._correctType(self.c_factor, self.minimalgrid)
        q_div = self._correctType(self.q_factor, self.minimalgrid)
        d_div = self._correctType(self.d_factor, self.minimalgrid)

        self.coarse.div = tuple(c_div)
        self.qpoints.div = tuple(q_div)
        self.dense.div = tuple(d_div)

    def load(self, database):
        db = pk.load(database, False)
        self.kdistance = db.dget(self.db_dict, 'kdistance')
        self.qdistance = db.dget(self.db_dict, 'qdistance')
        self.calc = db.dget(self.db_dict, 'calc')
        self.qpoints.load(database = database)        
        self.coarse.load(database = database)
        self.dense.load(database = database)

    def dump(self, database):
        db = pk.load(database, False)
        db.dadd(self.db_dict,('kdistance', self.kdistance))
        db.dadd(self.db_dict,('qdistance', self.qdistance))
        db.dadd(self.db_dict,('calc', self.calc))
        self.coarse.dump(database = database)
        self.dense.dump(database = database)       
        self.qpoints.dump(database = database)


class ICalculation():

    db_strategy = tools.IMakeDatabase()
    db_dict = 'calc_stage'
    calc_stage = None
    variables = [['routine'], ['qe_programs'], ['pseudo_folder'], ['np', 'nk']]
    grid_variable = None
    programs = None
    sufix = None

    def __init__(self, structure):
        assert self.sufix, 'Input sufix must be defined!'
        assert self.programs, 'Programs must be defined!'
        assert self.calc_stage, 'Calculation satages must be defined!'
        assert self.grid_variable, 'Variables must be defined!'

        self.cell = structure
        self.prefix = structure.name
        self.input_name = self.prefix + '.' + self.sufix + '.in'

        self.variables = ICalculation.variables
        self.variables.append(self.grid_variable)

        for prog in self.programs:
            self.variables.append(prog.variables)

    def calculate(self):
        assert self.database, 'Database must be defined'
        assert self.directory, 'Directory must be defined'
   
        current_dir = os.getcwd()

        os.chdir(self.directory.work)
        self._setVariables()
        self._setPrograms()
        self._runPrograms()

        os.chdir(current_dir)

    def _setVariables(self):
        assert False, '_setVariables() must be defined!'
    
    def _setPrograms(self):
        assert False, '_setPrograms() must be defined!'

    def _runPrograms(self):

        print('\n-> Run programs\n')

        for program in self.programs:
            program.set_Dir(dir= self.directory)

            if not(self.calc_stage.get('w_'+program.name)):
                os.chdir(self.directory.input)
                program.write()
                self.calc_stage.update({'w_'+program.name:True})
                self.dump(database= self.database.path)

            if not (self.calc_stage.get('r_'+program.name)):
                os.chdir(self.directory.calc)   
                program.run(mpi=self.mpi)
                self.calc_stage.update({'r_'+program.name:True})
                self.dump(database= self.database.path)
            else:
                print('Check: {} fineshed'.format(program.name))
        
    def makeDatabase(self, name):
        database = self.db_strategy.make_database(name = name)
        return database

    def load(self, database):   
        db = pk.load(database, False)
        self.calc_stage = db.dgetall(self.db_dict)

    def dump(self, database)
        db = pk.load(database, False)
        db.set(self.db_dict, self.calc_stage)
        db.dump()       

class TC(AttrDisplay, ICalculation):

    db_strategy = tools.TcDatabase()

    sufix = 'TC'

    calc_stage = {'calc_grids': False, 
     'w_scf1'  : False,  'w_scf2' : False, 
     'w_ph'    : False, 'w_q2r'   : False, 
     'w_matdyn': False, 'w_lambda': False, 
     'r_scf1'  : False,  'r_scf2' : False, 
     'r_ph'    : False, 'r_q2r'   : False, 
     'r_matdyn': False, 'r_lambda': False}

    variables = []

    grid_variable = ['coarse_div', 'coarse_off', 'qpoints_div', 'qpoints_off', 'dense_div', 
                     'dense_off', 'kdistance', 'qdistance', 'calc']

    pw1    = qe.Pwscf1()
    pw2    = qe.Pwscf()
    ph     = qe.Phonon()
    q2r    = qe.Q2r()
    matdyn = qe.Matdyn()
    lamb   = qe.Lambda()
    programs = [pw2, ph, q2r, matdyn, lamb]         
    
    def __init__(self, structure):
        ICalculation.__init__(self, structure = structure)
                              
    def _setVariables(self):
        
        self.pseudo = qe.Pseudo()   
        self.pseudo.load(database= self.database.path)

        self.grids = TcGridsWrapper(cell= self.cell)
        self.grids.load(database= self.database.path)

        if not(self.calc_stage.get('calc_grids')):
            self.grids.calculate()
            self.grids.dump(database= self.database.path)
            self.calc_stage.update({'calc_grids':True})
            self.dump(database= self.database.path)
        
        else:
            print('Check: Grids fineshed')

    def _setPrograms(self): 
      
        self.mpi    = qe.Mpi()
        self.pw1    = qe.Pwscf1(prefix= self.prefix, grid= self.grids.dense, 
                            cell= self.cell, pseudo= self.pseudo, name= 'scf1')
        self.pw2    = qe.Pwscf(prefix= self.prefix, grid= self.grids.coarse, 
                            cell= self.cell, pseudo= self.pseudo, name = 'scf2')
        self.ph     = qe.Phonon(prefix= self.prefix, grid= self.grids.qpoints)
        self.q2r    = qe.Q2r(prefix= self.prefix)
        self.matdyn = qe.Matdyn(prefix= self.prefix)
        self.lamb   = qe.Lambda(prefix= self.prefix, dyndir = self.directory.calc)

        for obj in [self.mpi, self.pw1, self.pw2, self.ph, self.q2r, self.matdyn, self.lamb]:
            obj.load(database= self.database.path)

        self.programs = [self.pw1, self.pw2, self.ph, self.q2r, self.matdyn, self.lamb]


class Energy(AttrDisplay, ICalculation):

    db_strategy = tools.EnergyDatabase()

    sufix = 'EN'

    calc_stage = {'w_scf'  : False, 
                  'r_scf'  : False}

    grid_variable = ['kpoints_div', 'kpoints_off']

    pw = qe.Pwscf()
    programs = [pw]

    def __init__(self, structure):
        ICalculation.__init__(self, structure = structure)
                              
    def _setVariables(self):
        
        self.pseudo = qe.Pseudo()   
        self.pseudo.load(database= self.database.path)

        self.grid = qe.Grid(name= 'kpoints')
        self.grid.load(database= self.database.path)

    def _setPrograms(self):  
        self.mpi    = qe.Mpi()
        self.pw     = qe.Pwscf(prefix= self.prefix, grid= self.grid, 
                            cell= self.cell, pseudo= self.pseudo)

        for obj in [self.mpi, self.pw]:
            obj.load(database= self.database.path)

        self.programs = [self.pw]


class PhononGamma(AttrDisplay, ICalculation):

    db_strategy = tools.PhononDatabase()

    calc_stage = {'w_scf': False, 'w_ph': False, 
                  'r_scf': False, 'r_ph': False}
    sufix = 'PH'
    
    grid_variable = ['kpoints_div', 'coarse_off', 'kpoints_off']

    pw    = qe.Pwscf1()
    ph    = qe.Phonon()
    programs = [pw, ph] 
    
    def __init__(self, structure):
        ICalculation.__init__(self, structure = structure)
        
    def _setVariables(self):
        
        self.pseudo = qe.Pseudo()
        self.pseudo.load(database= self.database.path)

        self.grid = qe.Grid(name= 'kpoints')
        self.grid.load(database= self.database.path)

    def _setPrograms(self): 
      
        self.mpi    = qe.Mpi()
        self.pw     = qe.Pwscf(prefix= self.prefix, grid= self.grid,
                            cell= self.cell, pseudo= self.pseudo)
        self.ph     = qe.Phonon(prefix= self.prefix, vector = ' 0.0 0.0 0.0')

        for obj in [self.mpi, self.pw, self.ph]:
            obj.load(database= self.database.path)

        self.programs =  [self.pw, self.ph]

################################################################################
##----------------------------------------------------------------------------##
################################################################################

if __name__ == '__main__':

    #Unit Tests
    cell_name = 'Nb.cif'
    prefix = cell_name.split('.')[0]
    cell_atoms = qe.CellStructure(filename= cell_name)

    #Test class TC
    print('\nTEST TC\n')
    tc_prefix = prefix + '_TC'
    tc_dbname = tc_prefix + '.db'
    
    print('\nTEST 1: criate a TC obj from .cif file')
    tc_test = TC(structure= cell_atoms)
    print(tc_test)

    print('\nTEST 2: criate a TC general database')
    tc_test.database = tc_test.makeDatabase(name= tc_dbname)
    print(tc_test.database)

    print('\nTEST 3: criate a TC directory')
    tc_directory = prg.Directory(workpath= os.path.join(os.getcwd(), tc_prefix))
    tc_directory.makedirs()
    tc_test.directory = tc_directory
    print(tc_test.directory)

    print('\nTEST 4: calculate TC')
    tc_test.calculate()
    
    #Test class PHONON_GAMMA
    print('\nTEST PHONON_GAMMA\n')
    ph_prefix = prefix + '_PH'
    ph_dbname = ph_prefix + '.db'
    
    print('\nTEST 1: criate a PHONON_GAMMA obj from .cif file')
    ph_test = PhononGamma(structure= cell_atoms)
    print(ph_test)

    print('\nTEST 2: criate a PHONON_GAMMA general database')
    ph_test.database = tc_test.makeDatabase(name= ph_dbname)
    print(ph_test.database)

    print('\nTEST 3: criate a PHONON_GAMMA directory')
    ph_directory = prg.Directory(workpath= os.path.join(os.getcwd(), ph_prefix))
    ph_directory.makedirs()
    ph_test.directory = ph_directory
    print(ph_test.directory)

    print('\nTEST 4: calculate PHONON_GAMMA')
    ph_test.calculate()

    #Test class ENERGY
    print('\nTEST ENERGY\n')
    en_prefix = prefix + '_EN'
    en_dbname = en_prefix + '.db'
    
    print('\nTEST 1: criate a ENERGY obj from .cif file')
    en_test = PhononGamma(structure= cell_atoms)
    print(en_test)

    print('\nTEST 2: criate a ENERGY general database')
    en_test.database = en_test.makeDatabase(name= en_dbname)
    print(en_test.database)

    print('\nTEST 3: criate a ENERGY directory')
    en_directory = prg.Directory(workpath= os.path.join(os.getcwd(), en_prefix))
    en_directory.makedirs()
    en_test.directory = en_directory
    print(en_test.directory)

    print('\nTEST 4: calculate ENERGY')
    en_test.calculate()
    