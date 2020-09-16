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
from models_qe import*
from models_program import*
from tools import WriteStrategy

class ICalculation():

    db_dict = 'calc_stage'
    calc_stage = None
    variables = None
    input_sufix = None

    def __init__(self, strtucture):
        assert self.calc_stage, 'Calculation satages must be defined!'
        assert self.variables, 'Variables must be defined!'
        assert self.variables, 'Input sufix must be defined!'

        self.cell = strtucture
        self.prefix = strtucture.name
        self.input_name = self.prefix + self.input_sufix

    def calculate(self):
        assert self.database, 'Database must be defined'
        assert self.database, 'Directory must be defined'

        os.chdir(self.dir.work)
        self._setVariables()
        self._setPrograms()
        self._runPrograms()

    def _setVariables(self):
        assert False, '_setVariables() must be defined!'
    
    def _setPrograms(self):
        assert False, '_setPrograms() must be defined!'

    def _runPrograms(self):
        assert False, '_runPrograms() must be defined!'

    def makeDatabase(self):
        assert False, 'generateDatabase() must be defined!'

    def load(self, database):
        db = pk.load(database, False)
        self.calc_stage = db.dgetall(self.db_dict)
    
    def dump(self, database):
        db = pk.load(database, False)
        db.set(self.db_dict, self.calc_stage)
        db.dump()       

class TcGridsWrapper(AttrDisplay):
    db_dict = 'grid'

    def __init__(self, cell , kdistance = 0, qdistance = 0):
        self.cell = cell
        self.kdistance = kdistance
        self.qdistance = qdistance
        self.coarse = Grid('coarse')
        self.dense = Grid('dense')
        self.qpoints = Grid('qpoints')

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

        a = cell.structure.cell
        
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

        self.coarse.setGrid(div = c_div)
        self.qpoints.setGrid(div = q_div)
        self.dense.setGrid(div = d_div)

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

class CalcTC(AttrDisplay, ICalculation):

    calc_stage = {'calc_grids': False, 
     'w_scf1'  : False,  'w_scf2'  : False, 
     'w_ph'    : False, 'w_q2r'   : False, 
     'w_matdyn': False, 'w_lambda': False, 
     'r_scf1'  : False,  'r_scf2'  : False, 
     'r_ph'    : False, 'r_q2r'   : False, 
     'r_matdyn': False, 'r_lambda': False}

    variables = [['routine'], ['qe_programs'], ['pseudo_folder'], ['np', 'nk'], ['coarse_div', 'coarse_off', 
              'qpoints_div', 'qpoints_off', 'dense_div', 'dense_off', 'kdistance', 'qdistance',
              'calc'], ['restart_mode', 'occupations', 'smearing', 'degauss', 'conv_thr'], 
              ['tr2_ph', 'ldisp', 'recover', 'electron_phonon', 'el_ph_sigma', 'el_ph_nsigma'],
              ['zasr', 'la2F'], ['asr', 'la2F', 'dos', 'nk1', 'nk2', 'nk3', 'ndos'], 
              ['sigma_omega', 'mu']]
    
    input_sufix = '.TC.in'

    def __init__(self, strtucture):
        ICalculation.__init__(self, strtucture = strtucture)
                              
    def _setVariables(self):
        
        self.pseudo = Pseudo()   
        self.pseudo.load(database= self.database.file)

        self.grids = TcGridsWrapper(cell= self.cell)
        self.grids.load(database= self.database.file)

        if not(self.calc_stage.get('calc_grids')):
            self.grids.calculate()
            self.grids.dump(database= self.database.file)
            self.calc_stage.update({'calc_grids':True})
            self.dump(database= self.database.file)
        
        else:
            print('Check: Grids fineshed')

    def _setPrograms(self): 
      
        self.mpi    = Mpi()
        self.pw1    = Pwscf1(prefix= self.prefix, grid= self.grids.dense, 
                            cell= self.cell, pseudo= self.pseudo)
        self.pw2    = Pwscf2(prefix= self.prefix, grid= self.grids.coarse, 
                            cell= self.cell, pseudo= self.pseudo)
        self.ph     = Phonon(prefix= self.prefix, grid= self.grids.qpoints)
        self.q2r    = Q2r(prefix= self.prefix)
        self.matdyn = Matdyn(prefix= self.prefix)
        self.lamb   = Lambda(prefix= self.prefix)

        for obj in [self.mpi, self.pw1, self.pw2, self.ph, self.q2r, self.matdyn, self.lamb]:
            obj.load(database= self.database.file)

    def _runPrograms(self):

        print('\n-> Run programs\n')

        for program in [self.pw1, self.pw2, self.ph, self.q2r, self.matdyn, self.lamb]:
            program.set_Dir(dir= self.dir)

            if not(self.calc_stage.get('w_'+program.name)):
                os.chdir(self.dir.input)
                program.write()
                self.calc_stage.update({'w_'+program.name:True})
                self.dump(database= self.database.file)

            if not (self.calc_stage.get('r_'+program.name)):
                os.chdir(self.dir.calc)   
                program.run(mpi=self.mpi)
                self.calc_stage.update({'r_'+program.name:True})
                self.dump(database= self.database.file)
            else:
                print('Check: {} fineshed'.format(program.name))

    def makeDatabase(self, name):
        
        _DATABASE = name

        db = pk.load(_DATABASE, False)

        #optional: in case quantum-espresso programs are not in linux $PATH
        #qe - quantum espresso programs' path
        db.dcreate ('calculation')
        db.dadd('calculation',('routine','TC') )
        
        #qe - quantum espresso programs' path
        db.dcreate ('qe')
        db.dadd('qe',('qe_programs','') )

        #pseudo - pseudo potentials' folder path
        db.dcreate('pseudo')
        db.dadd('pseudo',('pseudo_folder','/home/camila/Documentos/EMA/Program-TC/pseudo/USPP') )

        #mpi - parallen running description
        db.dcreate ('mpi')

        db.dadd('mpi',('np',4) )
        db.dadd('mpi',('nk',4) )

        #grids - points distance used to calculate grids for scf and phonons calculations
        db.dcreate ('grid')

        db.dadd('grid',('coarse_div',(9,9,9)) )
        db.dadd('grid',('coarse_off',(0,0,0)) )
        db.dadd('grid',('qpoints_div',(3,3,3)) )
        db.dadd('grid',('qpoints_off',(0,0,0)) )
        db.dadd('grid',('dense_div',(18,18,18)) )
        db.dadd('grid',('dense_off',(0,0,0)) )
        db.dadd('grid',('kdistance', 0.3) ) #distance in angstrons of two consecutive grid k-points
        db.dadd('grid',('qdistance', 0.9) ) #distance in angstrons of two consecutive grid q-points
        db.dadd('grid',('calc', 1) ) #habilitate the calculation from reciprocal vector

        #Quantum-espresso programs' parameters:
        db.dcreate ('pw_par')

        db.dadd('pw_par',('restart_mode','from_scratch') )
        db.dadd('pw_par',('occupations','smearing') )
        db.dadd('pw_par',('smearing','marzari-vanderbilt') )
        db.dadd('pw_par',('degauss',0.05) )
        db.dadd('pw_par',('conv_thr',1e-10) )

        db.dcreate ('ph_par')

        db.dadd('ph_par',('tr2_ph',1e-12) )
        db.dadd('ph_par',('ldisp','.true.') )
        db.dadd('ph_par',('recover','.false.') )
        db.dadd('ph_par',('electron_phonon','interpolated') )
        db.dadd('ph_par',('el_ph_sigma',0.005) )
        db.dadd('ph_par',('el_ph_nsigma',10) )

        db.dcreate ('q2r_par')

        db.dadd('q2r_par',('zasr','simple') )
        db.dadd('q2r_par',('la2F','.true.') )

        db.dcreate ('matdyn_par')

        db.dadd('matdyn_par',('asr','simple') )
        db.dadd('matdyn_par',('la2F','.true.') )
        db.dadd('matdyn_par',('dos','.true.') )
        db.dadd('matdyn_par',('nk1',10) )
        db.dadd('matdyn_par',('nk2',10) )
        db.dadd('matdyn_par',('nk3',10) )
        db.dadd('matdyn_par',('ndos',50) )

        db.dcreate ('lambda_par')

        db.dadd('lambda_par',('sigma_omega', 0.12 ) )
        db.dadd('lambda_par',('mu',0.16) )

        return db

class CalcEnergy(AttrDisplay, ICalculation):
    pass

class CalcPhonon(AttrDisplay, ICalculation):
    pass

################################################################################
##----------------------------------------------------------------------------##
################################################################################

if __name__ == '__main__':
    pass