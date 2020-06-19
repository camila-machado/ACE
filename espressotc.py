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
from qeprograms  import *

class Directory:
    db_dict = 'dir'

    def __init__(self, work_path):
        self.work = work_path
        self.calc = os.path.join(work_path, 'calc_dir')
        self.output = os.path.join(work_path, 'out_dir')
        self.input = os.path.join(work_path, 'input_dir')

    def makedirs(self):
        os.makedirs(self.calc)
        os.makedirs(self.output)
        os.makedirs(self.input)       
   
    def delete(self):
        shutil.rmtree(self.calc)
        shutil.rmtree(self.output)
        shutil.rmtree(self.input)

    def dump(self, database):
        pass

    def load(self, database):
        db = pk.load(database, False)
        self.qe = db.dget(db_dict, 'qe_programs')

class GridsWrapper(AttrDisplay):
    db_dict = 'grids'

    def __init__(self, cell , kdistance = 0, qdistance = 0):
        self.cell = cell
        self.kdistance = kdistance
        self.qdistance = qdistance
        self.coarse = Grid('coarse')
        self.dense = Grid('dense')
        self.qpoints = Grid('qpoints')

    def calculate(self):
        self.reciprocal  = self._calcReciprocalModule(cell= self.cell)
        self.coarse.div  = self._calcGrid(distance= self.kdistance)
        self.qpoints.div = self._calcGrid(distance= self.qdistance)
        self.dense.div   = self._calcDenseGrid(k_div= self.coarse.div, 
                                               q_div= self.qpoints.div)

    def _calcReciprocalModule(self, cell):
        a = np.zeros((3,3)) #direct vectors
        b = np.zeros((3,3)) #reciprocal vectors
        nb = [0,0,0]    #reciprocal vectors' modules

        for i in range (0,3):
            a[i] = cell.structure.cell[i]
        
        pi = np.pi
        b[0] = 2*pi*np.cross(a[1],a[2])/np.dot(a[0],np.cross(a[1],a[2]))
        b[1] = 2*pi*np.cross(a[2],a[0])/np.dot(a[1],np.cross(a[2],a[0]))
        b[2] = 2*pi*np.cross(a[0],a[1])/np.dot(a[2],np.cross(a[0],a[1]))
        
        for i in range (0,3):
            nb[i] = np.sqrt(np.dot(b[i],b[i]))

        return tuple(nb)

    def _calcGrid(self, distance): 
        grid = [0,0,0]
        for i in range (0,3):
            module = self.reciprocal[i]
            grid[i] = int((module/distance) + 1)

        return tuple(grid)

    def _calcDenseGrid (self, k_div, q_div):
        grid = [0,0,0]
        for i in range (0,3):
            grid[i] = int(np.lcm(int(q_div[i]),int(k_div[i])))
            if grid[i] == k_div[i]:
                grid[i] = 2*grid[i]

        return tuple(grid)

    def load(self, database):
        db = pk.load(database, False)
        self.kdistance = db.dget(self.db_dict, 'kdistance')
        self.qdistance = db.dget(self.db_dict, 'qdistance')
        self.qpoints.load(database = database)        
        self.coarse.load(database = database)
        self.dense.load(database = database)

    def dump(self, database):
        db = pk.load(database, False)
        db.dadd(self.db_dict,('kdistance', self.kdistance))
        db.dadd(self.db_dict,('qdistance', self.qdistance))
        self.qpoints.dump(database = database)
        self.coarse.dump(database = database)
        self.dense.dump(database = database)        

class CriticalTemp(AttrDisplay):
    db_dict = 'calc_stage'

    calc_stage = {'database': False, 
     'w_scf1'  : False,  'w_scf2'  : False, 
     'w_ph'    : False, 'w_q2r'   : False, 
     'w_matdyn': False, 'w_lamb': False, 
     'r_scf1'  : False,  'r_scf2'  : False, 
     'r_ph'    : False, 'r_q2r'   : False, 
     'r_matdyn': False, 'r_lamb': False,
     'calc_grids': False}

    operation_modes = {'n': 'new directory', 
     'c': 'continue calculation', 
     'w': 'overwrite'}
    
    def __init__(self, infile, control = 'n'):
        self.op = self._verCommand(op= control)
        self._setInputInfo(infile= infile)
        self._setDir()

    def calculate(self):
        current_dir = os.getcwd()
        os.chdir(self.dir.work)
        self.cell.copyFile()
        self._setDatabase()
        self._setPseudo()
        self._setGrids()
        self._setPrograms()
        self._runPrograms()
        os.chdir(current_dir)

    def _verCommand(self, op):
        valid_modes = list(self.operation_modes.keys())
        try:
            valid_modes.index(op)
        except ValueError as err:
            print('Wrong: Invalid command "{}". Default mode activated.'.format(op))
            op = 'n'
    
        print('Check: Operation mode "{mode} - {description}" activated'.format(
            mode = op, description = self.operation_modes.get(op)))

        return op    
    
    def _verPathExist(self, path):
        exist = os.path.exists(path)
        if exist == True:
            new_path = path + '_new'
            path = self._verPathExist(new_path)

        return path

    def _setInputInfo(self, infile):
        infile = PurePath(infile)

        if infile.suffix == '.cif':  
            self.prefix = infile.stem
            self.cell = CellStructure(infile)
            if self.op == 'c':
                self.database = Database(infile.with_name(self.prefix + '_config.db'))
                self.db = True
            else:
                self.db = False
            
        elif infile.suffix == '.db':
            self.prefix = infile.stem.replace('_config','')
            self.cell = CellStructure(infile.with_name(self.prefix + '.cif'))
            self.database = Database(infile)
            self.db = True

        else:
            raise AssertionError('Input file with incompatible format')

    def _setDir(self):
        db = self.db
        op = self.op

        if op == 'n': 
            work_path = os.path.join(os.getcwd(), self.prefix)
            work_path = self._verPathExist(path= work_path)
            self.dir = Directory(work_path= work_path)

        elif op == 'c':
            self.dir = Directory(work_path=  self.database.dir)
            return

        elif op == 'w' and db: 
            self.dir = Directory(work_path=  self.database.dir) 
            self.dir.delete()

        elif op == 'w' and not(db):
            work_path = os.path.join(os.getcwd(), self.prefix)
            self.dir = Directory(work_path= work_path)
            if os.path.exists(work_path):
                self.dir.delete()

        self.dir.makedirs()
        
    def _setDatabase(self):
        db = self.db
        op = self.op
        
        if op == 'c': 
            self.load(database = self.database.file)
            
        if op == 'n' and db: 
            self.database.copyFile()

        elif op == 'n' and not(db): 
            self.database= self._mkeDatabase()
            
        elif op == 'w' and not(db): 
            self.database= self._mkeDatabase()

        self.calc_stage.update({'database':True})
        print('Check: Database fineshed')

    def _mkeDatabase(self):

        db_name = self.prefix + '_config.db'

        try:
            subprocess.run(['write.py', db_name ])
        except OSError as err:
            print("Execution failed:", err, file=sys.stderr)
            raise
        else:
            print('Check: Database created')
        
        database = Database(os.path.join(self.dir.work, db_name))

        return database

    def _setGrids(self):
        self.grids = GridsWrapper(cell= self.cell)
        self.grids.load(database= self.database.file)

        if not(self.calc_stage.get('calc_grids')):
            self.grids.calculate()
            self.grids.dump(database= self.database.file)
            self.calc_stage.update({'calc_grids':True})
            self.dump(database= self.database.file)
        
        print(self.grids)

    def _setPseudo(self):
        self.pseudo = Pseudo()   
        self.pseudo.load(database= self.database.file)

    def _setPrograms(self):
        
        self.mpi    = Mpi()
        self.pw1    = Pwscf(prefix= self.prefix, grid= self.grids.dense, 
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
        print('-> Run programs')

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

    def load(self, database):
        db = pk.load(database, False)
        self.calc_stage = db.dgetall(self.db_dict)
    
    def dump(self, database):
        db = pk.load(database, False)
        db.set(self.db_dict, self.calc_stage)
        db.dump()

def get_input_data():
    user_input = sys.argv
    input_file = PurePath(user_input[1])

    return infile, command

################################################################################
##----------------------------------------------------------------------------##
################################################################################

if __name__ == '__main__':

    #-------test GridsWrapper--------------
    cell_Nb = CellStructure('/home/camila/Documentos/EMA/Program-TC/cif_database/Nb.cif')
    cell_H3S = CellStructure('/home/camila/Documentos/EMA/Program-TC/cif_database/H3S-200GPa.cif')

    print('Nb cell:')
    print(cell_Nb.structure.get_cell())

    print('H3S cell:')
    print(cell_H3S.structure.get_cell())

    grids = GridsWrapper(cell = cell_H3S, kdistance= 0.2, qdistance= 0.6)
    grids.calculate()
    print(grids)

    #-------test CriticalTemp--------------
    H3S_cif = '/home/camila/Documentos/EMA/Program-TC/cif_database/H3S-200GPa.cif'
    Nb_cif = '/home/camila/Documentos/EMA/Program-TC/cif_database/Nb.cif'

    TC = CriticalTemp(infile= Nb_cif , control= 'n')
    TC.calculate()

"""
    for cif in [Nb_cif, H3S_cif]:
        TC = CriticalTemp(infile= cif , control= 'w')
        TC.calculate()
    
    H3S_cif = '/home/camila/Documentos/EMA/Program-TC/teste/H3S-200GPa/H3S-200GPa.cif'
    Nb_cif = '/home/camila/Documentos/EMA/Program-TC/teste/Nb/Nb.cif'

    for cif in [H3S_cif, Nb_cif]:
        TC = CriticalTemp(infile= cif , control= 'c')
        TC.calculate()

    H3S_cif = '/home/camila/Documentos/EMA/Program-TC/teste/H3S-200GPa/H3S-200GPa_config.db'
    Nb_cif = '/home/camila/Documentos/EMA/Program-TC/teste/Nb/Nb_config.db'
    H3S_cif_new = '/home/camila/Documentos/EMA/Program-TC/teste/H3S-200GPa_new/H3S-200GPa.cif'

    for cif in [H3S_cif_new]:
        TC = CriticalTemp(infile= cif , control= 'w')
        TC.calculate()
"""



