#!/usr/bin/env python3

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

class CriticalTemp(AttrDisplay):
    db_dict = 'calc_stage'

    calc_stage = {'database': False, 
     'w_scf1'  : False,  'w_scf2'  : False, 
     'w_ph'    : False, 'w_q2r'   : False, 
     'w_matdyn': False, 'w_lamb': False, 
     'r_scf1'  : False,  'r_scf2'  : False, 
     'r_ph'    : False, 'r_q2r'   : False, 
     'r_matdyn': False, 'r_lamb': False}

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

    def _setPrograms(self):
        cell = self.cell
        pseudo = Pseudo()
        grid1  = Grid('dense')
        grid2  = Grid('coarse')

        for obj in [pseudo, grid1, grid2]:
            obj.load(database= self.database.file)

        self.pw1    = Pwscf(prefix= cell.name, grid= grid1, cell= cell, pseudo= pseudo)
        self.pw2    = Pwscf2(prefix= cell.name, grid= grid2, cell= cell, pseudo= pseudo)
        self.ph     = Phonon(prefix= cell.name)
        self.q2r    = Q2r(prefix= cell.name)
        self.matdyn = Matdyn(prefix= cell.name)
        self.lamb   = Lambda(prefix= cell.name)
        self.mpi    = Mpi()

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
    H3S_cif = '/home/camila/Documentos/EMA/Program-TC/cif_database/H3S-200GPa.cif'
    Nb_cif = '/home/camila/Documentos/EMA/Program-TC/cif_database/Nb.cif'

    for cif in [H3S_cif, Nb_cif]:
        TC = CriticalTemp(infile= cif , control= 'n')
        TC.calculate()

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




