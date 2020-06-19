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

from classtools import AttrDisplay
from pathlib import PurePath
import ase
import ase.io
import os
import shutil
import sys
import pickledb as pk
import subprocess
        
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

        db.dcreate (self.database)
        db.dadd(self.db_dict,('np', self.np) )
        db.dadd(self.db_dict,('nk', self.nk) )

        db.dump()

class File(AttrDisplay):
    """
    
    This class provide file interface used to manage database and
    structures files 
    """
    db_dict = 'file'

    def __init__(self, file_in):
        self.path = PurePath(file_in)
        self.file = str(self.path)
        self.name = self.path.stem
        self.format = self.path.suffix.strip('.')
        self.dir = str(self.path.parent)
        print(self.path)
        assert os.path.isfile(self.path), 'This is not a file directory'

    def load(self, database):
        db = pk.load(database, False)
        file_in = db.dget(self.db_dict,'dir')
        self.__init__(file_in = file_in)

    def dump(self, database):
        db = pk.load(database, False)

        db.dcreate(db_dict)
        db.dadd(self.db_dict,('dir', self.file) )
        db.dump()

    def copyFile(self):
        """

        This function copies the file represented by the File object 
        to the current work path and atualize the File object information
        accordinly with teh path of the copied file.
        """
        current_path = os.getcwd()
        try:
            new_file = shutil.copy(self.file, current_path)
        except shutil.SameFileError as e:
            new_file =  os.path.join(current_path, os.path.basename(self.file))
            print('Check: File already here')

        except IOError as e:
            print('Unable to copy file. %s' % e)
            raise
        except:
            print('Unexpected error:', sys.exc_info())
            raise
        else:
            print('Check: File copied')
            self.__init__(new_file)    

class CellStructure(File):
    """
    
    This class provides an interface for cell structure .cif files,
    easily reaching information useful for quantum espresso programs
    """
    db_dict = 'cell_structure'

    def __init__(self, file_str = ''):
        File.__init__(self, file_str)
        assert self.format == 'cif', 'Structure file must be .cif'
        self.structure = ase.io.read(filename = self.file,
                            format=self.format,
                            subtrans_included = False,
                            primitive_cell= True)

        symbols = self.structure.get_chemical_symbols()
        self.elements = list(dict.fromkeys(symbols))
    
    def copyFile(self):
        print('\n-> Copying structure file ...\n')
        File.copyFile(self)

class Database(File):
    """
    
    This class provides an interface for database .db files
    """
    def __init__(self, file_str = ''):
        File.__init__(self, file_str)
        assert self.format == 'db', 'Database file must be .db'
        #verifirar de está corrompido ou não (database fineshed = True)

    def copyFile(self):
        print('\n-> Copying database file ...\n')
        File.copyFile(self)

class Pseudo(AttrDisplay):
    db_dict = 'pseudo'

    def __init__(self, folder = None):
        self.folder = folder
    
    def load (self, database):
        db = pk.load(database, False)
        self.folder = db.dget(self.db_dict,'pseudo_folder')
        self._setFiles()
        db.dump()
    
    def dump(self, database):
        db = pk.load(database, False)
        db.dadd(self.db_dict,('pseudo_folder', self.folder))
        db.dump()

    def _setFiles(self):     
        pseudo_files = os.listdir(self.folder)
        self.files = {}
        for file_name in pseudo_files:
            element = file_name.split('.')[0]
            self.files.update({element:file_name})
    
    def selectFiles (self, elements):
        pseudo_dict = {}
        for element in elements:
            pseudo_file = self.files.get(element)
            pseudo_dict.update({element:pseudo_file})
        return pseudo_dict
            
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

class Program(AttrDisplay):
    name = None
    program = None
    file_order = None
    section = 'input'

    def __init__(self, prefix):  
        assert self.file_order, 'File parameters order must be defined!'
        assert self.name, "Files' name must be defined!"
        assert self.program, 'Program name must be defined!'
            
        self.db_dict = self.program.replace('.x','_par')
        self.input =  prefix + '.' + self.name + '.in'
        self.output =  prefix + '.' + self.name + '.out'
        self.parameters = {}
        self._setParameters(prefix)

    def _setParameters(self, prefix):
        assert False, '"setParameters must be defined!"'

    def set_Dir (self, dir):
        self.input = os.path.join(dir.input, self.input)
        self.output = os.path.join(dir.output, self.output)
        self.calc = dir.calc

    def command(self, mpi):
        command = 'mpiexec' + ' -np ' + mpi.np + ' ' + self.program +' -nk ' + mpi.nk + ' -in ' + self.input
        return command
    
    def load(self, database):
        db = pk.load(database, False)

        par = dict(db.dgetall(self.db_dict).items())
        keys = par.keys()
        
        for key in keys:
            self.parameters.update({key: par.get(key)})
    
    def run(self, mpi):
        command = self.command(mpi= mpi)

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
        first_line = '&'+self.section.upper()+'\n/\n'
        with open(self.input, 'w') as f:     
            f.write(first_line)

        self.file_order.reverse()
        for key in self.file_order:
            value = str(self.parameters.get(key))
            self._addParameter(key= key, value= value, section= self.section) 

    def _addParameter(self, key, value, section):

        with open(self.input, 'r') as f:
            content = f.readlines()

        value = self._treatValue(word = value)
        section_pos = self._getPosition(content= content, word= section)
        identation = self._setIdentation(word= key)

        new_pos = section_pos + 1
        new_info = '   ' + key + identation + '= ' + value + '\n'
        content.insert(new_pos,new_info)

        with open(self.input, 'w') as f:     
            for line in content:
                f.write(line)
    
    def _isnumber (self, number):
        number = number.replace('e','')
        number = number.replace('-','')
        number = number.replace('.','')
        test = number.isnumeric()
        return test
    
    def _setIdentation(self, word):
        identation = ''
        i = 0
        while i < 17 - len(word):
            identation += ' '
            i += 1
        return identation
    
    def _treatValue(self, word):
        if word.find('true') ==-1 and word.find('false') == -1:
            if self._isnumber(word) == False:
                word = "'"+word+"'"
        return word
    
    def _getPosition(self, content, word):
        for line in content:
            if line.find(word.upper() or word.lower()) > -1:
                position = content.index(line)
                break
        return position

class Pwscf (Program):
    name = 'scf1'
    program = 'pw.x'
    file_order = ['prefix','restart_mode','pseudo_dir','outdir','occupations',
                  'smearing','degauss','ecutwfc','ecutrho','conv_thr']
    
    def __init__(self, prefix, pseudo, cell, grid):
        Program.__init__(self, prefix= prefix)
        self.kgrid = grid.div
        self.koffset = grid.off
        self.cell = cell
        self.pseudo = pseudo.selectFiles(cell.elements)
        self.parameters.update({'pseudo_dir': pseudo.folder})
        self._setEcut()

    def _setParameters(self, prefix):
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
        ase.io.write(filename = self.input,
                     images = self.cell.structure,
                     format = 'espresso-in',
                     input_data = self.parameters,
                     pseudopotentials = self.pseudo,
                     kpts = self.kgrid,
                     koffset = self.koffset)
            
class Pwscf2(Pwscf):
    name = 'scf2'

    def write(self):
        Pwscf.write(self) 
        self._addParameter(key= 'la2F', value = '.true.', 
                               section = 'system')
                               
class Phonon (Program):
    name  = 'ph'
    program = 'ph.x'
    file_order = ['prefix', 'outdir', 'tr2_ph', 'fildyn','ldisp', 'nq1',
                  'nq2', 'nq3', 'electron_phonon','fildvscf','recover', 'el_ph_sigma', 
                  'el_ph_nsigma']
    section = 'inputph'

    def __init__(self, prefix, grid):
        Program.__init__(self, prefix= prefix)
        self.parameters.update({'nq1': grid.div[0]})
        self.parameters.update({'nq2': grid.div[1]})
        self.parameters.update({'nq3': grid.div[2]})

    def _setParameters(self, prefix):
        self.parameters.update({'prefix': prefix})
        self.parameters.update({'outdir': './out'})
        self.parameters.update({'fildyn': prefix+'.dyn'})
        self.parameters.update({'fildvscf': prefix+'.dv'})

class Q2r (Program):
    name  = 'q2r'
    program = 'q2r.x'
    file_order = ['zasr', 'fildyn', 'flfrc', 'la2F']
    
    def _setParameters(self, prefix):
        self.parameters.update({'fildyn': prefix+'.dyn'})
        self.parameters.update({'flfrc': prefix+'.frc'})
    
    def command(self, mpi):
        command = self.program + ' < ' + self.input
        return command

class Matdyn (Program):
    name = 'matdyn'
    program = 'matdyn.x'
    file_order = ['asr', 'flfrc', 'flfrq', 'la2F', 'dos', 'fldos', 'nk1',
                  'nk2','nk3','ndos']

    def _setParameters(self, prefix):
        self.parameters.update({'flfrc': prefix+'.frc'})
        self.parameters.update({'flfrq': prefix+'.freq'})
        self.parameters.update({'fldos': prefix+'.phonon.dos'})

    def command(self, mpi):
        command = 'mpiexec'+' -np '+ mpi.np +' '+ self.program +' -in '+ self.input
        return command

class Lambda (Program):
    name = 'lambda'
    program = 'lambda.x'
    file_order = ['zasr', 'fildyn', 'flfrc', 'la2F']

    def _setParameters(self, prefix):
        self.parameters.update({'fildyn': prefix+'.dyn'})

    def command(self, mpi):
        command = self.program + ' < ' + self.input
        return command

    def write(self):

        current_dir = os.getcwd()
        os.chdir(os.path.join(current_dir, self.calc))
        q_info, q_num = self._getQPoints()
        dyn_info = self._readFiles(prefix= self.parameters.get('fildyn'),
                                         sufix_list= range(1, q_num+1))
        os.chdir(current_dir)

        multiplicity = self._getMultplicity(dyn_files = dyn_info)
        up_freq = self._getUpperFreq(dyn_files = dyn_info)

        part_1 = str(up_freq)+'  '+str(self.parameters.get('sigma_omega'))+'  1\n'

        part_2 = q_info
        for i in range(1,q_num+1):
            sufix = '  '+str(multiplicity[i-1]) +'\n'
            part_2[i] = part_2[i].replace('\n', sufix)

        part_3 = []
        aux = 'elph_dir/elph.inp_lambda.'
        for num in range(1,q_num+1):
            line = aux + str(num) + '\n'
            part_3.append(line)
        part_3.append(str(self.parameters.get('mu'))+'\n')

        lambda_info =[]
        lambda_info.append(part_1)
        lambda_info.extend(part_2)
        lambda_info.extend(part_3)

        with open(self.input, 'w') as f:     
            for line in lambda_info:
                f.write(line)

    def _getMultplicity(self, dyn_files):
        multiplicity_values = []
        for element in dyn_files:
            multiplicity = 0
            for line in element:
                if line.find('Dynamical  Matrix in cartesian axes') > -1:
                    multiplicity += 1
            multiplicity_values.append(multiplicity)

        return multiplicity_values

    def _getUpperFreq(self, dyn_files):
        freq_info = []
        for element in dyn_files:
            for line in element:
                if line.find('freq') > -1:
                    freq = line.split('=')[1]
                    freq = freq.replace('[THz]','').strip()
                    freq_info.append(float(freq))
            
        max_freq = max(freq_info)
        upper_freq = int(max_freq+1)

        return upper_freq

    def _readFiles(self, prefix, sufix_list):
        info = []
        for suf in sufix_list:
            file_name = prefix + str(suf)
            with open(file_name) as f:
                content = f.readlines()
            info.append(content)

        return info
    
    def _getQPoints(self):
        file_name = self.parameters.get('fildyn') + '0'
        with open(file_name, 'r') as f:
            content = f.readlines()
            content.pop(0)
        q_info = content 
        q_num = int(q_info[0])

        return q_info, q_num

################################################################################
##----------------------------------------------------------------------------##
################################################################################

if __name__ == '__main__':    
    
    database = '/home/camila/Documentos/EMA/Program-TC/H3S-200GPa/H3S-200GPa_config.db'
    prefix = 'H3S'
    cell = CellStructure('/home/camila/Documentos/EMA/Program-TC/cif_database/H3S-200GPa.cif')

    mpi    = Mpi()
    pseudo = Pseudo()
    grid1  = Grid('dense')
    grid2  = Grid('coarse')
    qgrid = Grid('qpoints')

    for obj in [pseudo, grid1, grid2, qgrid, mpi]:
        obj.load(database= database)

    pw1    = Pwscf(prefix= prefix, grid= grid1, cell= cell, pseudo= pseudo)
    pw2    = Pwscf2(prefix= prefix, grid= grid2, cell= cell, pseudo= pseudo)
    ph     = Phonon(prefix= prefix, grid= qgrid)
    q2r    = Q2r(prefix= prefix)
    matdyn = Matdyn(prefix= prefix)
    lamb   = Lambda(prefix=prefix)

    for obj in [pw1, pw2, ph, q2r, matdyn, lamb]:
        obj.load(database= database)
    
    print("-> Test object's construction \n")
    for obj in [mpi, cell, pseudo, grid1, grid2, pw1, pw2, ph, q2r, matdyn, lamb]:        
        print(obj, '\n')

    print('-> Test CellStructure methods \n')
    cell.copyFile()

    print('-> Test Programs methods \n')
    for program in [pw1, pw2, ph, q2r, matdyn, lamb]:
        program.write()
        program.run(mpi=mpi)