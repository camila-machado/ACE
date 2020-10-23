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

import ase
import ase.io
import os
import shutil
import sys
import pickledb as pk
from pathlib import PurePath
from classtools import AttrDisplay

class Directory:

    def __init__(self, workpath):
        self.work = workpath
        self.calc = os.path.join(workpath, 'calc')
        self.output = os.path.join(workpath, 'output')
        self.input = os.path.join(workpath, 'input')

    def makedirs(self):
        for dir in [self.calc, self.output, self.input]:
            try:
                os.makedirs(dir)
            except FileExistsError:
                print('directory {} already exist'.format(dir))       
   
    def delete(self):
        for dir in [self.calc, self.output, self.input]:
            try:
                shutil.rmtree(dir)
            except FileNotFoundError:
                print('directory {} do not exist'.format(dir))   

class IFile(AttrDisplay):
    """
    
    This class provide file interface used to manage database and
    structures files 
    """
    def __init__(self, filename):
        filepath = os.path.join(os.getcwd(), filename)
        filepath = os.path.normpath(filepath)
        assert os.path.isfile(filepath), 'This is not a file directory'

        path = PurePath(filepath)
        self.name = path.stem
        self.format = path.suffix.strip('.')
        self.dir = str(path.parent)
        self.path = str(path)

    def copy(self):
        """

        This function copies the file represented by the File object 
        to the current work path and atualize the File object information
        accordinly with the copied file new path.
        """
        current_path = os.getcwd()
        try:
            filenew = shutil.copy(self.path, current_path)
        except shutil.SameFileError as e:
            filenew =  os.path.join(current_path, os.path.basename(self.path))
            print('Check: File already here')
        except IOError as e:
            print('Unable to copy file. %s' % e)
            raise
        except:
            print('Unexpected error:', sys.exc_info())
            raise
        else:
            print('Check: File copied')
            self.__init__(filename= filenew)    

    def rename(self, name):
        filenew = os.path.join(self.dir, name)
        os.rename(src= self.path, dst=  filenew)
        self.__init__(filename= filenew)

class CellStructure(IFile):
    """
    
    This class provides an interface for cell structure .cif files,
    easily reaching information useful for quantum espresso programs
    """
    db_dict = 'cell_structure'

    def __init__(self, filename = ''):
        IFile.__init__(self, filename)
        assert self.format == 'cif', 'Structure file must be .cif'

    def copy(self):
        print('\n-> Copying structure file ...\n')
        IFile.copy(self)

    def read(self):
        atoms_structure = ase.io.read(filename = self.path,
                                    format=self.format,
                                    subtrans_included = False,
                                    primitive_cell= True)
        return atoms_structure
    
    def elements(self):
        atoms = self.read()
        symbols = atoms.get_chemical_symbols()
        elements = list(dict.fromkeys(symbols))

        return elements

class Database(IFile):
    """
    
    This class provides an interface for database .db files
    """
    def __init__(self, filename = ''):
        IFile.__init__(self, filename)
        assert self.format == 'db', 'Database file must be .db'

    def copy(self):
        print('\n-> Copying database file ...\n')
        IFile.copy(self)

    def read(self):
        db = pk.load(self.path, False)

        sections = list(db.getall())

        parameters = []
        keys = []
        for section in sections:
            parameters.append(dict(db.dgetall(section))) 
            keys.append(list(db.dkeys(section)))

        return sections, parameters, keys

class InputVariables(IFile):
    """
    
    This class provides an interface for input variables file
    """
    def __init__(self, filename = ''):
        IFile.__init__(self, filename)
        assert self.format == 'in', 'Input file must be .in'

    def copy(self):
        print('\n-> Copying input variables file ...\n')
        IFile.copy(self)

    def read(self):

        with open(self.path, 'r') as f:
            file_content = f.readlines()

        sections = []
        parameters = []
        keys = []

        sec_key = []
        sec_par = {}
        for line in file_content:

            if line.count('&'):
                section = line.strip('&\n ').lower()
                sections.append(section)
                if sec_par:
                    parameters.append(sec_par.copy())
                    sec_par.clear()               
            else:
                key, value = self._parseLine(line)
                if key:
                    sec_par.update({key:value})
                    sec_key.append(key)

        parameters.append(sec_par)
        keys.append(sec_key)

        return sections, parameters, keys

    def _parseLine(self, line):
        eq = line.find('=')
        if eq == -1:
            return 0, 0
        key = line[:eq].strip()
        value = line[eq+1:].strip()

        return key, self._parseValue(value)

    def _parseValue(self, expr):
        try:
            return eval(expr)
        except:
            return expr
        else:
            return expr

################################################################################
##----------------------------------------------------------------------------##
################################################################################

if __name__ == '__main__':    
     
    #Unit Tests
         
    #Test class Directory
    print('\nTEST DIRECTORY\n')

    dir = Directory(workpath= os.getcwd())
    print('\nTEST 1: criate a Directory obj\n', dir)

    print('\nTEST 2: criate subdirectories\n')
    dir.makedirs()

    print('\nTEST 3: delete subdirectories\n')
    dir.delete()

    #Test class CellStructure
    print('\nTEST CELLSTRUCTURE\n')

    cell1 = CellStructure(filename= '/home/camila/Documentos/EMA/Program-TC/cif_database/H3S.cif')
    print('\nTEST 1: criate a CellStructure obj from .cif file from absolute path\n', cell1)

    cell2 = CellStructure(filename='../models_qe/Nb.cif')
    print('\nTEST 2: criate a CellStructure obj from .cif file from relative path\n', cell2)

    print('\nTEST 3: copy CellStructure .cif file\n')
    cell1.copy()
    
    print('\nTEST 4: read a CellStructure obj\n')
    atoms = cell1.read()
    print(atoms)

    print('\nTEST 4: get elements from a CellStructure obj\n')
    elements = cell1.elements()
    print(elements)
    
    #Test class Database
    print('\nTEST DATABASE\n')

    database = Database(filename= 'Nb_config.db')
    print('\nTEST 1: criate a Database obj from .db file\n', database)

    print('\nTEST 2: read a Database obj\n')

    sec, par , key = database.read()

    print('database', database.name,":\n")
    print('sections:', sec, '\n')
    print('keys:', key, '\n')
    print('parameters:', par, '\n')
    
    #Test class InputVariables
    print('\nTEST INPUTVARIABLES\n')

    inputvar = InputVariables(filename= 'inputTC.in')
    print('\nTEST 1: criate a Inputvar obj from .in file\n', inputvar)

    print('\nTEST 2: read a Inputvar obj\n')

    sec, par , key = inputvar.read()

    print('inputvar', inputvar.name,":\n")
    print('sections:', sec, '\n')
    print('keys:', key, '\n')
    print('parameters:', par, '\n')
    
    
