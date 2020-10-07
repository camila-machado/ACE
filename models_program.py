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

class IFile(AttrDisplay):
    """
    
    This class provide file interface used to manage database and
    structures files 
    """
    
    def __init__(self, infile):

        if not(os.path.isfile(infile)):
            infile = os.path.join(os.getcwd(), infile)
            assert os.path.isfile(infile), 'This is not a file directory'

        path = PurePath(infile)
        self.name = path.stem
        self.format = path.suffix.strip('.')
        self.dir = str(path.parent)
        self.file = str(path)

    def copy(self):
        """

        This function copies the file represented by the File object 
        to the current work path and atualize the File object information
        accordinly with the copied file new path.
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
            self.__init__(infile= new_file)    

    def rename(self, name):
        new_file = os.path.join(self.dir, name)
        os.rename(src= self.file, dst=  new_file)
        self.__init__(infile= new_file)

class CellStructure(IFile):
    """
    
    This class provides an interface for cell structure .cif files,
    easily reaching information useful for quantum espresso programs
    """
    db_dict = 'cell_structure'

    def __init__(self, infile = ''):
        IFile.__init__(self, infile)
        assert self.format == 'cif', 'Structure file must be .cif'
        self.structure = ase.io.read(filename = self.file,
                                    format=self.format,
                                    subtrans_included = False,
                                    primitive_cell= True)
        symbols = self.structure.get_chemical_symbols()
        self.elements = list(dict.fromkeys(symbols))
    
    def copy(self):
        print('\n-> Copying structure file ...\n')
        IFile.copy(self)

class Database(IFile):
    """
    
    This class provides an interface for database .db files
    """
    def __init__(self, infile = ''):
        IFile.__init__(self, infile)
        assert self.format == 'db', 'Database file must be .db'

    def copy(self):
        print('\n-> Copying database file ...\n')
        IFile.copy(self)

    def read(self):
        
        db = pk.load(self.file, False)

        sections = list(db.getall())

        parameters = []
        keys = []
        for section in sections:
            parameters.append(dict(db.dgetall(section))) 
            keys.append(list(db.dgetall(section)))

        return sections, parameters, keys

class InputVariables(IFile):
    """
    
    This class provides an interface for input variables file
    """
    def __init__(self, infile = ''):
        IFile.__init__(self, infile)

        assert self.format == 'in', 'Input file must be .in'
        sec, par = self.read(infile= infile)
        self.sections = sec
        self.parameters = par

    def copy(self):
        print('\n-> Copying input variables file ...\n')
        IFile.copy(self)

    def read(self, infile):
        with open(infile, 'r') as f:
            file_content = f.readlines()

        sections = []
        parameters =[]
        
        par = {}
        for line in file_content:

            if line.count('&'):
                section = line.strip('&\n ').lower()
                sections.append(section)
                if par:
                    parameters.append(par.copy())
                    par.clear()               
            else:
                key, value = self._parseLine(line)
                if key:
                    par.update({key:value})

        parameters.append(par)

        return sections, parameters

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

