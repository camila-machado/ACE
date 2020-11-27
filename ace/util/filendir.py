# Copyright 2020, Camila Machado de Araújo
# (see accompanying license files for details).

"""

This module provide class interfaces for quantum espresso programs
"""
import os
import shutil
import sys
from pathlib import PurePath

import ase
import ase.io

from util.classtools import AttrDisplay
import util.io as io

class Directory(AttrDisplay):

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
            self._delete(dir)
    
    def clean(self):
        for dir in [self.calc, self.output, self.input]:
            if self._isempty(dir):
                self._delete(dir)

    def _isempty(self, dir):
        if os.path.exists(dir) and os.path.isdir(dir):
            if not os.listdir(dir):
                return True
            else:
                return False

    def _delete(self, dir):
        try:
            shutil.rmtree(dir)
        except FileNotFoundError:
            print('directory {} do not exist'.format(dir))   


class File(AttrDisplay):

    read_strategy = io.IReadStrategy()
    """
    
    This class provide file interface used to manage database and
    structures files 
    """
    def __init__(self, filename, fileformat, read_strategy):
             
        if os.path.isdir(filename):
            filepath = os.path.normpath(filename)
        elif os.path.isfile(filename):
            filepath = os.path.join(os.getcwd(), filename)
            filepath = os.path.normpath(filepath)
        else:
            assert (os.path.isfile(filepath) and os.path.isdir(filepath)), 'This is not a file directory'

        path = PurePath(filepath)
        self.file = os.path.basename(filepath)
        self.name = path.stem
        self.format = path.suffix.strip('.')
        self.dir = str(path.parent)
        self.path = str(path)
        self.read_strategy = read_strategy

        assert self.format == fileformat, 'File must be of .{} format'.format(fileformat)

    def rename(self, name):
        filenew = os.path.join(self.dir, name)
        os.rename(src= self.path, dst=  filenew)
        self.__init__(filename= filenew,
                      fileformat= self.format, 
                      read_strategy= self.read_strategy)

    def read(self):

        content = self.read_strategy.read(filepath= self.path)
        return content

    def copy(self):
        print('\n-> Copying {} file ...\n'.format(self.file))
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
            self.__init__(filename= filenew, 
                          fileformat= self.format, 
                          read_strategy= self.read_strategy)   


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

    cell1 = File(filename= '/home/camila/Documentos/EMA/REPOSITORY/exemples/cif_db/H3S.cif',
                 fileformat= 'cif', read_strategy= io.ReadCif())
    print('\nTEST 1: criate a CellStructure obj from .cif file from absolute path\n', cell1)

    cell2 = File(filename='../models_qe/Nb.cif',
                 fileformat= 'cif', read_strategy= io.ReadCif())
    print('\nTEST 2: criate a CellStructure obj from .cif file from relative path\n', cell2)

    print('\nTEST 3: copy CellStructure .cif file\n')
    cell1.copy()
    
    print('\nTEST 4: read a CellStructure obj\n')
    atoms = cell1.read()
    print(atoms)

    print('\nTEST 4: get elements from a CellStructure obj\n')
    elements = cell1.elements()
    print(elements)
    
    
    
