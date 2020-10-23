#!/usr/bin/env python3
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

"""
from classtools import AttrDisplay
from models_qe import CellStructure
from models_program import Database
from models_program import Directory
from models_program import InputVariables

import calculations_qe as calc_qe
import tools
from pathlib import PurePath

import os
import pickledb as pk
import sys


class DatabaseFactory():     

    def make(self, op, prefix, inputvar, calc):
        
        db_name = prefix + '_config.db'

        if op == 'c': 
            database= self._find(db_name= db_name)

        if op == 'w' or op == 'n': 
            if os.path.exists(db_name):
                os.remove(db_name)
        
            if inputvar:
                database= self._new(inputvar= inputvar, db_name= db_name, 
                                            keys= calc.variables)
            else:
                database = calc.makeDatabase(name = db_name)

        return database

    def _find(self, db_name):

        work_dir = os.getcwd()
        db_dir = os.path.join(work_dir, db_name)

        if os.path.exists(db_dir):
            database = Database(db_dir)
        else:
            raise AssertionError('No database file found')

        return database

    def _new(self, inputvar, db_name, keys):

        sections, parameters = inputvar.read()

        _DATABASE = db_name

        db = pk.load(_DATABASE, True)

        for i in range(len(sections)):
            db.dcreate(sections[i])
            for key in keys[i]:
                value = parameters[i].get(key, False)
                if value:
                    db.dadd(sections[i], (key, value))
        
        db_dir = os.path.join(os.getcwd(), db_name)

        database = Database(infile= db_dir)

        return database
        

class DirectoryFactory:

    def make(self, op, prefix, infile, calc):

        dir_name = prefix + '_' + calc.sufix

        if op == 'c':
            directory = self._copy(infile= infile)
#           return directory

        elif op == 'n': 
            directory = self._new(name = dir_name)


        elif op == 'w':
            directory = self._overwrite(name = dir_name)
        
        directory.makedirs()

        return directory

    def _copy(self, infile):
        return Directory(workpath= infile.dir) 

    def _new(self, name):
        work_path = os.path.join(os.getcwd(), name)
        work_path = self._ver_path_exist(path= work_path)
         
        return Directory(workpath= work_path)

    def _overwrite(self, name):
        work_path = os.path.join(os.getcwd(), name)
        directory = Directory(workpath= work_path)
        if os.path.exists(work_path):
            directory.delete()
        
        return directory

    def _ver_path_exist(self, path):
        exist = os.path.exists(path)
        if exist == True:
            new_path = path + '_new'
            path = self._ver_path_exist(new_path)

        return path

class ControllerProgram(AttrDisplay):
    
    write_method = tools.WriteSections()
    drfactory = DirectoryFactory()
    dbfactory = DatabaseFactory()

    operation_modes = {'default':'n',
     'n': 'new directory', 
     'c': 'continue calculation', 
     'w': 'overwrite'}

    calculation_modes = {'default':'TC',
     'TC': 'superconductivity critical temperature',
     'EN': 'cell structure energy',
     'PH': 'phonons distribution in gamma'}


    def __init__ (self, clmode):
                                            
        infile, opmode, inputvar = self._get_input(commands_valid= self.operation_modes)
 
        self.cell = CellStructure(infile)
        self.prefix = PurePath(infile).stem
        self.opmode = self._ver_command(command= opmode, 
                                        commands_valid= self.operation_modes)
        self.clmode = self._ver_command(command= clmode, 
                                        commands_valid= self.calculation_modes)
        self.routine = self._set_calc_routine(clmode= self.clmode)
        
        if inputvar:
            self.inputvar = InputVariables(inputvar)
        else:
            self.inputvar = None
            
    def _set_calc_routine(self, clmode):

        if clmode == 'TC':
            routine = calc_qe.TC(structure= self.cell)
        elif clmode == 'EN':
            routine = calc_qe.Energy(structure= self.cell)
        elif clmode == 'PH':
            routine = calc_qe.PhononGamma(structure= self.cell)
        return routine

    def calculate(self):
    
        self.routine.directory = self.drfactory.make(op=self.opmode, 
                                               prefix= self.prefix, 
                                               infile= self.cell,
                                               calc= self.routine)
        
        os.chdir(self.routine.directory.work)
        self.cell.copy()
        self.routine.database= self.dbfactory.make(op=self.opmode, 
                                                   prefix= self.prefix, 
                                                   inputvar= self.inputvar, 
                                                   calc= self.routine)
                                                   
        if self.inputvar:
            self.inputvar.copy()
            self.inputvar.rename(name= self.routine.input_name)
        else:
            self._make_inputvar(database= self.routine.database, 
                            name= self.routine.input_name)
        
        try:       
            self.routine.load(database= self.routine.database.path)
        except KeyError:
            pass
        except Exception:
            raise(Exception)
   
        self.routine.calculate()

    def _make_inputvar(self, database, name):
        sections, parameters, keys = database.read()
        self.write_method.write(section= sections, filename= name, 
                               file_order=keys, parameters= parameters)
    

    def _get_input(self, commands_valid):
        """ 
        This function read input for CriticalTemp instantiation from user shell
        @return: control (str): command for critical temperature calculation
                infile (str) : path of input file
        """
        user_input = sys.argv

        infile = user_input[1]

        if len(user_input) == 2:
            operation = 'n'
            inputvar = None

        elif len(user_input) == 3:
            if  self._is_opmode(str= user_input[2]):
                operation = user_input[2]
                inputvar = None
            else:
                operation = 'n'
                inputvar = user_input[2]

        elif len(user_input) == 4:
            if  self._is_opmode(str= user_input[2]):
                operation = user_input[2]
                inputvar = user_input[3]
            else:
                operation = user_input[3]
                inputvar = user_input[2]
        
        else:
            raise AssertionError('Bad input argument')
        
        return  infile, operation, inputvar 

    def _is_opmode(self, str):
        return list(self.operation_modes.keys()).count(str)

    def _ver_command(self, command, commands_valid):
        commands_all = list(commands_valid.keys())
        try:
            commands_all.index(command)
        except ValueError as err:
            print('Wrong: Invalid command "{}". Default mode activated.'.format(command))
            command = commands_valid.get('default')
    
        print('Check: Operation mode "{mode} - {description}" activated'.format(
            mode = command, description = commands_valid.get(command)))

        return command  

################################################################################
##----------------------------------------------------------------------------##
################################################################################

if __name__ == '__main__':
    
    pass