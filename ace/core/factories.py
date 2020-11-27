# Copyright 2020, Camila Machado de Araújo
# (see accompanying license files for details).

"""
"""
import os
import sys
from pathlib import PurePath

import pickledb as pk

from util.classtools import AttrDisplay
from util.filendir import File
from util.filendir import Directory
import util.io as io
import core.calculations as calculations

def calcfactory(clmode, cell, prefix):

    if clmode == 'tc':
        routine = calculations.TC(structure= cell, 
                                  prefix= prefix)
    elif clmode == 'eos':
        routine = calculations.EquationOfState(structure= cell,
                                               prefix= prefix)
    elif clmode == 'ph':
        routine = calculations.PhononGamma(structure= cell, 
                                           prefix= prefix)
    else:
        assert False, 'Calculation not supported yet'

    return routine

def dbfactory(opmode, prefix, filvar, calc):
    
    db_name = prefix + '_config.db'

    if opmode == 'c': 
        database= _dbfind(db_name= db_name)

    if opmode == 'w' or opmode == 'n': 
        if os.path.exists(db_name):
            os.remove(db_name)
    
        if filvar:
            database= _dbnew(filvar= filvar, db_name= db_name, 
                                        keys= calc.variables)
        else:
            dir_database = calc.makeDatabase(name = db_name)
            database = File(filename= dir_database, 
                            fileformat= 'db',
                            read_strategy= io.ReadDatabase())

    return database

def _dbfind(db_name):

    dir_work = os.getcwd()
    dir_database = os.path.join(dir_work, db_name)

    if os.path.exists(db_dir):
        database = File(filename= dir_database, 
                        fileformat= 'db',
                        read_strategy= io.ReadDatabase())
    else:
        raise AssertionError('No database file found')

    return database

def _dbnew(filvar, db_name, keys):

    sections, parameters = filvar.read()

    _DATABASE = db_name

    db = pk.load(_DATABASE, True)

    for i in range(len(sections)):
        db.dcreate(sections[i])
        for key in keys[i]:
            value = parameters[i].get(key, False)
            if value:
                db.dadd(sections[i], (key, value))
    
    dir_database = os.path.join(os.getcwd(), db_name)

    database = File(filename= dir_database, 
                    fileformat= 'db',
                    read_strategy= io.ReadDatabase())

    return database
        

def dirfactory(opmode, prefix, filcell, clmode):

    dir_name = prefix + '_' + clmode.upper()

    if opmode == 'c':
        directory = _dircopy(filcell= infile)

    elif opmode == 'n': 
        directory = _dirnew(name = dir_name)

    elif opmode == 'w':
        directory = _diroverwrite(name = dir_name)
    
    directory.makedirs()

    return directory

def _dircopy(filcell):
    return Directory(workpath= filcell.dir) 

def _dirnew(name):
    work_path = os.path.join(os.getcwd(), name)
    work_path = _ver_path_exist(path= work_path)
        
    return Directory(workpath= work_path)

def _diroverwrite(name):
    work_path = os.path.join(os.getcwd(), name)
    directory = Directory(workpath= work_path)
    if os.path.exists(work_path):
        directory.delete()
    
    return directory

def _ver_path_exist(path):
    exist = os.path.exists(path)
    if exist == True:
        new_path = path + '_new'
        path = _ver_path_exist(new_path)

    return path


