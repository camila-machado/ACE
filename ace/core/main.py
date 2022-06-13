# Copyright 2020, Camila Machado de Araújo
# (see accompanying license files for details).

import os
import sys

import core.factories as fct
import util.io as io
from util.filendir import File

operations = {
    'default':'n',
    'n': 'new directory', 
    'c': 'continue calculation', 
    'w': 'overwrite'}

calculations = {
    'default':'tc',
    'tc': 'superconductivity critical temperature',
    'eos': 'cell equation of state',
    'ph': 'phonons distribution in gamma'}

def main():
    
    inputinfo = _get_input()

    clmode = inputinfo.get("clmode")
    opmode = inputinfo.get("opmode")
    filcell = inputinfo.get("filcell")
    filvar = inputinfo.get("filvar")

    _showcommand(command= clmode, allcommands= calculations)
    _showcommand(command= opmode, allcommands= operations)

    if filcell == None:
        assert False, "No cristalographic input file"

    cell = filcell.read()
    prefix = filcell.name
    print(prefix)
    
    routine = fct.calcfactory(clmode= clmode, cell= cell, prefix= prefix)

    #build environment
    routine.directory = fct.dirfactory(opmode= opmode, prefix= prefix, 
                                       filcell= filcell, clmode= clmode)
                        
    os.chdir(routine.directory.work)

    filcell.copy()

    routine.database = fct.dbfactory(opmode= opmode, prefix= prefix, 
                                     filvar= filvar, calc= routine)
    #set inputvar file                                            
    if filvar:
        filvar.copy()
        filvar.rename(name= routine.input_name)
    else:
        sections, parameters, keys = routine.database.read()
        io.write_sections(sections= sections, filename= routine.input_name, 
                          file_order=keys, parameters= parameters)

    #calculate routine
    try:       
        routine.load(database= routine.database.path)
    except KeyError:
        pass
    except Exception:
        raise(Exception)

    routine.calculate()
    routine.directory.clean()   

def _get_input():
    """ 
    This function read input for Calculation instantiation from user shell
    @return: control (str): command for critical temperature calculation
             infile (str) : path of input file
    """

    inputinfo = {
    "opmode":operations.get("default"),
    "clmode":calculations.get("default"),
    "filcell":None,
    "filvar":None}

    user_input = sys.argv[1:]

    flag_cell= True
    flag_var= True
    flag_op= True
    flag_cl = True

    if len(user_input) < 1 or len(user_input) > 4:
        raise AssertionError('Bad number of input arguments')

    for i in user_input:
        if _iscommand(command = i, allcommands= operations) and flag_op:
            inputinfo.update({"opmode":i})
            flag_op = False
    
        if _iscommand(command = i, allcommands= calculations) and flag_cl:
            inputinfo.update({"clmode":i})
            flag_cl = False

        if _isfile(fil = i, format= 'cif') and flag_cell:

            filcell = File(filename= i, fileformat= 'cif', read_strategy= io.ReadCif())

            inputinfo.update({"filcell": filcell})
            flag_cell = False

        if _isfile(fil = i, format= 'in') and flag_var:

            filvar = File(filename= i, fileformat= 'in', read_strategy= io.ReadInputvar())

            inputinfo.update({"filvar": filvar})
            flag_var = False

    return inputinfo

def _iscommand(command, allcommands):
    return list(allcommands.keys()).count(command)

def _isfile(fil, format):

    try:
        filformat = fil.split('.')[-1]
    except IndexError:
        return False
    else: 
        if os.path.isfile(fil) and filformat == format:
            return True
        else: 
            return False

def _showcommand(command, allcommands):
    print('Check: Operation mode "{mode} - {description}" activated'.format(
        mode = command, description = allcommands.get(command)))

    return command  




    
    