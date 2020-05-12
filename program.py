#!/usr/bin/env python3

##------------------------------------------------------------------------------
# CNPEM - National Brazilian Center for Research in Energy and Materials
# Sirius - Research Group EMA
#
# Code project: TC_caculus
# 
# Objectif:
# Automate the use of Quantum Espresso (QE) for the calculus of 
# Superconductivity Critical Temperature (Tc) of diferent molecules and 
# cell structures
# 
##------------------------------------------------------------------------------

#MODULES

import numpy as np
import ase
import ase.io
import pickledb as pk
import sys
import time
import subprocess
import os
import shutil

from pathlib import PurePath

#FUNCTIONS

def get_input_data():

    # Read the input from user and save the correspondent file information 
    # (file directory, file format, and prefix) in the database
    #  
    #  @param: - 'data': (str) input from user
    #  @return: no return. the function modify the database directly

    #get input information
    user_input = sys.argv

    #treat input command
    command = ver_command(user_input = user_input)

    #treat input file
    input_file = user_input[1]

    path_file = PurePath(input_file)

    file_name = path_file.stem
    file_path = str(path_file)
    file_format = path_file.suffix.strip('.')

    assert os.path.isfile(path_file), 'Entrance must be a file directory'    
    
    #construct input_info
    if file_format == 'cif' and command != 'c':
        input_info = {
        'prefix': file_name,
        'cell_file': file_path,
        'cell_format':'cif',
        'tag_config':False
        }
    
    elif file_format == 'cif' and command == 'c':
        prefix = file_name
        path_parent = str(path_file.parent)
        path_config = os.path.join(path_parent, prefix+'_config.db')

        assert os.path.isfile(path_config), 'Config file not found' 
        assert ver_finish_stage( stage = 'database', 
                                 database = path_config), 'Config entrance corrupted'
        print('Check: Database file found')

        input_info = {
        'prefix': file_name,
        'cell_file': file_path,
        'cell_format': 'cif',
        'tag_config':True,
        'config_dir': path_parent,
        'config_file': path_config,
        'config_name': prefix+'_config.db'
        }

    elif file_format == 'db':
        
        assert file_name.find('_config') != -1, 'Config entrance must be a database'
        assert ver_finish_stage( stage = 'database', 
                                 database = file_path), 'Config entrance corrupted'
        print('Check: Database file found')
        
        prefix = file_name.replace('_config','')
        path_parent = str(path_file.parent)
        path_structure = os.path.join(path_parent, prefix+'.cif')

        assert os.path.isfile(path_structure), 'Could not find .cif file in {}'.format(path_parent)    

        input_info = {
        'prefix': prefix,
        'cell_file': path_structure,
        'cell_format':'cif',
        'tag_config':True,
        'config_dir': path_parent,
        'config_file': file_path,
        'config_name':file_name + '.db'
        }
                   
    else:
        raise AssertionError('Input file with incompatible format')
    
    input_info.update({'op_mode':command})

    return input_info

##----------------------------------------------------------------------------##
def ver_path_exist(path):

    exist = os.path.exists(path)
    if exist == True:
        new_path = path + '_new'
        path = ver_path_exist(new_path)

    return path

##----------------------------------------------------------------------------##
def ver_command(user_input):

    #get input command
    if len(user_input) > 2:
        command = user_input[2]
    else:
        command = 'n'

    #verify operation mode
    op_mode = command

    operation_modes = {
    'n': 'new directory', 
    'c': 'continue calculation', 
    'w': 'overwrite' 
    }
    
    valid_modes = list(operation_modes.keys())
    
    try:
        valid_modes.index(op_mode)
    except ValueError as err:
        print('Wrong: Invalid command "{}". Default mode activated.'.format(op_mode))
        op_mode = 'n'
 
    value = operation_modes.get(op_mode)

    print('Check: Operation mode "{mode} - {extended}" activated'.format(
        mode = op_mode, extended = value))

    return op_mode
##----------------------------------------------------------------------------##

def mk_work_environment(input_info):

    # 
    #  
    #  @param:  
    #  @return: 

    prefix = input_info.get('prefix')
    test_config = input_info.get('tag_config')
    op_mode = input_info.get('op_mode')
    current_path = os.getcwd()
    new_dir = False
    flush_dir = False

    if test_config and op_mode == 'c': 
        work_path = input_info.get('config_dir')

    elif test_config and op_mode == 'n': 
        work_path = os.path.join(current_path, prefix)
        work_path = ver_path_exist(work_path)
        new_dir = True

    elif test_config and op_mode == 'w': 
        work_path = input_info.get('config_dir')
        flush_dir = True

    elif not(test_config) and op_mode == 'c': 
        work_path = input_info.get('config_dir')

    elif not(test_config) and op_mode == 'n': 
        work_path = os.path.join(current_path, prefix)
        work_path = ver_path_exist(work_path)
        new_dir = True

    elif not(test_config) and op_mode == 'w': 
        work_path = os.path.join(current_path, prefix)
        exist = os.path.exists(work_path)
        if exist:
            flush_dir = True
        else:
            new_dir = True
 
    dir_names = ['calc_dir', 'out_dir', 'input_dir']
    dir_info = {}
    for name in dir_names:
        path = os.path.join(work_path, name)
        dir_info.update({name:path})
        if new_dir:
            os.makedirs(path)
        elif flush_dir:
            shutil.rmtree(path)
            os.makedirs(path)

    dir_info.update({'config_dir':work_path})

    os.chdir(work_path)

    return dir_info
##----------------------------------------------------------------------------##
def finish_stage (stage, database):

    #Verifica se está no diretório do banco de dados e faz a alteração
    
    curent_dir = os.getcwd()
    dir_name = os.path.basename(curent_dir)
    if dir_name.find('_dir') > -1:
        work_dir = os.path.dirname(curent_dir)
        os.chdir(work_dir)

    db = pk.load(database, False)
    db.dadd('calc_stage', (stage, True))
    db.dump()
    
    os.chdir(curent_dir)

    return 
##----------------------------------------------------------------------------##
def ver_finish_stage (stage, database):

    curent_dir = os.getcwd()
    dir_name = os.path.basename(curent_dir)
    if dir_name.find('_dir') > -1:
        new_dir = os.path.dirname(curent_dir)
        os.chdir(new_dir)

    db = pk.load(database, False)
    test = db.dget('calc_stage', stage)
    db.dump()
    
    os.chdir(curent_dir)

    return test
##----------------------------------------------------------------------------##
def construct_database(input_info, dir_info):

    test_config = input_info.get('tag_config')
    op_mode = input_info.get('op_mode')

    if op_mode == 'c': 
        database = input_info.get('config_name')
        return database

    elif test_config and op_mode == 'n': 
        database = copy_data_base(input_info = input_info)

    elif test_config and op_mode == 'w':
        database = input_info.get('config_name')

    elif not(test_config) and op_mode == 'n': 
        database = mk_data_base(input_info = input_info)
        
    elif not(test_config) and op_mode == 'w': 
        database = mk_data_base(input_info = input_info)

    set_dir_info(dir_info = dir_info, database = database)
    set_calc_stage (database = database)
    finish_stage (stage = 'database', database = database)

    print('Check: Database fineshed')

    return database
##----------------------------------------------------------------------------##

def mk_data_base(input_info):
    
    # Order of the 'set' functions matters
    #
    #  @param:  
    #  @return:

    prefix =  input_info.get('prefix')
    database = prefix + '_config.db'

    #Construct database
    try:
        subprocess.run(['write.py', database])
    except OSError as err:
        print("Execution failed:", err, file=sys.stderr)
        raise
    else:
        print('Check: Database created')
    
    set_cell_structure(input_info = input_info, database = database)

    set_fixed_values(prefix = prefix, database = database)

    try:
        set_pseudo_potencials(database = database)
    except KeyError as err:
        print('Pseudopotentials failed:', err, file=sys.stderr)
        raise
    else:
        print('Check: Pseudopotentials found')

    set_program_parameters(prefix = prefix, database = database)
    
    return database

##----------------------------------------------------------------------------##
def copy_data_base(input_info):
    
    # Order of the 'set' functions matters
    #
    #  @param:  
    #  @return:

    database = input_info.get('config_name')
    config_file = input_info.get('config_file')
    current_path = os.getcwd()

    # adding exception handling
    try:
        shutil.copy(config_file, current_path)
    except IOError as e:
        print('Unable to copy config.db file. %s' % e)
        raise
    except:
        print('Unexpected error:', sys.exc_info())
        raise
    else: 
        print('Check: Database file copied')

    set_cell_structure(input_info = input_info, database = database)

    return database

##----------------------------------------------------------------------------##
def copy_cell_structure(input_info):

    file_dir = input_info.get('cell_file')

    current_path = os.getcwd()

    # adding exception handling
    try:
        file_path = shutil.copy(file_dir, current_path)
    except shutil.SameFileError as e:
        file_path =  os.path.join(current_path, os.path.basename(file_dir))
    except IOError as e:
        print('Unable to copy cell structure file. %s' % e)
        raise
    except:
        print('Unexpected error:', sys.exc_info())
        raise
    
    print('Check: Cell structure file found')
    
    return file_path
##----------------------------------------------------------------------------##
def set_cell_structure(input_info, database):

    file_format = input_info.get('cell_format')
    file_path = copy_cell_structure(input_info = input_info)

    db = pk.load(database, False)

    #Cell structures information
    db.dcreate ('cell_structure')
    db.dadd('cell_structure',('dir',file_path) )
    db.dadd('cell_structure',('format',file_format) )

    db.dump()
    
    return
##----------------------------------------------------------------------------##

def set_dir_info (dir_info, database):

    calc_dir = dir_info.get('calc_dir')
    out_dir = dir_info.get('out_dir')
    input_dir = dir_info.get('input_dir')
    work_dir = dir_info.get('config_dir')

    db = pk.load(database, False)

    #Directories information
    db.dadd('dir',('calc_dir',calc_dir) )
    db.dadd('dir',('out_dir',out_dir) )
    db.dadd('dir',('input_dir',input_dir) )
    db.dadd('dir',('config_dir',work_dir) )

    db.dump()

    return
##----------------------------------------------------------------------------##
def set_calc_stage (database):

    db = pk.load(database, False)

    #Calculation stage
    db.dcreate('calc_stage')
    db.dadd('calc_stage',('database', False ))
    db.dadd('calc_stage',('w_scf1', False ))
    db.dadd('calc_stage',('w_scf2', False))
    db.dadd('calc_stage',('w_ph', False))
    db.dadd('calc_stage',('w_q2r', False))
    db.dadd('calc_stage',('w_matdyn', False))
    db.dadd('calc_stage',('w_lambda', False))
    db.dadd('calc_stage',('r_scf1', False ))
    db.dadd('calc_stage',('r_scf2', False))
    db.dadd('calc_stage',('r_ph', False))
    db.dadd('calc_stage',('r_q2r', False))
    db.dadd('calc_stage',('r_matdyn', False))
    db.dadd('calc_stage',('r_lambda', False))

    db.dump()

    return
##----------------------------------------------------------------------------##
def set_fixed_values(prefix, database):

    db = pk.load(database, False)

    #Input files' names
    db.dcreate('input_name')
    db.dadd('input_name',('scf1', prefix + '.scf1.in' ))
    db.dadd('input_name',('scf2', prefix + '.scf2.in' ))
    db.dadd('input_name',('ph', prefix + '.ph.in' ))
    db.dadd('input_name',('q2r', prefix + '.qr2.in' ))
    db.dadd('input_name',('matdyn', prefix + '.matdyn.in' ))
    db.dadd('input_name',('lambda', prefix + '.lambda.in' ))

    #Output files' names
    db.dcreate('output_name')
    db.dadd('output_name',('scf1', prefix + '.scf1.out' ))
    db.dadd('output_name',('scf2', prefix + '.scf2.out' ))
    db.dadd('output_name',('ph', prefix + '.ph.out' ))
    db.dadd('output_name',('q2r', prefix + '.qr2.out' ))
    db.dadd('output_name',('matdyn', prefix + '.matdyn.out' ))
    db.dadd('output_name',('lambda', prefix + '.lambda.out' ))

    #Files' strutures
    db.dcreate ('file_order')
    db.dadd('file_order',('pw',[
            'prefix','restart_mode','pseudo_dir','outdir','occupations',
            'smearing','degauss','ecutwfc','ecutrho','conv_thr']))
    db.dadd('file_order',('ph',[
            'prefix', 'outdir', 'tr2_ph', 'fildyn','ldisp', 'nq1',
            'nq2', 'nq3', 'electron_phonon','fildvscf','recover', 'el_ph_sigma', 
            'el_ph_nsigma']))
    db.dadd('file_order',('q2r',[
            'zasr', 'fildyn', 'flfrc', 'la2F']))
    db.dadd('file_order',('matdyn',[
            'asr', 'flfrc', 'flfrq', 'la2F', 'dos', 'fldos', 'nk1',
            'nk2','nk3','ndos']))

    #Programs' input files section names
    db.dcreate ('section_name')
    db.dadd('section_name',('ph','inputph'))
    db.dadd('section_name',('q2r','input'))
    db.dadd('section_name',('matdyn','input'))

    #Get all pseudo files
    pseudo_path = db.dget('dir','pseudo_folder')
    pseudo_files = os.listdir(pseudo_path)

    db.dcreate ('pseudo_dirs')
    for file_name in pseudo_files:
        element = file_name.split('.')[0]
        db.dadd('pseudo_dirs',(element,file_name))

    db.dump()

    return

##----------------------------------------------------------------------------##

def set_program_parameters(prefix, database):

    db = pk.load(database, False)

    #CREATE PREFIX DEPENDEND PARAMETERS
    #pw

    pseudo_path = db.dget('dir', 'pseudo_folder')

    try:
        ecutwfc, ecutrho = set_ecut(database = database)
    except ValueError as err:
        print('Ecut failed:', err, file=sys.stderr)
        raise
    else:
        print('Check: Cutoff energies found')

    #ph
    fildyn = prefix + '.dyn'
    fildvscf = prefix + '.dv'
    outdir = './out'

    #q2r
    flfrc = prefix + '.frc'

    #matdyn    
    flfrq= prefix+'.freq'    
    fldos= prefix+'.phonon.dos'

    #ADD PROGRAMS' PARAMETERS
    db.dadd('pw_par',('prefix',prefix) )
    db.dadd('pw_par',('outdir',outdir) )
    db.dadd('pw_par',('pseudo_dir', pseudo_path) )
    db.dadd('pw_par',('ecutwfc',ecutwfc) )
    db.dadd('pw_par',('ecutrho',ecutrho) )

    db.dadd('ph_par',('prefix',prefix) )
    db.dadd('ph_par',('fildyn',fildyn) )
    db.dadd('ph_par',('fildvscf',fildvscf) )
    db.dadd('ph_par',('outdir',outdir) )

    db.dadd('q2r_par',('flfrc',flfrc) )
    db.dadd('q2r_par',('fildyn',fildyn) )

    db.dadd('matdyn_par',('flfrc',flfrc) )
    db.dadd('matdyn_par',('flfrq',flfrq) )
    db.dadd('matdyn_par',('fldos',fldos) )

    db.dump()

    return

##----------------------------------------------------------------------------##

def set_pseudo_potencials(database):

    # Prepaire the pseudo potencials dir used to write th espresso-in scf 
    # input file. 
    #
    # - Read the cell_elements from atoms object and save them in a list
    # - Remove cell_elements list repetition 
    # - Get pseudopotentials directories for each cell element from database
    # - Group pseudopotencials of interest in a dictionary
    #  
    #  @param:  - 'cell_structure': (atoms obeject) that contains all cell_elements 
    #  @return: - 'pseudo_dir'    : (dict) cell_elements and correspondent 
    #                               directories of pseudopotencial files 
    #                               to be used to create pw (scf) input files
    
    db = pk.load(database, False)
    
    #get cell elements
    structure = get_cell_structure(database = database)
    cell_elements = structure.get_chemical_symbols()
    cell_elements = list( dict.fromkeys(cell_elements) )#remove repeated

    #save pseudo_files of interest
    db.dcreate ('pseudo_pw')
    for element in cell_elements:
        pseudo_file = db.dget('pseudo_dirs', element)
        db.dadd('pseudo_pw',(element, pseudo_file))

    db.dump()

    return 
##----------------------------------------------------------------------------##

def set_ecut (database):

    # Get ecutwfc and ecutrho values based on pseudopotencial files suggestions
    #
    # - Read all pseudopotencial files and save their content in the list
    # 'pseuto_content'
    # - Read each pseudopotencial content and ge ecut information saved in the 
    # lists: 'ecutwfc_values' and 'ecutrho_values'
    # - Get the maximmum values suggested for each variable
    #  
    #  @param:  - 'pseudo' : (dict) contain the directory of each element's 
    #                        pseudopotenail file
    #  @return: - 'ecutwfc': (int) maximum value suggested for ecutwfc variable
    #             'ecutrho': (int) maximum value suggested for ecutrho variable

    db = pk.load(database, False)
    pseudo = db.dgetall('pseudo_pw')
    pseudo_path = db.dget('dir', 'pseudo_folder')

    pseudo_dir = pseudo.values()

    pseudo_files =[]
    for file_name in pseudo_dir:
        file_directory = os.path.join(pseudo_path, file_name)
        with open(file_directory, 'r') as f:
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
    
    return ecutwfc, ecutrho
##----------------------------------------------------------------------------##

def get_cell_structure(database):

    # Get info for Writting input pw.x input files. 
    #
    # - Read cell_structure user input file and save atoms object
    # - Set pseudopotencials acoordinly with cell_structure cell_elements 
    # - Get ecut energies from pseudopotential files suggestions 
    #
    #  @param:  - 
    #  @return: -
  
    db = pk.load(database, False)

    #get cell structure from .cif file
    file_dir = db.dget('cell_structure','dir')
    cell_structure_format = db.dget('cell_structure','format')

    structure = ase.io.read(filename = file_dir,
                            format=cell_structure_format,
                            subtrans_included = False,
                            primitive_cell= True)

    return structure
##----------------------------------------------------------------------------##

def isnumber (number):

    # Test if a string represents a number
    #
    #  @param:  - 'number': (str) to be tested     
    #  @return: - 'test': (bool) See example     
    #
    # Ex: a = '0.05', b = '10e-4', c = '-100', d = '2', e = 'abc'
    #     - a,b,c and d return 'True'
    #     - e return 'False 
           
    number = number.replace('e','')
    number = number.replace('-','')
    number = number.replace('.','')

    test = number.isnumeric()

    return test
##----------------------------------------------------------------------------##    

def add_parameter_qe (file_dir, parameter, section):

    # Add a parameter to a input Quantum Espresso file at the indicated 
    # section.
    #
    # - File content is saved in a list
    # - Find position of the section in the file list
    # - Add "" if needed 
    # - Set identation of new information to be added 
    # - Insert information in the list at corect section
    # - Rewrite the content list in the file.
    #  
    #  @param: - 'file_dir' : (str) directory of input file to be modified
    #          - 'parameter': (str tuple) key and value to be added
    #          - 'section'  : (str) section name of the QE input file where the parameter 
    #                         should be added
    #  @return: no return.        

    with open(file_dir, 'r') as f:
       content = f.readlines()

    for line in content:
        if line.find(section.upper() or section.lower()) > -1:
            section_pos = content.index(line)
            break
    
    key = str(parameter[0])
    value = str(parameter[1])

    if value.find('true') ==-1 and value.find('false') == -1:
        if isnumber(value) == False:
            value = "'"+value+"'"

    identation = ''
    i = 0
    while i < 17 - len(key):
        identation += ' '
        i += 1

    new_pos = section_pos + 1
    new_info = '   ' + key + identation + '= ' + value + '\n'
    
    content.insert(new_pos,new_info)

    with open(file_dir, 'w') as f:     
            for line in content:
                f.write(line)

    return
##----------------------------------------------------------------------------##

def write_espresso_in_pw ():

    # Write input pw.x file for scf = scf1 or scf = scf2.
    #
    # - Get the parametrs from database accordinly with scf tag
    # - Write output file
    # 
    #  @param: - scf: (str) indicates if it will be genarated scf1 or 
    #                  scf2 file. Must be 'scf1' or 'scf2'.
    #  @return: no return.
    
    structure = get_cell_structure(database = _DATABASE)

    #get data from database
    db = pk.load(_DATABASE, False)


    prefix = db.dget('pw_par','prefix')
    pseudo = db.dgetall('pseudo_pw')
    pw_parameters = db.dgetall('pw_par')
    input_dir = db.dget('dir', 'input_dir')
    work_dir = db.dget('dir', 'config_dir')
    #scf1
    kgrid1 = db.dget('grids','kdense_div')
    koffset1 = db.dget('grids','kdense_off')
    file_name1 = db.dget('input_name', 'scf1')
    #scf2
    kgrid2 = db.dget('grids','kcoarse_div')
    koffset2 = db.dget('grids','kcoarse_off')
    file_name2 = db.dget('input_name', 'scf2')

    #set work directory    
    os.chdir(input_dir)

    #make scf1 input QE file

    if (not ver_finish_stage(stage = 'w_scf1', database = _DATABASE)):
        ase.io.write(filename = file_name1,
                    images = structure,
                    format = 'espresso-in',
                    input_data = pw_parameters,
                    pseudopotentials = pseudo,
                    kpts = kgrid1,
                    koffset = koffset1)

        add_parameter_qe(file_dir = file_name1, 
                        parameter = ('la2F','.true.'), 
                        section = 'system')

        finish_stage(stage = 'w_scf1', database = _DATABASE)
    else:
        print('input scf1 found!')

    if (not ver_finish_stage(stage = 'w_scf2', database = _DATABASE)):
        ase.io.write(filename = file_name2,
                    images = structure,
                    format = 'espresso-in',
                    input_data = pw_parameters,
                    pseudopotentials = pseudo,
                    kpts = kgrid2,
                    koffset = koffset2)

        finish_stage(stage = 'w_scf2', database = _DATABASE) 
    else:
        print('input scf2 found!')
    
    print('write scf!!')

    #return work directory    
    os.chdir(work_dir)

    return
##----------------------------------------------------------------------------##

def write_espresso_in(program):

    # Write espresso_in file accordinly with the program
    # 
    # - Get database info
    # - Create input file 
    #
    #  @param: 'program': (str) identifier of the program input file to be 
    #                     written. Vslues possible: 'ph', 'matdyn', 'q2r'. 
    #  @return: no return.
    
    if ver_finish_stage(stage = 'w_'+program, database = _DATABASE):
        text = 'input {program_txt} found!'
        print(text.format(program_txt = program))
        return

    qe_dict = program + '_par'
   
    #get data from database
    db = pk.load(_DATABASE, False)

    prefix = db.dget('pw_par','prefix')
    parameters = dict(db.dgetall(qe_dict).items())
    section_name = db.dget('section_name', program)
    file_order = db.dget('file_order', program)
    file_name = db.dget('input_name', program)
    input_dir = db.dget('dir', 'input_dir')
    work_dir = db.dget('dir', 'config_dir')

    #set work directory  
    os.chdir(input_dir)

    #create and write file
    first_line = '&'+section_name.upper()+'\n/\n'
    with open(file_name, 'w') as f:     
        f.write(first_line)

    file_order.reverse()
    for key in file_order:
        value = str(parameters.get(key))
        add_parameter_qe(file_dir = file_name, 
                        parameter = (key,value), 
                        section = section_name) 

    finish_stage(stage = 'w_'+program, database = _DATABASE) 
    
    #return work directory    
    os.chdir(work_dir)  

    return 
##----------------------------------------------------------------------------##

def run_espresso_in(program):

    if ver_finish_stage(stage = 'r_'+program, database = _DATABASE):
        text = '{program_txt} previously run!'
        print(text.format(program_txt = program))
        return

    # Run espresso_in file
    # 
    #  @param: 'program': (str) identifier of the program input file to be 
    #                     run. Vslues possible: 'ph', 'matdyn', 'q2r'. 
    #  @return: no return.

    #Get info from database
    db = pk.load(_DATABASE, False)
    
    input_dir = db.dget('dir', 'input_dir')
    work_dir = db.dget('dir', 'config_dir')
    calc_dir = db.dget('dir', 'calc_dir')
    out_dir = db.dget('dir', 'out_dir')   
    dir_qe = str(db.dget('dir','qe_programs'))    
    input_name = str(db.dget('input_name', program))
    output_name = str(db.dget('output_name', program))
    np = str(db.dget('mpi','np'))  
    nk = str(db.dget('mpi','nk'))
    
    #set work directory    
    os.chdir(calc_dir)

    input_name = os.path.join('..', input_dir, input_name)

    #set and verify QE directory    
    if program.find('scf') > -1:
        dir_program = os.path.join(dir_qe, 'pw.x')
    else:
        dir_program = os.path.join(dir_qe, program + '.x')

    if dir_qe:
        assert os.path.isdir(dir_program), 'Could not find Quantum Espresso'
    else:
        p = subprocess.run(['which', dir_program], stdout = subprocess.PIPE, 
                           stderr = subprocess.PIPE)
        assert str(p.stdout).split("'")[1], 'Could not find Quantum Espresso'        
    
    #make command
    if program.find('lambda' or 'q2r') > -1:
        command = dir_program + ' < ' + input_name
    elif program.find('matdyn') > -1:
        command = 'mpiexec'+' -np '+np+' '+dir_program+' -in '+input_name
    else:
        command = 'mpiexec'+' -np '+np+' '+dir_program+' -nk '+nk+' -in '+input_name

    #run program
    text = '\nrunning {program_txt}...'
    print(text.format(program_txt = program))
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
        #set work directory    
        os.chdir(out_dir)

        with open(output_name, 'w') as f:
            f.write(str(process.stdout))
        text = '{program_txt} fineshed! :D'
        print(text.format(program_txt = program))
        
        finish_stage(stage = 'r_'+program, database = _DATABASE) 

    #return work directory    
    os.chdir(work_dir)
    
    return 
##----------------------------------------------------------------------------##

def get_multplicity(dyn_files):

    # Get multiplicity values from dyn files. (ph.x output)
    # 
    #  @param:  - 'dyn_files'          : (str list list) 
    #  @return: - 'multiplicity_values': (int list)
        
    multiplicity_values = []
    for element in dyn_files:
        multiplicity = 0
        for line in element:
            if line.find('Dynamical  Matrix in cartesian axes') > -1:
                multiplicity += 1
        multiplicity_values.append(multiplicity)

    return multiplicity_values
##----------------------------------------------------------------------------##

def get_upper_freq(dyn_files):

    # Get upper frequence from dyn files. (ph.x output)
    # 
    #  @param:  - 'dyn_files' : (str list list)
    #  @return: - 'upper_freq': (int)
        
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
##----------------------------------------------------------------------------##

def write_espresso_in_lambda():

    if ver_finish_stage(stage = 'w_lambda', database = _DATABASE):
        print('\nlambda.in file found!')
        return

    # Write espresso_in lambda file
    # 
    # - Load data from database 
    # - Get data from dyn files
    # - Construct lambda.in file
    #  
    #  @param: - no parameters
    #  @return: - no return
    
    #get data from database
    db = pk.load(_DATABASE, False)

    sigma = db.dget('lambda_par','sigma_omega')
    input_name = db.dget('input_name', 'lambda')
    mu =  db.dget('lambda_par','mu')
    file_dyn_prefix = db.dget('ph_par','fildyn')
    prefix = db.dget('pw_par','prefix')
    calc_dir = db.dget('dir', 'calc_dir')
    input_dir = db.dget('dir', 'input_dir')
    work_dir = db.dget('dir', 'config_dir')
    
    #set work directory    
    os.chdir(calc_dir)

    #get q points information
    #file .dyn0 contains q points information used in lambda.in file
    file_name = file_dyn_prefix + '0'
    with open(file_name, 'r') as f:
       content = f.readlines()
    content.pop(0)

    q_points_info = content 
    num_q_points = int(q_points_info[0])

    #read/get information from dyn files
    dyn_info = []
    for num in range(1,num_q_points+1):
        file_name = file_dyn_prefix + str(num)
        with open(file_name) as f:
            content = f.readlines()
        dyn_info.append(content)

    multiplicity = get_multplicity(dyn_files = dyn_info)
    upper_freq = get_upper_freq(dyn_files = dyn_info)

    #set lambda.in file parts
    part_1 = str(upper_freq)+'  '+str(sigma)+'  1\n'

    part_2 = q_points_info
    for i in range(1,num_q_points+1):
        sufix = '  '+str(multiplicity[i-1]) +'\n'
        part_2[i] = part_2[i].replace('\n', sufix)

    part_3 = []
    aux = 'elph_dir/elph.inp_lambda.'
    for num in range(1,num_q_points+1):
        line = aux + str(num) + '\n'
        part_3.append(line)
    part_3.append(str(mu)+'\n')

    #write lambda.in
    lambda_info =[]
    lambda_info.append(part_1)
    lambda_info.extend(part_2)
    lambda_info.extend(part_3)

    #set work directory    
    os.chdir(input_dir)

    with open(input_name, 'w') as f:     
        for line in lambda_info:
            f.write(line)

    finish_stage(stage = 'w_lambda', database = _DATABASE)

    print('\nlambda.in file wrote!')

    #return work directory    
    os.chdir(work_dir)  

    return
##----------------------------------------------------------------------------##

################################################################################
##----------------------------------------------------------------------------##
################################################################################

#treat input data and prepaire work environment
start_time = time.time()

input_info = get_input_data()
dir_info   = mk_work_environment(input_info = input_info)

_DATABASE = construct_database(input_info = input_info, 
                               dir_info = dir_info)

#write input files
write_espresso_in_pw()
write_espresso_in(program = 'ph')
write_espresso_in(program = 'q2r')
write_espresso_in(program = 'matdyn')

#run programs pw, ph, q2r and matdyn
run_espresso_in(program = 'scf1')
run_espresso_in(program = 'scf2')

print('\nTotal Execution time:', time.time()-start_time)
run_espresso_in(program = 'ph')
print('\nTotal Execution time:', time.time()-start_time)

run_espresso_in(program = 'q2r')
run_espresso_in(program = 'matdyn')

#Tc calculus
write_espresso_in_lambda()
run_espresso_in(program = 'lambda')

print('\nTotal Execution time:', time.time()-start_time)