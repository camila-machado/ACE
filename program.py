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

def get_input_data(data):

    # Read the input from user and save the correspondent file information 
    # (file directory, file format, and prefix) in the database
    #  
    #  @param: - 'data': (str) input from user
    #  @return: no return. the function modify the database directly

    #get input information
    file_dir = os.path.normpath(data)
    assert os.path.isfile(file_dir), 'Entrance must be a file directory'    
    path = PurePath(file_dir)
       
    prefix = path.stem
    file_format = path.suffix.strip('.')

    assert file_format == 'cif', 'Entrance file must be .cif type'

    input_info = {
    'prefix': prefix,
    'file_dir':file_dir,
    'file_format':file_format
    }

    return input_info
##----------------------------------------------------------------------------##
def mk_work_environment(input_info):

    # 
    #  
    #  @param:  
    #  @return: 

    prefix = input_info.get('prefix')

    old_path = os.getcwd()
    new_path = os.path.join(old_path, prefix)

    #@ trecho a ser mudado -> "continuar de onde parou"
    exist = os.path.exists(new_path)
    if exist == True:
        shutil.rmtree(new_path)

    dir_names = ['calc_dir', 'out_dir', 'input_dir']
    dir_info = {}
    for name in dir_names:
        path = os.path.join(new_path, name)
        dir_info.update({name:path})
        os.makedirs(path)
    dir_info.update({'work_dir':new_path})

    os.chdir(new_path)

    return dir_info
##----------------------------------------------------------------------------##

def mk_data_base(input_info, dir_info):
    
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
        print('Database created')
    
    set_cell_structure(input_info = input_info, database = database)

    set_fixed_values(prefix = prefix, 
                     dir_info = dir_info, database = database)

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
def set_cell_structure(input_info, database):

    file_dir = input_info.get('file_dir')
    file_format = input_info.get('file_format')

    current_path = os.getcwd()

    # adding exception handling
    try:
        file_path = shutil.copy(file_dir, current_path)
    except IOError as e:
        print('Unable to copy cell structure file. %s' % e)
        raise
    except:
        print('Unexpected error:', sys.exc_info())
        raise
    
    print('Check: Cell structure file found')
    
    db = pk.load(database, False)

    #Cell structures information
    db.dcreate ('cell_structure')
    db.dadd('cell_structure',('dir',file_path) )
    db.dadd('cell_structure',('format',file_format) )

    db.dump()
    
    return
##----------------------------------------------------------------------------##

def set_fixed_values(prefix, dir_info, database):

    calc_dir = dir_info.get('calc_dir')
    out_dir = dir_info.get('out_dir')
    input_dir = dir_info.get('input_dir')
    work_dir = dir_info.get('work_dir')

    db = pk.load(database, False)

    #Directories information
    db.dadd('dir',('calc_dir',calc_dir) )
    db.dadd('dir',('out_dir',out_dir) )
    db.dadd('dir',('input_dir',input_dir) )
    db.dadd('dir',('work_dir',work_dir) )
    
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
            'nq2', 'nq3', 'electron_phonon','fildvscf', 'el_ph_sigma', 
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
    work_dir = db.dget('dir', 'work_dir')
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

    #make scf2 input QE file
    ase.io.write(filename = file_name2,
                images = structure,
                format = 'espresso-in',
                input_data = pw_parameters,
                pseudopotentials = pseudo,
                kpts = kgrid2,
                koffset = koffset2)    

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
    
    qe_dict = program + '_par'
   
    #get data from database
    db = pk.load(_DATABASE, False)

    prefix = db.dget('pw_par','prefix')
    parameters = dict(db.dgetall(qe_dict).items())
    section_name = db.dget('section_name', program)
    file_order = db.dget('file_order', program)
    file_name = db.dget('input_name', program)
    input_dir = db.dget('dir', 'input_dir')
    work_dir = db.dget('dir', 'work_dir')

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

    #return work directory    
    os.chdir(work_dir)  

    return 
##----------------------------------------------------------------------------##

def run_espresso_in(program):

    # Run espresso_in file
    # 
    #  @param: 'program': (str) identifier of the program input file to be 
    #                     run. Vslues possible: 'ph', 'matdyn', 'q2r'. 
    #  @return: no return.

    #Get info from database
    db = pk.load(_DATABASE, False)
    
    input_dir = db.dget('dir', 'input_dir')
    work_dir = db.dget('dir', 'work_dir')
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
    work_dir = db.dget('dir', 'work_dir')
    
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

user_input = sys.argv[1]
input_info = get_input_data(data = user_input)
dir_info   = mk_work_environment(input_info = input_info)

_DATABASE = mk_data_base(input_info = input_info, 
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
run_espresso_in(program= 'lambda')

print('\nTotal Execution time:', time.time()-start_time)
