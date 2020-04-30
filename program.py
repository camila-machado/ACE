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

#!/usr/bin/env python

#MODULES

import numpy as np
import ase
import ase.io
import pickledb as pk
import sys

import subprocess
import os
import shutil

from pathlib import PurePath

#DEFINES

_DATABASE = 'config.db'

#FUNCTIONS

def get_input_data(data):

    # Read the input from user and save the correspondent file information 
    # (file directory, file format, and prefix) in the database
    #  
    #  @param: - 'data': (str) input from user
    #  @return: no return. the function modify the database directly

    #get input information
    file_dir = os.path.normpath(data)
    path = PurePath(file_dir)
       
    prefix = path.stem
    file_format = path.suffix.strip('.')

    input_info = {
    'prefix': prefix,
    'file_dir':file_dir,
    'file_format':file_format
    }

    return input_info
##----------------------------------------------------------------------------##
def mk_work_enviroment(input_info):

    # 
    #  
    #  @param:  
    #  @return: 

    prefix = input_info.get('prefix')

    old_path = os.getcwd()
    print(old_path)

    new_path = os.path.join(old_path, prefix, 'calc_dir')
    print(new_path)

    exist = os.path.exists(new_path)
    if exist == True:
        shutil.rmtree(new_path)

    os.makedirs(new_path)
    os.chdir(new_path)

    return
##----------------------------------------------------------------------------##
def mk_data_base(input_info):
    
    # Construct the files structure to be used in 'write_espresso_in' function
    #
    #  
    #  @param:  
    #  @return:

    prefix = input_info.get('prefix')
    file_dir = input_info.get('file_dir')
    file_format = input_info.get('file_format')
    
    os.system('python3 ../../write.py')

    #create prefix dependent parameters
    #ph
    fildyn = prefix + '.dyn'
    fildvscf = prefix + '.dv'
    outdir = './out'

    #q2r
    flfrc = prefix + '.frc'

    #matdyn    
    flfrq= prefix+'.freq'    
    fldos= prefix+'.phonon.dos'

    #Programs' parameters
    db = pk.load(_DATABASE, False)

    db.dcreate ('cell_structure')
    db.dadd('cell_structure',('dir',file_dir) )
    db.dadd('cell_structure',('format',file_format) )

    db.dadd('pw_par',('prefix',prefix) )
    db.dadd('pw_par',('outdir',outdir) )
    db.dadd('pw_par',('pseudo_dir','../../pseudo') )

    db.dadd('ph_par',('prefix',prefix) )
    db.dadd('ph_par',('fildyn',fildyn) )
    db.dadd('ph_par',('fildvscf',fildvscf) )
    db.dadd('ph_par',('outdir',outdir) )

    db.dadd('q2r_par',('flfrc',flfrc) )
    db.dadd('q2r_par',('fildyn',fildyn) )

    db.dadd('matdyn_par',('flfrc',flfrc) )
    db.dadd('matdyn_par',('flfrq',flfrq) )
    db.dadd('matdyn_par',('fldos',fldos) )

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

    db.dump()

    return
##----------------------------------------------------------------------------##

def get_file (file_dir):

    # Open file and save it's content in a list of strings
    #  
    #  @param:  - 'file_dir': (str) directory of the file to be read
    #  @return: - 'content' : (str list) each line of the file is saved as a
    #                         element of 'content'         

    f = open(file_dir, 'r')
    content = f.readlines()
    f.close()

    return content
##----------------------------------------------------------------------------##

def get_word_position (word, content):

    # Get the str element posiotion of where 'word' was found in the 'content' 
    # list of strings. Only returns the first occurance.
    #  
    #  @param:  - 'word'    : (str) to be searched
    #           - 'content' : (str list) where 'word' will be searched
    #  @return: - 'position': (int) index position of the str element that 
    #                         contains 'word'
    #
                    
    for element in content:
        if element.find(word.upper() or word.lower()) > -1:
            position = content.index(element)
            break

    return position
##----------------------------------------------------------------------------##

def write_file (file_dir, content):

    # open/create file and rewrite it's content
    #  
    #  @param:  - 'file_dir': (str) directory of the file to be oppend/created
    #           - 'content' : (str list) list of strings that constains 
    #                         each line of the file to be rewritten
    #  @return: no return.
           
    f = open(file_dir, 'w')

    for line in content:
        f.write(line)

    f.close()

    return
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

    content = get_file(file_dir = file_dir)

    section_pos = get_word_position (word= section, 
                                    content = content)
    
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

    write_file (file_dir = file_dir, 
                        content = content)

    return
##----------------------------------------------------------------------------##

def remove_repeated (list_input):

    # Remove all the repeated values of a list. 
    # 
    # When the list is used as keys to create a dictionare, the duplicats 
    # are automatically removed because keys are unique. 
    # Later, the keys from the dictionary are transformed into a list again.
    #  
    #  @param: - 'list_input': (list) to be removed repetitions
    #  @return: - 'list_output': (list) without repetition

    list_output = list( dict.fromkeys(list_input) )

    return list_output
##----------------------------------------------------------------------------##

def set_pseudo_potencials (cell_structure):

    # Prepaire the pseudo potencials dir used to write th espresso-in scf 
    # input file. 
    #
    # - Read the elements from atoms object and save them in a list
    # - Remove elements list repetition 
    # - Get pseudopotentials directories for each cell element from database
    # - Group pseudopotencials of interest in a dictionary
    #  
    #  @param:  - 'cell_structure': (atoms obeject) that contains all elements 
    #  @return: - 'pseudo_dir'    : (dict) elements and correspondent 
    #                               directories of pseudopotencial files 
    #                               to be used to create pw (scf) input files
    
    elements = cell_structure.get_chemical_symbols()
    elements = remove_repeated(elements)

    db = pk.load(_DATABASE, False)
    
    pseudo_dir = {}
    for element in elements:
        pseudo = db.dget('pseudo',element)
        pseudo_dir.update({element:pseudo})

    db.dump()

    return pseudo_dir
##----------------------------------------------------------------------------##

def set_ecut (pseudo):

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

    pseudo_dir = pseudo.values()

    pseudo_content =[]
    for directory in pseudo_dir:
        file_name = '../../pseudo/' + directory
        file_content = get_file(file_name)
        pseudo_content.append(file_content)
    

    ecutwfc_values = []
    ecutrho_values = []
    for content in pseudo_content:
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


def save_file_name (program, name):

##------------------------------------------------------------------------------
#  
##------------------------------------------------------------------------------
    db = pk.load(_DATABASE, False)

    db.dadd('dir',(program, name))

    db.dump()
    
    return 

##----------------------------------------------------------------------------##

def write_espresso_in_pw (scf1):

    # Write input pw.x file for scf1 or scf2 calculus. If scf1 = True, it writs 
    # the scf1 file. Otherwise, it writes scf2 file
    #
    # - Read cell_structure user input file and save atoms object
    # - Set pseudopotencials acoordinly with cell_structure elements 
    # - Get ecut energies from pseudopotential files suggestions 
    # - Get the rest of parametrs from database accordinly with scf1 tag
    # - Write output file
    # 
    #  @param: - scf1: (bool) indicates if it wiol be genarated scf1 or 
    #                  scf2 file
    #  @return: no return.
    
    db = pk.load(_DATABASE, False)

    #get cell structure
    file_dir = db.dget('cell_structure','dir')
    cell_structure_format = db.dget('cell_structure','format')    
    structure = ase.io.read(filename = file_dir,
                            format=cell_structure_format,
                            subtrans_included = False,
                            primitive_cell= True)
    
    #set variables
    pseudo = set_pseudo_potencials(structure)
    ecutwfc, ecutrho = set_ecut(pseudo)
    db.dadd('pw_par',('ecutwfc',ecutwfc) )
    db.dadd('pw_par',('ecutrho',ecutrho) )
    pw_parameters = db.dgetall('pw_par')
    prefix = db.dget('pw_par','prefix')

    if(scf1 == True):
        kgrid = db.dget('grids','kdense_div')
        koffset = db.dget('grids','kdense_off')
        scf_file_name = prefix+'.scf1.in'
        db.dump()
        save_file_name(program = 'scf1', name = scf_file_name)
        
    else:
        kgrid = db.dget('grids','kcoarse_div')
        koffset = db.dget('grids','kcoarse_off')
        scf_file_name = prefix+'.scf2.in'
        db.dump()
        save_file_name(program = 'scf2', name = scf_file_name)

    #make scf input QE files
    ase.io.write(filename= scf_file_name,
                images = structure,
                format='espresso-in',
                input_data=pw_parameters,
                pseudopotentials=pseudo,
                kpts=kgrid,
                koffset=koffset)

    if(scf1 == True):
        add_parameter_qe(file_dir = scf_file_name, 
                        parameter = ('la2F','.true.'), 
                        section = 'system')

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

    db.dump()

    #create and write file
    file_name = prefix + '.' + program + '.in'
    save_file_name(program = program, name = file_name)

    line = ['&'+section_name.upper()+'\n/\n']
    write_file(file_dir = file_name, content = line)

    file_order.reverse()
    for key in file_order:
        value = str(parameters.get(key))
        add_parameter_qe(file_dir = file_name, 
                        parameter = (key,value), 
                        section = section_name) 

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

    dir_qe = str(db.dget('dir','qe_programs'))    
    file_name = str(db.dget('dir', program))
    prefix = (db.dget('pw_par', 'prefix'))
    np = str(db.dget('mpi','np'))  
    nk = str(db.dget('mpi','nk'))

    db.dump()    

    #set variables
    output_name = prefix + '.' + program + '.out'
    
    if program.find('scf') > -1:
        dir_program = dir_qe + 'pw' +'.x'
    else:
        dir_program = dir_qe + program +'.x'

    #make command str
    if program.find('lambda' or 'q2r') > -1:
        command = dir_program +' < ' + file_name + ' > ' + output_name
    elif program.find('matdyn') > -1:
        command = ('mpiexec' + ' -np ' + np + ' ' + dir_program + ' -in ' + file_name 
                   + ' > ' + output_name)
    else:
        command = ('mpiexec' + ' -np ' + np + ' ' + dir_program + ' -nk ' + nk + ' -in ' 
                  + file_name + ' > ' + output_name)

    print (command)

    #run program
    text1 = '\n running {program_txt}... \n'
    print(text1.format(program_txt = program))
    os.system(command) 
    text2 = '\n {program_txt} fineshed! :D \n'
    print(text2.format(program_txt = program))

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
    mu =  db.dget('lambda_par','mu')
    file_dyn_prefix = db.dget('ph_par','fildyn')
    prefix = db.dget('pw_par','prefix')

    db.dump()

    #get q points information
    #file .dyn0 contains q points information used lambda.in file
    file_name = file_dyn_prefix + '0'
    content = get_file(file_dir = file_name)
    content.pop(0)

    q_points_info = content 
    num_q_points = int(q_points_info[0])

    #read/get information from dyn files
    dyn_info = []
    for num in range(1,num_q_points+1):
        file_name = file_dyn_prefix + str(num)
        content = get_file(file_dir = file_name)
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

    lambda_file_name = prefix + '.lambda.in'
    save_file_name(program = 'lambda',name = lambda_file_name)

    write_file (file_dir = lambda_file_name, content=lambda_info)

    print('\n lambda.in file wrote \n')

    return
##----------------------------------------------------------------------------##

################################################################################
##----------------------------------------------------------------------------##
################################################################################

#treat input data and prepaire 
user_input = sys.argv[1]

input_info = get_input_data(data = user_input)
mk_work_enviroment(input_info = input_info)
mk_data_base(input_info = input_info)

#write input files
write_espresso_in_pw(scf1 = True)
write_espresso_in_pw(scf1= False)
write_espresso_in(program = 'ph')
write_espresso_in(program = 'q2r')
write_espresso_in(program = 'matdyn')

#run programs pw, ph, q2r and matdyn
run_espresso_in(program= 'scf1')
run_espresso_in(program= 'scf2')
run_espresso_in(program= 'ph')
run_espresso_in(program= 'q2r')
run_espresso_in(program= 'matdyn')

write_espresso_in_lambda()

run_espresso_in(program= 'lambda')