##------------------------------------------------------------------------------
# CNPEM - Centro Nacional de Pesquisa em Energia e Materiais, Brazil
# Sirius - Research Group EMA
#
# Code project: TC_caculus
# 
# Objectif:
# Automate the use of Quantum Espresso (QE) for the calculus of critical 
# superconductivity temperature of diferent molecules and cell structures
# 
##------------------------------------------------------------------------------

#MODULES

import numpy as np
import ase
import ase.io
import pickledb as pk
import sys

import subprocess
import os

#DEFINES

_DATABASE = 'config.db'

#FUNCTIONS

def save_cell_structure (cell_structure_dir):

##------------------------------------------------------------------------------
#  Read the input line and save the correspondent file data (file directory, 
# file format, and prefix) in the config database
#  
#  @param: - 'cell_structure_dir': directory of the cell structure file
#  @return: no return. the function modify the database directly
#
##------------------------------------------------------------------------------

    cell_structure_dir = cell_structure_dir.strip()
    
    #get entrance file format
    aux = cell_structure_dir.split("/")
    aux_len = len(aux)
    aux2 = aux[aux_len-1].split(".")
    aux2_len = len(aux2)
    file_format = aux2[aux2_len - 1]

    #set prefix name of future output from input file name
    prefix =''
    for i in range(0,aux2_len-1):
        prefix += aux2[i]
                                                                                    
    #check-----------------------
    print('\nPrefix:\n'+prefix)
    print('\nIput file directory:\n'+cell_structure_dir)
    print('\nInput file format:\n'+file_format)
    #----------------------------

    #save data obtaind in database
    db = pk.load(_DATABASE, False)

    db.dcreate ('cell_structure')

    db.dadd('cell_structure',('dir',cell_structure_dir) )
    db.dadd('cell_structure',('format',file_format) )
    db.dadd('pw_par',('prefix',prefix) )

    db.dump()
    return


##----------------------------------------------------------------------------##

def get_file_content (file_input):

##------------------------------------------------------------------------------
# open file and save content
#           
##------------------------------------------------------------------------------

    file_in = open(file_input, 'r')
    content = file_in.readlines()
    file_in.close()

    return content

##----------------------------------------------------------------------------##

def get_str_list_position (str_input, list_input):

##------------------------------------------------------------------------------
#get 'str' line posiotion in the file
#           
##------------------------------------------------------------------------------
    position = 0
    for x in list_input:
        if x.find(str_input.upper() or str_input.lower()) > -1:
            break
        position +=1

    return position
##----------------------------------------------------------------------------##

def write_file_content (file_input, content):

##------------------------------------------------------------------------------
#rewrite content in file
#           
##------------------------------------------------------------------------------
    file_in = open(file_input, 'w')

    for line in content:
        file_in.write(line)

    file_in.close()

    return

##----------------------------------------------------------------------------##    

def add_parameter_qe (file_input, parameter, section):

##------------------------------------------------------------------------------
#  Add a parameter to a input Quantum Espresso file at the indicated section.
#
# First, the file content is saved in a list, second, the position of the 
# section in the list that represents the file is found. 
# The idented new information is inserted in the list and corect section and the
# content is rewritten in the file.
#  
#  @param: - 'file_input': directory of input file to be modified
#          - 'parameter':  tuple of key, value to be added
#          - 'section':    section name of the QE input file where the parameter 
#                          should be added
#  @return: no return, the function modify and save the document 
#           
##------------------------------------------------------------------------------

    content = get_file_content (file_input = file_input)

    section_pos = get_str_list_position (str_input= section, 
                                              list_input= content)
        
    #count spaces to write with beautiful QE identation
    aux = len(parameter[0])
    space = ''
    
    i = 0
    while i < 17-aux:
        space += ' '
        i += 1

    #transform parameter in str
    parameter0 = str(parameter[0])
    parameter1 = str(parameter[1])

    new_position = section_pos+1
    new_info = '   '+parameter0+space+'= '+parameter1+ '\n'
    
    new_content = content
    new_content.insert(new_position,new_info)

    write_file_content (file_input = file_input, content = new_content)

    return

##----------------------------------------------------------------------------##

def remove_repeated (list_input):

##------------------------------------------------------------------------------
#  Remove all the repeated values of a list. When the list is used as keys to 
# create a dictionare, the duplicats are automatically removed because keys as 
# unique. Later, the keys from teh dictionary are transformed into a list again.
#  
#  @param: - 'list_input': variable list name to be removed repetiions
#  @return: - 'list_output': variable list without repetition
##------------------------------------------------------------------------------
    list_output = list( dict.fromkeys(list_input) )
    return list_output

##----------------------------------------------------------------------------##

def set_pseudo_potencials (molecule_structure):

##------------------------------------------------------------------------------
#  Prepaire the pseudo potencials dir used to write th espresso-in scf input 
# file. 
# - Read the atoms object elements
# - Remove repetition 
# - Open database 
# - Search in databasa for the pseudopotentials directory of each element 
# - Group the pseudopotencials of interest in a dictionary
#  
#  @param: - 'molecule_structure': atoms obeject that contains all elements 
#  @return: - 'pseudo_dir': list of the directory of pseudopotencial files
#                                 to be used to describe the molecule structure
##------------------------------------------------------------------------------
    
    #getting elements from cell structure input file
    elements = molecule_structure.get_chemical_symbols()
    elements = remove_repeated(elements)
    #check----------------------------------------
    #print('\nElements from cell structure file:')
    #print(elements)
    #---------------------------------------------

    db = pk.load(_DATABASE, False)
    
    #setting pseudo potencials files to be used
    pseudo_dir = {}
    for x in elements:
        pseudo_dir.update({x:db.dget('pseudo',x)})
    #check----------------------------------------
    #print('\nPseudo potential files:')
    #print(pseudo)
    #---------------------------------------------  

    db.dump()

    return pseudo_dir

##----------------------------------------------------------------------------##

def get_pseudo_ecut (pseudo):

##------------------------------------------------------------------------------
#  
##------------------------------------------------------------------------------
    pseudo_dir = pseudo.values()
    
    pseudo_content =[]
    for dir_element in pseudo_dir:
        pseudo_file_name = 'pseudo/' + dir_element
        file_content = get_file_content (pseudo_file_name)
        pseudo_content.append(file_content)
    

    ecutwfc_values = []
    ecutrho_values = []
    for content in pseudo_content:
        for line in content:
            if line.find('Suggested minimum cutoff for wavefunctions') > -1:
                aux = line.split(':')
                aux2 = aux[1].strip(' .Ry\n')
                ecutwfc_values.append(int(aux2))
            if line.find('Suggested minimum cutoff for charge density') > -1:
                aux = line.split(':')
                aux2 = aux[1].strip(' .Ry\n')
                ecutrho_values.append(int(aux2))
        
    ecutwfc = max(ecutwfc_values)
    ecutrho = max(ecutrho_values)
    
    return ecutwfc, ecutrho

##----------------------------------------------------------------------------##

def write_espresso_in_scf_file (scf1):

##------------------------------------------------------------------------------
#  
#If scf1 = True, it writs the first scf1 file. Otherwise, it writes teh second
# scf input file
# 
#  @param: - 
#  @return: - 
##------------------------------------------------------------------------------
    
    #-----import configuration data
    db = pk.load(_DATABASE, False)

    #-----variables for ase.io.read
    cell_structure_dir = db.dget('cell_structure','dir')
    cell_structure_format = db.dget('cell_structure','format')    

    #-----import cell structure from input file
    structure = ase.io.read(filename = cell_structure_dir,
                            format=cell_structure_format,
                            subtrans_included = False,
                            primitive_cell= True)
    #check----------------------------------------
    #print('\nCell structure obtained:')
    #print(structure.get_cell_lengths_and_angles())
    #---------------------------------------------

    #-----variables for ase.io.write 
    prefix = db.dget('pw_par','prefix')

    if(scf1 == True):
        kgrid = db.dget('grids','kdense_div')
        koffset = db.dget('grids','kdense_off')
        scf_file_name = prefix+'.scf1.in'
        
    else:
        kgrid = db.dget('grids','kcoarse_div')
        koffset = db.dget('grids','kcoarse_off')
        scf_file_name = prefix+'.scf2.in'

    pseudo = set_pseudo_potencials(structure)
    ecut = get_pseudo_ecut(pseudo)

    db.dadd('pw_par',('ecutwfc',ecut[0]) )
    db.dadd('pw_par',('ecutrho',ecut[1]) )

    pw_parameters = db.dgetall('pw_par')

    db.dump()

    #-----exporting scf input QE files
    ase.io.write(filename= scf_file_name,
                images = structure,
                format='espresso-in',
                input_data=pw_parameters,
                pseudopotentials=pseudo,
                kpts=kgrid,
                koffset=koffset)

    if(scf1 == True):
        add_parameter_qe(file_input = scf_file_name, 
                    parameter = ('la2F','.true.'), 
                    section = 'system')

    return 

##----------------------------------------------------------------------------##

def write_espresso_in_file(qe_dict, file_order, file_sufix, section_name):

##------------------------------------------------------------------------------
# Write espresso_in file
# 
#  @param: - qe_par - str of dictionary parameters name
#          - file_order - list with keys file order 
#          - file_sufix - str accordinly defined by the type of file generated
#          - section_name - str with the final file section name 
#  @return: - 
##------------------------------------------------------------------------------
    
    #save/make variables 
    db = pk.load(_DATABASE, False)

    prefix = db.dget('pw_par','prefix')
    parameters = dict(db.dgetall(qe_dict).items())

    db.dump()

    #Define file order
    file_order.reverse()

    #create and write file
    file_name = prefix+file_sufix
    qe_file = open(file_name, 'w')
    qe_file.write('&'+section_name.upper()+'\n/\n')
    qe_file.close()

    for key in file_order:
        key_value = str(parameters.get(key))
        add_parameter_qe(file_input = file_name, 
                        parameter = (key,key_value), 
                        section = section_name) 

    return 
##----------------------------------------------------------------------------##

def write_espresso_in_ph():

##------------------------------------------------------------------------------
# Write espresso_in el-phon file
# 
#  @param: - 
#  @return: - 
##------------------------------------------------------------------------------
    
    #save/make ph variables 
    db = pk.load(_DATABASE, False)

    #pw dependent variables
    prefix = db.dget('pw_par','prefix')
    outdir = db.dget('pw_par','outdir')
    
    fildyn = "'" +prefix + '.dyn'+ "'"
    fildvscf = "'" +prefix + '.dv'+ "'"
    prefix = "'"+prefix+"'"
    outdir = "'"+outdir+"'"

    db.dadd('ph_par',('prefix',prefix) )
    db.dadd('ph_par',('outdir',outdir) )
    db.dadd('ph_par',('fildyn',fildyn) )
    db.dadd('ph_par',('fildvscf',fildvscf) )

    db.dump()

    ph_file_order = ['prefix', 'outdir', 'tr2_ph', 'fildyn', 
                        'ldisp', 'nq1', 'nq2', 'nq3', 'electron_phonon', 
                        'fildvscf', 'el_ph_sigma', 'el_ph_nsigma']

    write_espresso_in_file(qe_dict='ph_par', 
                           file_order = ph_file_order, 
                           file_sufix = '.ph.in', 
                           section_name = 'inputph')

    return 
##----------------------------------------------------------------------------##

def write_espresso_in_q2r():

##------------------------------------------------------------------------------
# Write espresso_in q2r file
# 
#  @param: - 
#  @return: - 
##------------------------------------------------------------------------------
    
    #save/make ph variables 
    db = pk.load(_DATABASE, False)

    #pw dependent variables
    prefix = db.dget('pw_par','prefix')
    fildyn =  db.dget('ph_par','fildyn')
    
    flfrc = "'"+prefix + '.frc'+"'" #acrescentar rede de q
    prefix = "'"+prefix+"'"

    db.dadd('q2r_par',('flfrc',flfrc) )
    db.dadd('q2r_par',('fildyn',fildyn) )

    db.dump()

    #Define file order
    q2r_file_order = ['zasr', 'fildyn', 'flfrc', 'la2F']
 
    write_espresso_in_file(qe_dict='q2r_par', 
                           file_order = q2r_file_order, 
                           file_sufix = '.q2r.in', 
                           section_name = 'input')
    return 
##----------------------------------------------------------------------------##

def write_espresso_in_matdyn():

##------------------------------------------------------------------------------
# Write espresso_in matdyn file
# 
#  @param: - 
#  @return: - 
##------------------------------------------------------------------------------
    
    #save/make ph variables 
    db = pk.load(_DATABASE, False)

    #pw dependent variables
    prefix = db.dget('pw_par','prefix')
    flfrc = db.dget('q2r_par','flfrc')
      
    flfrq= "'"+prefix+'.freq'+"'"      
    fldos= "'"+prefix+'.phonon.dos'+"'"

    db.dadd('matdyn_par',('flfrc',flfrc) )
    db.dadd('matdyn_par',('flfrq',flfrq) )
    db.dadd('matdyn_par',('fldos',fldos) )
    
    db.dump()

    #Define file order
    matdyn_file_order = ['asr', 'flfrc', 'flfrq', 'la2F', 'dos', 'fldos', 'nk1', 'nk2',
                        'nk3','ndos']
    
    write_espresso_in_file(qe_dict='matdyn_par', 
                           file_order = matdyn_file_order, 
                           file_sufix = '.matdyn.in', 
                           section_name = 'input')
    return 
##----------------------------------------------------------------------------##

def run_espresso_in(dir_qe, program, dir_file, output_name ):

##------------------------------------------------------------------------------
# Run espresso_in file
# 
#  @param: - 
#  @return: - 
##------------------------------------------------------------------------------
    
    
    command = [dir_qe + program,'<',dir_file,'>',output_name]
    command_str = command[0]+' '+command[1]+' '+command[2]+' '+command[3]+' '+command[4]+' ' 
    print(command_str)
    os.system(command_str)

    return 
##----------------------------------------------------------------------------##
def get_dyn_file_info(dyn_file):

##------------------------------------------------------------------------------
# Get multiplicity and upper freq
# 
#  @param: - 
#  @return: - 
##------------------------------------------------------------------------------
    
    freq_info = []
    multiplicity = 0

    for line in dyn_file:
        if line.find('Dynamical  Matrix in cartesian axes') > -1:
            multiplicity += 1
        if line.find('freq') > -1:
            aux = line.split('=')
            aux2 = aux[1].split(' ')
            freq_position = 0
            for i in aux2:
                if i.find('THz') > -1:
                    freq_position = aux2.index(i) - 1
                    freq_info.append(float(aux2[freq_position]))
    
    max_freq = max(freq_info)

    return max_freq, multiplicity
##----------------------------------------------------------------------------##

def write_espresso_in_lambda():

##------------------------------------------------------------------------------
# Write espresso_in lambda file
# 
#  @param: - 
#  @return: - 
##------------------------------------------------------------------------------
    
    #get/make lambda variables 
    db = pk.load(_DATABASE, False)

    sigma = db.dget('lambda_par','sigma_omega')
    mu =  db.dget('lambda_par','mu')
    file_dyn_prefix = db.dget('ph_par','fildyn').replace("'",'')
    prefix = db.dget('pw_par','prefix')

    db.dump()

    #get q points
    file_name = file_dyn_prefix + '0'
    content = get_file_content(file_input = file_name)
    content.pop(0)

    q_points_info = content 
    num_q_points = int(q_points_info[0])

    #read .dyn files
    dyn_info = []
    for num in range(1,num_q_points+1):
        file_in = file_dyn_prefix + str(num)
        freq_file = get_file_content(file_input = file_in)
        dyn_info.append(freq_file)

    #filter freq information from freq file
    multiplicity = []
    freq = []
    for element in dyn_info:
        info = get_dyn_file_info(dyn_file = element)
        freq.append(info[0])
        multiplicity.append(info[1])

    upper_freq = int((float(max(freq))+1))

    part_1 = str(upper_freq)+'  '+str(sigma)+'  1\n'
    
    part_2 = q_points_info
    for i in range(1,num_q_points+1):
        aux = '  '+str(multiplicity[i-1]) +'\n'
        part_2[i] = part_2[i].replace('\n', aux)
    
    part_3 = []
    aux = 'elph_dir/elph.inp_lambda.'
    for j in range(1,num_q_points+1):
        line = aux + str(j) + '\n'
        part_3.append(line)
    part_3.append(str(mu)+'\n')

    lambda_info =[]
    lambda_info.append(part_1)
    lambda_info.extend(part_2)
    lambda_info.extend(part_3)

    lambda_file_name = prefix + '.lambda.in'

    write_file_content (file_input = lambda_file_name, content=lambda_info)

    return
##----------------------------------------------------------------------------##



################################################################################
##----------------------------------------------------------------------------##
################################################################################

#treat input data
file_in = sys.argv[1]
save_cell_structure(cell_structure_dir = file_in)

#write espresso_in files
write_espresso_in_scf_file(scf1 = True)
write_espresso_in_scf_file(scf1= False)
write_espresso_in_ph()
write_espresso_in_q2r()
write_espresso_in_matdyn()

#run programs pw, ph, q2r and matdyn

print('\n running scf1... \n')

run_espresso_in(dir_qe= '/home/ABTLUS/camila.araujo/Downloads/Programs/qe-6.5/bin/', 
                program= 'pw.x', 
                dir_file= '/home/ABTLUS/camila.araujo/Documents/Programa/H3S-200GPa.scf1.in', 
                output_name= 'scf1.out')

print('\n scf1 fineshed! :D \n')

print('\n running scf2... \n')
run_espresso_in(dir_qe= '/home/ABTLUS/camila.araujo/Downloads/Programs/qe-6.5/bin/', 
                program= 'pw.x', 
                dir_file= '/home/ABTLUS/camila.araujo/Documents/Programa/H3S-200GPa.scf2.in', 
                output_name= 'scf2.out')
print('\n scf2 fineshed! :D \n')


print('\n running phonons... \n')
#run_espresso_in(dir_qe= '/home/ABTLUS/camila.araujo/Downloads/Programs/qe-6.5/bin/', 
#                program= 'ph.x', 
#                dir_file= '/home/ABTLUS/camila.araujo/Documents/Programa/H3S-200GPa.ph.in', 
#                output_name= 'ph.out')
print('\n phonons fineshed! :D \n')



print('\n running q2r... \n')
run_espresso_in(dir_qe= '/home/ABTLUS/camila.araujo/Downloads/Programs/qe-6.5/bin/', 
                program= 'q2r.x', 
                dir_file= '/home/ABTLUS/camila.araujo/Documents/Programa/H3S-200GPa.q2r.in', 
                output_name= 'q2r.out')
print('\n q2r fineshed! :D \n')


print('\n running matdyn... \n')
run_espresso_in(dir_qe= '/home/ABTLUS/camila.araujo/Downloads/Programs/qe-6.5/bin/', 
                program= 'matdyn.x', 
                dir_file= '/home/ABTLUS/camila.araujo/Documents/Programa/H3S-200GPa.matdyn.in', 
                output_name= 'matdyn.out')
print('\n matdyn fineshed! :D \n')


write_espresso_in_lambda()
print('\n lambda.in file wrote \n')

print('\n running lambda... almost there \n')
run_espresso_in(dir_qe= '/home/ABTLUS/camila.araujo/Downloads/Programs/qe-6.5/bin/', 
                program= 'lambda.x', 
                dir_file= '/home/ABTLUS/camila.araujo/Documents/Programa/H3S-200GPa.lambda.in', 
                output_name= 'lambda.out')
print("\n lambda fineshed! :D \n Go get your Tc's --->")









