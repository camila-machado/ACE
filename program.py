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

    db.dcreate ('structure')

    db.dadd('structure',('dir',cell_structure_dir) )
    db.dadd('structure',('format',file_format) )
    db.dadd('pw_par',('prefix',prefix) )

    db.dump()
    return


##----------------------------------------------------------------------------##

def add_parameter (file_input, parameter, section):

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
    #open file and save content
    file_in = open(file_input, 'r')
    content = file_in.readlines()
    file_in.close()

    #get 'section' posiotion in the file
    section_position = 0
    for line in content:
        if line.find(section.upper() or section.lower()) > -1:
            break
        section_position +=1
    
    #count spaces to write with correct identation
    aux = len(parameter[0])
    space = ''
    
    i = 0
    while i < 17-aux:
        space += ' '
        i += 1

    #transform parameter in str
    parameter0 = str(parameter[0])
    parameter1 = str(parameter[1])

    content.insert(section_position+1, '   '+parameter0+space+'= '+parameter1+ '\n')


    #rewrite content in file
    file_in = open(file_input, 'w')

    for line in content:
        file_in.write(line)

    file_in.close()
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
    cell_structure_dir = db.dget('structure','dir')
    cell_structure_format = db.dget('structure','format')
    #-----variables for ase.io.write 
    pw_parameters = db.dgetall('pw_par')
    prefix = db.dget('pw_par','prefix')

    if(scf1 == True):
        kgrid = db.dget('grids','kdense_div')
        koffset = db.dget('grids','kdense_off')
        scf_file_name = prefix+'.scf1.in'
        
    else:
        kgrid = db.dget('grids','kcoarse_div')
        koffset = db.dget('grids','kcoarse_off')
        scf_file_name = prefix+'.scf2.in'
        
    db.dump()

    #-----import cell structure from input file
    structure = ase.io.read(filename = cell_structure_dir,
                            format=cell_structure_format,
                            subtrans_included = False,
                            primitive_cell= True)
    #check----------------------------------------
    #print('\nCell structure obtained:')
    #print(structure.get_cell_lengths_and_angles())
    #---------------------------------------------

    pseudo = set_pseudo_potencials(structure)

    #-----exporting scf input QE files
    ase.io.write(filename= scf_file_name,
                images = structure,
                format='espresso-in',
                input_data=pw_parameters,
                pseudopotentials=pseudo,
                kpts=kgrid,
                koffset=koffset)

    if(scf1 == True):
        add_parameter(file_input = scf_file_name, 
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
    qe_file.write('&'+section_name.upper()+'\n/')
    qe_file.close()

    for key in file_order:
        key_value = parameters.get(key)
        add_parameter(file_input = file_name, 
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
    
    fildyn = prefix + '.dyn'
    fildvscf = prefix + '.dv'
    db.dadd('ph_par',('prefix',prefix) )
    db.dadd('ph_par',('outdir',outdir) )
    db.dadd('ph_par',('fildyn',fildyn) )
    db.dadd('ph_par',('fildvscf',fildvscf) )

    db.dump()

    ph_file_order = ['prefix', 'outdir', 'trans', 'tr2_ph', 'fildyn', 
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
    
    flfrc = prefix + '.frc' #acrescentar rede de q
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
    flfrq= prefix+'.freq'      
    fldos= prefix+'.phonon.dos'

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



