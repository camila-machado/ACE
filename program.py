
import numpy as np
import ase
import ase.io
import pickledb as pk

#leitura e tratamento de dados

file_in = input("Enter the .cif cell structure file directory:")
file_in = file_in.strip()

dir_aux = file_in.split("/")
dir_aux_s = len(dir_aux)

dir_in = dir_aux[dir_aux_s-1].split(".")
dir_in_s = len(dir_in)

format_in = dir_in[dir_in_s - 1]

prefix_in =''
for i in range(0,dir_in_s-1):
    prefix_in += dir_in[i]

#checagem--------------------
print(prefix_in)
print(file_in)
print(format_in)
#----------------------------

db = pk.load('config.db', False)

#estrutura.cif
db.dcreate ('structure')

###entrada do programa
db.dadd('structure',('dir',file_in) )
db.dadd('structure',('format',format_in) )
db.dadd('pw_par',('prefix',prefix_in) )

#-----Variáveis do ase.io.read
structure_dir = db.dget('structure','dir')
structure_format = db.dget('structure','format')

#-----Variáveis do ase.io.write (4)
kgrid = db.dget('grids','kdense_div')
koffset = db.dget('grids','kdense_off')

#!!!!!!!!!! precis de função que leia os pw especificos
pseudo = db.dgetall('pseudo')

pw_par = db.dgetall('pw_par')

db.dump()

