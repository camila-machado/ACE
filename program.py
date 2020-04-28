
import pickledb as pk

#importar dados anteriores
db = pk.load('config.db', False)

#-----Variáveis do ase.io.read
structure_dir = db.dget('structure','dir')
structure_format = db.dget('structure','format')

#-----Variáveis do ase.io.write (4)
kgrid = db.dget('grids','kdense_div')
koffset = db.dget('grids','kdense_off')

#!!!!!!!!!! precis de função que leia os pw especificos
pseudo = db.dgetall('pseudo')

pw_par = db.dgetall('pw_par')
