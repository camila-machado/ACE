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

This program creates a database .db file
spatially designed to provide all information needed by 
quantum espresso programs to calculate the critical 
superocnductivity temperature of a material
"""

import pickledb as pk
import sys

#database name is given as input
_DATABASE = sys.argv[1]

db = pk.load(_DATABASE, False)

#optional: in case quantum-espresso programs are not in linux $PATH
#qe - quantum espresso programs' path
db.dcreate ('qe')
db.dadd('qe',('qe_programs','') )

#pseudo - pseudo potentials' folder path
db.dcreate('pseudo')
db.dadd('pseudo',('pseudo_folder','/home/camila/Documentos/EMA/Program-TC/pseudo') )

#mpi - parallen running description
db.dcreate ('mpi')

db.dadd('mpi',('np',4) )
db.dadd('mpi',('nk',4) )

#grids - k-points grids used in scf calculations
db.dcreate ('grids')

db.dadd('grids',('kdistance', 0.2) ) #distance in angstrons of two consecutive grid k-points
db.dadd('grids',('qdistance', 0.6) ) #distance in angstrons of two consecutive grid q-points

#grids - grids used in scf and phonons calculations
db.dcreate('grid')
db.dadd('grid',('coarse_div',(9,9,9)) )
db.dadd('grid',('coarse_off',(0,0,0)) )
db.dadd('grid',('dense_div',(18,18,18)) )
db.dadd('grid',('dense_off',(0,0,0)) )
db.dadd('grid',('qpoints_div',(3,3,3)) )
db.dadd('grid',('qpoints_off',(0,0,0)) )

#Quantum-espresso programs' parameters:
db.dcreate ('pw_par')

db.dadd('pw_par',('restart_mode','from_scratch') )
db.dadd('pw_par',('occupations','smearing') )
db.dadd('pw_par',('smearing','marzari-vanderbilt') )
db.dadd('pw_par',('degauss',0.05) )
db.dadd('pw_par',('conv_thr',1e-10) )

db.dcreate ('ph_par')

db.dadd('ph_par',('tr2_ph',1e-12) )
db.dadd('ph_par',('ldisp','.true.') )
db.dadd('ph_par',('recover','.false.') )
db.dadd('ph_par',('electron_phonon','interpolated') )
db.dadd('ph_par',('el_ph_sigma',0.005) )
db.dadd('ph_par',('el_ph_nsigma',10) )

db.dcreate ('q2r_par')

db.dadd('q2r_par',('zasr','simple') )
db.dadd('q2r_par',('la2F','.true.') )

db.dcreate ('matdyn_par')

db.dadd('matdyn_par',('asr','simple') )
db.dadd('matdyn_par',('la2F','.true.') )
db.dadd('matdyn_par',('dos','.true.') )
db.dadd('matdyn_par',('nk1',10) )
db.dadd('matdyn_par',('nk2',10) )
db.dadd('matdyn_par',('nk3',10) )
db.dadd('matdyn_par',('ndos',50) )

db.dcreate ('lambda_par')

db.dadd('lambda_par',('sigma_omega', 0.12 ) )
db.dadd('lambda_par',('mu',0.16) )

db.dump()
