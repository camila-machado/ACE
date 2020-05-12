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
import pickledb as pk
import sys

#user input
_DATABASE = sys.argv[1]

#create database
db = pk.load(_DATABASE, False)

#------------------Dicionarios de Dados------------------#
##------------------------------------------------------------------------------

#opcional: caso os programas do quantum-espresso não estejam no $PATH do Linux
#dir - 
db.dcreate ('dir')
db.dadd('dir',('qe_programs','') )
db.dadd('dir',('pseudo_folder','/home/camila/Documentos/EMA/Program-TC/pseudo') )

#mpi - 
db.dcreate ('mpi')

db.dadd('mpi',('np',4) )
db.dadd('mpi',('nk',4) )

#qe_pw_par - entrada da funcao ase.write
db.dcreate ('pw_par')

db.dadd('pw_par',('restart_mode','from_scratch') )
db.dadd('pw_par',('occupations','smearing') )
db.dadd('pw_par',('smearing','marzari-vanderbilt') )
db.dadd('pw_par',('degauss',0.05) )
db.dadd('pw_par',('conv_thr',1e-10) )

#grids
db.dcreate ('grids')

db.dadd('grids',('kcoarse_div',(9,9,9)) )
db.dadd('grids',('kcoarse_off',(0,0,0)) )
db.dadd('grids',('kdense_div',(18,18,18)) )
db.dadd('grids',('kdense_off',(0,0,0)) )

#qe_ph_par
db.dcreate ('ph_par')

db.dadd('ph_par',('tr2_ph',1e-12) )
db.dadd('ph_par',('ldisp','.true.') )
db.dadd('ph_par',('recover','.false.') )
db.dadd('ph_par',('nq1',3) )
db.dadd('ph_par',('nq2',3) )
db.dadd('ph_par',('nq3',3) )
db.dadd('ph_par',('electron_phonon','interpolated') )
db.dadd('ph_par',('el_ph_sigma',0.005) )
db.dadd('ph_par',('el_ph_nsigma',10) )

#qe_q2r_par
db.dcreate ('q2r_par')

db.dadd('q2r_par',('zasr','simple') )
db.dadd('q2r_par',('la2F','.true.') )

#qe_matdyn_par
db.dcreate ('matdyn_par')

db.dadd('matdyn_par',('asr','simple') )
db.dadd('matdyn_par',('la2F','.true.') )
db.dadd('matdyn_par',('dos','.true.') )
db.dadd('matdyn_par',('nk1',10) )
db.dadd('matdyn_par',('nk2',10) )
db.dadd('matdyn_par',('nk3',10) )
db.dadd('matdyn_par',('ndos',50) )

#qe_lambda_par
db.dcreate ('lambda_par')

db.dadd('lambda_par',('sigma_omega', 0.12 ) )
db.dadd('lambda_par',('mu',0.16) )

#Escrever dados
db.dump()
