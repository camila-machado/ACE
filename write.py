import pickledb as pk

#importar dados anteriores
db = pk.load('config.db', False)

#------------------Dicionarios de Dados------------------#
##------------------------------------------------------------------------------

#dir - 
db.dcreate ('dir')

db.dadd('dir',('qe_programs','/home/ABTLUS/camila.araujo/Downloads/Programs/qe-6.5/bin/') )
db.dadd('dir',('input','/home/ABTLUS/camila.araujo/Documents/Programa/') )

#mpi - 
db.dcreate ('mpi')

db.dadd('mpi',('np',4) )
db.dadd('mpi',('nk',4) )

#qe_pw_par - entrada da funcao ase.write
db.dcreate ('pw_par')

db.dadd('pw_par',('restart_mode','from_scratch') )
db.dadd('pw_par',('pseudo_dir','./pseudo') )
db.dadd('pw_par',('outdir','./out') )
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

db.dadd('ph_par',('outdir','./out') )
db.dadd('ph_par',('tr2_ph','1e-12') )
db.dadd('ph_par',('ldisp','.true.') )
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

#pseudo potenciais
db.dcreate ('pseudo')

db.dadd('pseudo',('Ba','Ba.pbe-spn-kjpaw_psl.1.0.0.UPF') )
db.dadd('pseudo',('Bi','Bi.pbe-dn-kjpaw_psl.1.0.0.UPF') )
db.dadd('pseudo',('Ca','Ca.pbe-spn-kjpaw_psl.1.0.0.UPF') )
db.dadd('pseudo',('Cu','Cu.pbe-dn-kjpaw_psl.1.0.0.UPF') )
db.dadd('pseudo',('H' ,'H.pbe-kjpaw_psl.1.0.0.UPF') )
db.dadd('pseudo',('Hg','Hg.pbe-n-kjpaw_psl.1.0.0.UPF') )
db.dadd('pseudo',('O' ,'O.pbe-n-kjpaw_psl.1.0.0.UPF') )
db.dadd('pseudo',('Pd','Pd.pbe-n-kjpaw_psl.1.0.0.UPF') )
db.dadd('pseudo',('Pr','Pr.pbe-spdn-kjpaw_psl.1.0.0.UPF') )
db.dadd('pseudo',('S' ,'S.pbe-n-kjpaw_psl.1.0.0.UPF') )
db.dadd('pseudo',('Sr','Sr.pbe-spn-kjpaw_psl.1.0.0.UPF') )
db.dadd('pseudo',('Y' ,'Y.pbe-spn-kjpaw_psl.1.0.0.UPF') )

#Escrever dados
db.dump()