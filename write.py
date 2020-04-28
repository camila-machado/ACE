import pickledb as pk

#importar dados anteriores
db = pk.load('config.db', False)

#------------------Dicionarios de Dados------------------#

#qe_par - entrada da funcao ase.write
db.dcreate ('pw_par')

db.dadd('pw_par',('restart_mode','from_scratch') )
db.dadd('pw_par',('prefix','H3S') )
db.dadd('pw_par',('pseudo_dir','./pseudo') )
db.dadd('pw_par',('out','./out') )
db.dadd('pw_par',('occupations','smearing') )
db.dadd('pw_par',('smearing','methfessel-paxton') )
db.dadd('pw_par',('degauss',0.03) )
db.dadd('pw_par',('ecutwfc',10.0) )
db.dadd('pw_par',('ecutrho',120.0) )

#pseudopotenciais
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

#grids
db.dcreate ('grids')

db.dadd('grids',('kcoarse_div',(10,10,10)) )
db.dadd('grids',('kcoarse_off',(0,0,0)) )
db.dadd('grids',('kdense_div',(10,10,10)) )
db.dadd('grids',('kdense_off',(0,0,0)) )

#estrutura.cif
db.dcreate ('structure')

###entrada do programa
db.dadd('structure',('dir','/home/ABTLUS/camila.araujo/Documents/Programa/cif_database/H3S-200GPa.cif') )
db.dadd('structure',('format','cif') )

#Escrever dados
db.dump()