# Copyright 2020, Camila Machado de Araújo
# Copyright 2021, Lucas Henrique Francisco 
# (see accompanying license files for details).

import os
import pickledb as pk

import pseudo

class IMakeDatabase:

    def make_database(self, name):

        dir_original = os.getcwd()
        dir_database = os.path.join(dir_original, name)

        db = self._make_db(name = name)
        db.dump()

        return dir_database
    
    def _make_db(self, name):
        pass

class TcDatabase (IMakeDatabase):

    def _make_db(self, name):

        _DATABASE = name

        db = pk.load(_DATABASE, False)

        #optional: in case quantum-espresso programs are not in linux $PATH
        #qe - quantum espresso programs' path
        db.dcreate ('calculation')
        db.dadd('calculation',('routine','TC') )

        #pseudo - pseudo potentials' folder path
        db.dcreate('pseudo')
        db.dadd('pseudo',('pseudo_folder', pseudo.DIR+'/USPP') )

        #mpi - parallen running description
        db.dcreate ('mpi')

        db.dadd('mpi',('np',4) )
        db.dadd('mpi',('nk',4) )

        #grids - points distance used to calculate grids for scf and phonons calculations
        db.dcreate ('grid')

        db.dadd('grid',('coarse_div',(9,9,9)) )
        db.dadd('grid',('coarse_off',(0,0,0)) )
        db.dadd('grid',('qpoints_div',(3,3,3)) )
        db.dadd('grid',('qpoints_off',(0,0,0)) )
        db.dadd('grid',('dense_div',(18,18,18)) )
        db.dadd('grid',('dense_off',(0,0,0)) )
        db.dadd('grid',('kdistance', 0.3) ) #distance in angstrons of two consecutive grid k-points
        db.dadd('grid',('qdistance', 0.9) ) #distance in angstrons of two consecutive grid q-points
        db.dadd('grid',('calc', 1) ) #habilitate the calculation from reciprocal vector

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

        return db

class EosDatabase (IMakeDatabase):

    def _make_db(self, name):

        _DATABASE = name

        db = pk.load(_DATABASE, False)

        #optional: in case quantum-espresso programs are not in linux $PATH
        #qe - quantum espresso programs' path
        db.dcreate ('calculation')
        db.dadd('calculation',('routine','EQUATION OF STATE') )
        
        #pseudo - pseudo potentials' folder path
        db.dcreate('pseudo')
        db.dadd('pseudo',('pseudo_folder', pseudo.DIR+'/USPP') )

        #mpi - parallen running description
        db.dcreate ('mpi')

        db.dadd('mpi',('np',4) )
        db.dadd('mpi',('nk',4) )

        #grids - points distance used to calculate grids for scf and phonons calculations
        db.dcreate ('grid')

        db.dadd('grid',('kpoints_div',(9,9,9)) )
        db.dadd('grid',('kpoints_off',(0,0,0)) )

        #compression - parameters for cell compression and calculation of energy
        #only three parameters are needed, the fourth must be zero. 
        # If all of them are given, the program ignores parameter "step"

        db.dcreate ('compression')

        db.dadd('compression',('n', 20) )
        db.dadd('compression',('init', 0.8) )
        db.dadd('compression',('final', 1.2) )

        #Quantum-espresso programs' parameters:
        db.dcreate ('pw_par')

        db.dadd('pw_par',('restart_mode','from_scratch') )
        db.dadd('pw_par',('occupations','smearing') )
        db.dadd('pw_par',('smearing','marzari-vanderbilt') )
        db.dadd('pw_par',('degauss',0.05) )
        db.dadd('pw_par',('conv_thr',1e-10) )

        return db

class PhononDatabase (IMakeDatabase):

    def _make_db(self, name):

        _DATABASE = name

        db = pk.load(_DATABASE, False)

        #optional: in case quantum-espresso programs are not in linux $PATH
        #qe - quantum espresso programs' path
        db.dcreate ('calculation')
        db.dadd('calculation',('routine','PH') )

        #pseudo - pseudo potentials' folder path
        db.dcreate('pseudo')
        db.dadd('pseudo',('pseudo_folder',pseudo.DIR+'/USPP') )

        #mpi - parallen running description
        db.dcreate ('mpi')

        db.dadd('mpi',('np',4) )
        db.dadd('mpi',('nk',4) )

        #grids - points distance used to calculate grids for scf and phonons calculations
        db.dcreate ('grid')

        db.dadd('grid',('kpoints_div',(9,9,9)) )
        db.dadd('grid',('kpoints_off',(0,0,0)) )

        #Quantum-espresso programs' parameters:
        db.dcreate ('pw_par')

        db.dadd('pw_par',('restart_mode','from_scratch') )
        db.dadd('pw_par',('occupations','smearing') )
        db.dadd('pw_par',('smearing','marzari-vanderbilt') )
        db.dadd('pw_par',('degauss',0.05) )
        db.dadd('pw_par',('conv_thr',1e-10) )

        db.dcreate ('ph_par')

        db.dadd('ph_par',('tr2_ph',1e-12) )
        db.dadd('ph_par',('recover','.false.') )
        db.dadd('ph_par',('ldisp','.false.') )

        db.dcreate ('dynmat_par')

        db.dadd('dynmat_par',('asr','crystal') )

        return db

        ################################################################################
##----------------------------------------------------------------------------##
################################################################################

if __name__ == '__main__':

    class Test():
        
        db_method = IMakeDatabase()

        def make_database(self, name):
            self.db_method.make_database(name = name)

    #Test TcDatabase()
    E = Test()
    E.db_method = TcDatabase()
    E.make_database(name= 'tc.db')

    #Test PhononDatabase()
    F = Test()
    F.db_method = PhononDatabase()
    F.make_database(name= 'ph.db')

    #Test EosDatabase()
    G = Test()
    G.db_method = EosDatabase()
    G.make_database(name= 'eos.db')