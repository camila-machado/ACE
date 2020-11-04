import os

import ase
import ase.io
import pickledb as pk

class IWriteStrategy:

    def write(**keyargs):
        pass

class WriteSection(IWriteStrategy):
    """

    Provides an inheritable write input method that creates an 
    input file compatible with Qauntum Espresso Programs and
    qeprograms module    
    """   

    def write(self, section, filename, file_order, parameters):           
        first_line = '&'+section.upper()+'\n/\n'
        with open(filename, 'w') as f:     
            f.write(first_line)

        self._addParameters(parameters = parameters,  section= section, 
                              filename = filename, file_order = file_order)

    def _addParameters(self, parameters, section, filename, file_order):
        file_order.reverse()
        for key in file_order:
            value = parameters.get(key, False)
            if value:
                value = str(value)
                self.addParameter(key = key, value= value, section= section, 
                                  filename= filename)

    def addParameter(self, key, value, section, filename):
        line = self._mke_key_value_line(key= key, value = value)
        self.addLine(line= line, position= section, filename= filename)

    def addLine(self, line, position, filename):
        with open(filename, 'r') as f:
            content = f.readlines()

        line_pos = self._getPosition(content= content, word= position) + 1 

        content.insert(line_pos,line)

        with open(filename, 'w') as f:     
            for line in content:
                f.write(line)
    
    def _getPosition(self, content, word):
        for line in content:
            if line.find(word.upper() or word.lower()) > -1:
                position = content.index(line)
                break
        return position

    def _mke_key_value_line(self, key, value):
        value = self._treatValue(word = value)
        identation = self._setIdentation(word= key)

        line = '   ' + key + identation + '= ' + value + '\n'
        
        return line

    def _setIdentation(self, word):
        identation = ''
        i = 0
        while i < 17 - len(word):
            identation += ' '
            i += 1
        return identation
    
    def _treatValue(self, word):
        
        if not(self._isbool(word)):
            if  not(self._isnumber(word)):
                word = "'"+word+"'"
        return word
    
    def _isbool(self, word):
        return (word.find('true') !=-1) or (word.find('false') != -1)

    def _isnumber (self, number):
        number = number.replace('e','')
        number = number.replace('-','')
        number = number.replace('.','')
        number = number.replace(',','')
        number = number.replace(' ','')
        number = number.strip('()')
        test = number.isnumeric()
        return test

class WriteSections(WriteSection):

    def write(self, section, filename, file_order, parameters):

        assert type(section) == list, 'Section input must be a list'

        section_lines = []
        for sec in section:
            section_lines.append('&'+sec.upper()+'\n/\n')

        with open(filename, 'w') as f:
            for line in section_lines:     
                f.write(line)

        for i in range(0, len(section)):
            self._addParameters(parameters = parameters[i],  section= section[i], 
                                  filename = filename, file_order = file_order[i])

class WriteScf(WriteSection):

    def write(self, filename, images, input_data, pseudopotentials, kpts, koffset):
        ase.io.write(filename = filename,
                     images = images,
                     format = 'espresso-in',
                     input_data = input_data,
                     pseudopotentials = pseudopotentials,
                     kpts = kpts,
                     koffset = koffset)

class WriteLambda(IWriteStrategy):

    def write(self, filename, dyndir, sigma_omega, mu, fildyn):

        qpts_info, qpts_num, dynfiles = self._getDynInfo(dyndir= dyndir, fildyn= fildyn)
        multiplicity = self._getMultplicity(dynfiles = dynfiles)
        upperfreq = self._getUpperFreq(dynfiles = dynfiles)

        part_1 = str(upperfreq)+'  '+str(sigma_omega)+'  1\n'

        part_2 = qpts_info
        for i in range(1,qpts_num+1):
            sufix = '  '+str(multiplicity[i-1]) +'\n'
            part_2[i] = part_2[i].replace('\n', sufix)

        part_3 = []
        aux = 'elph_dir/elph.inp_lambda.'
        for num in range(1,qpts_num+1):
            line = aux + str(num) + '\n'
            part_3.append(line)
        part_3.append(str(mu)+'\n')

        lambda_info =[]
        lambda_info.append(part_1)
        lambda_info.extend(part_2)
        lambda_info.extend(part_3)

        with open(filename, 'w') as f:     
            for line in lambda_info:
                f.write(line)

    def _getDynInfo(self, dyndir, fildyn):

        current_dir = os.getcwd()
        os.chdir(os.path.join(current_dir, dyndir))

        #read .dyn0 file
        filename = fildyn + '0'
        with open(filename, 'r') as f:
            content = f.readlines()
            content.pop(0)
        qpts_info = content 
        qpts_num = int(qpts_info[0])
        
        #read .dyn other files
        dyn_info = []
        for num in range(1, qpts_num+1):
            filename = fildyn + str(num)
            with open(filename) as f:
                content = f.readlines()
            dyn_info.append(content)
        
        os.chdir(current_dir)

        return qpts_info, qpts_num, dyn_info

    def _getMultplicity(self, dynfiles):
        multiplicity_values = []
        for element in dynfiles:
            multiplicity = 0
            for line in element:
                if line.find('Dynamical  Matrix in cartesian axes') > -1:
                    multiplicity += 1
            multiplicity_values.append(multiplicity)

        return multiplicity_values

    def _getUpperFreq(self, dynfiles):
        freqpts_info = []
        for element in dynfiles:
            for line in element:
                if line.find('freq') > -1:
                    freq = line.split('=')[1]
                    freq = freq.replace('[THz]','').strip()
                    freqpts_info.append(float(freq))
            
        max_freq = max(freqpts_info)
        upper_freq = int(max_freq+1)

        return upper_freq

################################################################################
##----------------------------------------------------------------------------##
################################################################################

class IReadStrategy:

    def read(**keyargs):
        pass

class ReadCif:

    def read(self, filepath):
        atoms_structure = ase.io.read(filename = filepath,
                                    format= 'cif',
                                    subtrans_included = False,
                                    primitive_cell= True)
        return atoms_structure

class ReadDatabase:

    def read(self, filepath):
        db = pk.load(filepath, False)

        sections = list(db.getall())

        parameters = []
        keys = []
        for section in sections:
            parameters.append(dict(db.dgetall(section))) 
            keys.append(list(db.dkeys(section)))

        return sections, parameters, keys

class ReadInputvar:

    def read(self, filepath):

        with open(filepath, 'r') as f:
            file_content = f.readlines()

        sections = []
        parameters = []
        keys = []

        sec_key = []
        sec_par = {}
        for line in file_content:

            if line.count('&'):
                section = line.strip('&\n ').lower()
                sections.append(section)
                if sec_par:
                    parameters.append(sec_par.copy())
                    sec_par.clear()               
            else:
                key, value = self._parseLine(line)
                if key:
                    sec_par.update({key:value})
                    sec_key.append(key)

        parameters.append(sec_par)
        keys.append(sec_key)

        return sections, parameters, keys

    def _parseLine(self, line):
        eq = line.find('=')
        if eq == -1:
            return 0, 0
        key = line[:eq].strip()
        value = line[eq+1:].strip()

        return key, self._parseValue(value)

    def _parseValue(self, expr):
        try:
            return eval(expr)
        except:
            return expr
        else:
            return expr


################################################################################
##----------------------------------------------------------------------------##
################################################################################

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
        
        #qe - quantum espresso programs' path
        db.dcreate ('qe')
        db.dadd('qe',('qe_programs','') )

        #pseudo - pseudo potentials' folder path
        db.dcreate('pseudo')
        db.dadd('pseudo',('pseudo_folder','../pseudo/USPP') )

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
        
        #qe - quantum espresso programs' path
        db.dcreate ('qe')
        db.dadd('qe',('qe_programs','') )

        #pseudo - pseudo potentials' folder path
        db.dcreate('pseudo')
        db.dadd('pseudo',('pseudo_folder','../pseudo/USPP') )

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
        
        #qe - quantum espresso programs' path
        db.dcreate ('qe')
        db.dadd('qe',('qe_programs','') )

        #pseudo - pseudo potentials' folder path
        db.dcreate('pseudo')
        db.dadd('pseudo',('pseudo_folder','../pseudo/USPP') )

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

        return db

################################################################################
##----------------------------------------------------------------------------##
################################################################################

if __name__ == '__main__':

    class Test:

        write_method = IWriteStrategy()

        def write(self, **keyargs):
            self.write_method.write(**keyargs)

    #Test WriteSection()
    A = Test()
    A.write_method = WriteSection()
    A.write(section = 'input', 
            filename = 'test1.in',
            file_order = ['var1', 'var2', 'var3'],
            parameters = {'var1':1, 'var2':2, 'var3':3})

    #Test WriteSections()
    B = Test()
    B.write_method = WriteSections()
    B.write(section = ['input_num', 'input_abc'], 
            filename = 'test2.in',
            file_order = [['var1', 'var2', 'var3'],['varA', 'varB', 'varC']],
            parameters = [{'var1':1, 'var2':2, 'var3':3}, 
                          {'varA':'a', 'varB':'b', 'varC':'c'}])

    #Test WriteScf()
    image = ase.io.read(filename = 'Nb.cif',
                         format='cif',
                         subtrans_included = False,
                         primitive_cell= True)

    C = Test()
    C.write_method = WriteScf()
    C.write(filename= 'test3.in', 
            images= image, 
            input_data= {'prefix':'Nb','pseudo_dir':'/','outdir':'/out',
                         'occupations':'smearing', 'ecutwfc':50,'ecutrho':50}, 
            pseudopotentials= {'Nb':'Nb.pbe-spn-rrkjus_psl.1.0.0.UPF'}, 
            kpts= (3,3,3),
            koffset= (0,0,0))

    #Test WriteLambda()
    D = Test()
    D.write_method = WriteLambda()
    D.write(filename= 'test4.in', 
            dyndir= 'dyn/',
            sigma_omega= 0.12, 
            mu= 0.16, 
            fildyn= 'Nb.dyn')


    class Test2():
        
        db_method = IMakeDatabase()

        def make_database(self, name):
            self.db_method.make_database(name = name)

    #Test TcDatabase()
    E = Test2()
    E.db_method = TcDatabase()
    E.make_database(name= 'tc.db')

    #Test PhononDatabase()
    F = Test2()
    F.db_method = PhononDatabase()
    F.make_database(name= 'ph.db')

    #Test EosDatabase()
    G = Test2()
    G.db_method = EosDatabase()
    G.make_database(name= 'eos.db')

    class Test3():
        
        read_method = IReadStrategy()

        def read(self):
            assert self.filepath, 'filepath must be defined'
            content = self.read_method.read(filepath = self.filepath)
            return content
    
    #Test ReadCif()
    H = Test3()
    H.read_method = ReadCif()
    H.filepath = 'Nb.cif'
    print('read:', H.filepath, '\n\n', H.read(), '\n')

    #Test ReadDatabase()
    I = Test3()
    I.read_method = ReadDatabase()
    I.filepath = 'tc.db'
    sec, par, key = I.read()
    print('read:', I.filepath)
    print('\n', 'sections:', sec, '\n\n','parameters:', par, '\n\n', 'keys:', key, '\n')

    #Test ReadInputvar()
    J = Test3()
    J.read_method = ReadInputvar()
    J.filepath = 'test3.in'
    sec, par, key = J.read()
    print('read:', J.filepath)
    print('\n', 'sections:', sec, '\n\n','parameters:', par, '\n\n', 'keys:', key), '\n'