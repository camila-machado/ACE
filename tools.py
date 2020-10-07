import ase
import ase.io
import os

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
            value = str(parameters.get(key, False))
            if value:
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
            pseudopotentioals= {'Nb':'Nb.pbe-spn-rrkjus_psl.1.0.0.UPF'}, 
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
