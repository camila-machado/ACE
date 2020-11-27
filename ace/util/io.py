# Copyright 2020, Camila Machado de Araújo
# (see accompanying license files for details).

import os

import ase
import ase.io
import pickledb as pk

"""  
"""   
def write_sections(sections, filename, file_order, parameters):

    assert type(sections) == list, 'Sections input must be a list'

    #write sections names
    section_lines = []
    for s in sections:
        section_lines.append('&'+s.upper()+'\n/\n')

    with open(filename, 'w') as f:
        for line in section_lines:     
            f.write(line)

    #add respective parameteres above each section name
    for i in range(0, len(sections)):
        _add_parameters(parameters = parameters[i],  section= sections[i], 
                        filename = filename, file_order = file_order[i])

def write_section(section, filename, file_order, parameters):           
    
    section_line = '&'+section.upper()+'\n/\n'
    with open(filename, 'w') as f:     
        f.write(section_line)

    _add_parameters(parameters = parameters,  section= section, 
                    filename = filename, file_order = file_order)

def _add_parameters(parameters, section, filename, file_order):
    file_order.reverse()
    for key in file_order:
        value = parameters.get(key, False)
        if value:
            value = str(value)
            add_parameter(key = key, value= value, section= section, 
                          filename= filename)

def add_parameter(key, value, section, filename):
    line = _make_line(key= key, value = value)
    add_line(line= line, position= section, filename= filename)

def add_line(line, position, filename):
    with open(filename, 'r') as f:
        content = f.readlines()

    line_pos = _get_position(content= content, word= position) + 1 

    content.insert(line_pos,line)

    with open(filename, 'w') as f:     
        for line in content:
            f.write(line)

def _get_position(content, word):
    for line in content:
        if line.find(word.upper() or word.lower()) > -1:
            position = content.index(line)
            break
    return position

def _make_line(key, value):
    value = _treat_value(word = value)
    identation = _set_identation(word= key)

    line = '   ' + key + identation + '= ' + value + '\n'
    
    return line

def _set_identation(word):
    identation = ''
    i = 0
    while i < 17 - len(word):
        identation += ' '
        i += 1
    return identation

def _treat_value(word):
    
    if not(_isbool(word)):
        if  not(_isnumber(word)):
            word = "'"+word+"'"
    return word

def _isbool(word):
    return (word.find('true') !=-1) or (word.find('false') != -1)

def _isnumber (number):
    number = number.replace('e','')
    number = number.replace('-','')
    number = number.replace('.','')
    number = number.replace(',','')
    number = number.replace(' ','')
    number = number.strip('()')
    test = number.isnumeric()
    return test


def write_lambda(filename, dyndir, sigma_omega, mu, fildyn):

    qpts_info, qpts_num, dynfiles = _get_dyninfo(dyndir= dyndir, fildyn= fildyn)
    multiplicity = _get_multplicity(dynfiles = dynfiles)
    upperfreq = _get_upperfreq(dynfiles = dynfiles)

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

def _get_dyninfo(dyndir, fildyn):

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

def _get_multplicity(dynfiles):
    multiplicity_values = []
    for element in dynfiles:
        multiplicity = 0
        for line in element:
            if line.find('Dynamical  Matrix in cartesian axes') > -1:
                multiplicity += 1
        multiplicity_values.append(multiplicity)

    return multiplicity_values

def _get_upperfreq(dynfiles):
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

class ReadCif(IReadStrategy):

    def read(self, filepath):
        atoms = ase.io.read(filename = filepath,
                            format= 'cif',
                            subtrans_included = False,
                            primitive_cell= True)
        return atoms

class ReadDatabase(IReadStrategy):

    def read(self, filepath):
        db = pk.load(filepath, False)

        sections = list(db.getall())

        parameters = []
        keys = []
        for section in sections:
            parameters.append(dict(db.dgetall(section))) 
            keys.append(list(db.dkeys(section)))

        return sections, parameters, keys

class ReadInputvar(IReadStrategy):

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
                key, value = self._parse_line(line)
                if key:
                    sec_par.update({key:value})
                    sec_key.append(key)

        parameters.append(sec_par)
        keys.append(sec_key)

        return sections, parameters, keys

    def _parse_line(self, line):
        eq = line.find('=')
        if eq == -1:
            return 0, 0
        key = line[:eq].strip()
        value = line[eq+1:].strip()

        return key, self._parse_value(value)

    def _parse_value(self, expr):
        try:
            return eval(expr)
        except:
            return expr
        else:
            return expr


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
        
        read_method = IReadStrategy()

        def read(self):
            assert self.filepath, 'filepath must be defined'
            content = self.read_method.read(filepath = self.filepath)
            return content
    
    #Test ReadCif()
    H = Test2()
    H.read_method = ReadCif()
    H.filepath = 'Nb.cif'
    print('read:', H.filepath, '\n\n', H.read(), '\n')

    #Test ReadDatabase()
    I = Test2()
    I.read_method = ReadDatabase()
    I.filepath = 'tc.db'
    sec, par, key = I.read()
    print('read:', I.filepath)
    print('\n', 'sections:', sec, '\n\n','parameters:', par, '\n\n', 'keys:', key, '\n')

    #Test ReadInputvar()
    J = Test2()
    J.read_method = ReadInputvar()
    J.filepath = 'test3.in'
    sec, par, key = J.read()
    print('read:', J.filepath)
    print('\n', 'sections:', sec, '\n\n','parameters:', par, '\n\n', 'keys:', key), '\n'