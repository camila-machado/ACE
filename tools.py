class WriteStrategy:
    """

    Provides an inheritable write input method that creates an 
    input file compatible with Qauntum Espresso Programs and
    qeprograms module    
    """   

    def section_single(self, section, f_input, file_order, parameters):           
        first_line = '&'+section.upper()+'\n/\n'
        with open(f_input, 'w') as f:     
            f.write(first_line)

        self._addParameters(parameters = parameters,  section= section, 
                              f_input = f_input, file_order = file_order)

    def sections(self, section, f_input, file_order, parameters):

        section_lines = []
        for sec in section:
            section_lines.append('&'+sec.upper()+'\n/\n')

        with open(f_input, 'w') as f:
            for line in section_lines:     
                f.write(line)

        for i in range(0, len(section)):
            self._addParameters(parameters = parameters[i],  section= section[i], 
                                  f_input = f_input, file_order = file_order[i])

    def _addParameters(self, parameters, section, f_input, file_order):

        file_order.reverse()
        for key in file_order:
            value = str(parameters.get(key))
            self.addParameter(key = key, value= value, section= section, f_input= f_input)

    def addParameter(self, key, value, section, f_input):
        line = self._mke_key_value_line(key= key, value = value)
        self.addLine(line= line, position= section, f_input= f_input)

    def addLine(self, line, position, f_input):
        with open(f_input, 'r') as f:
            content = f.readlines()

        line_pos = self._getPosition(content= content, word= position) + 1 

        content.insert(line_pos,line)

        with open(f_input, 'w') as f:     
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
    

################################################################################
##----------------------------------------------------------------------------##
################################################################################

if __name__ == '__main__':
    
    class WriteTest(WriteStrategy):

        def __init__(self, section = 'input', f_input = 'input.in',
                    file_order = ['var1', 'var2', 'var3'],
                    parameters = {'var1':1, 'var2':2, 'var3':3}):
            self.section = section 
            self.f_input = f_input 
            self.file_order = file_order 
            self.parameters = parameters

        def write(self):
            self.section_single(section = self.section, 
                                f_input = self.f_input, 
                                file_order = self.file_order, 
                                parameters = self.parameters)

    Z = WriteTest()
    Z.write()

    class WriteTest2(WriteStrategy):

        def __init__(self, section = ['input_num', 'input_abc'], f_input = 'input2.in',
                    file_order = [['var1', 'var2', 'var3'],['varA', 'varB', 'varC']],
                    parameters = [{'var1':1, 'var2':2, 'var3':3}, 
                                  {'varA':'a', 'varB':'b', 'varC':'c'}]):
            self.section = section
            self.f_input = f_input 
            self.file_order = file_order 
            self.parameters = parameters

        def write(self):
            self.sections(section = self.section, 
                          f_input = self.f_input, 
                          file_order = self.file_order, 
                          parameters = self.parameters)


    T = WriteTest2()
    T.write()