#!/usr/bin/env python3
# File classtools.py
# "Assorted class utilities and tools"

class AttrDisplay:    
    """

    Provides an inheritable print overload method that displays    
    instances with their class names and a name=value pair for    
    each attribute stored on the instance itself (but not attrs    
    inherited from its classes). Can be mixed into any class,    
    and will work on any instance.    
    """        
    def gatherAttrs(self):        
        attrs = []        
        for key in sorted(self.__dict__):            
            attrs.append('%s = %s' % (key, getattr(self, key)))        
        return ', '.join(attrs)
    def __str__(self):        
        return '[%s: %s]' % (self.__class__.__name__, self.gatherAttrs())


class WriteInput:
    """

    Provides an inheritable write input method that creates an 
    input file compatible with Qauntum Espresso Programs and
    qeprograms module    
    """   

    def _writeInput(self, section, f_input, file_order, parameters):           
        first_line = '&'+section.upper()+'\n/\n'
        with open(f_input, 'w') as f:     
            f.write(first_line)

        file_order.reverse()
        for key in file_order:
            value = str(parameters.get(key))
            self._addParameter(key= key, value= value, section= section, f_input = f_input) 

    def _addParameter(self, key, value, section, f_input):

        with open(f_input, 'r') as f:
            content = f.readlines()

        value = self._treatValue(word = value)
        section_pos = self._getPosition(content= content, word= section)
        identation = self._setIdentation(word= key)

        new_pos = section_pos + 1
        new_info = '   ' + key + identation + '= ' + value + '\n'
        content.insert(new_pos,new_info)

        with open(f_input, 'w') as f:     
            for line in content:
                f.write(line)
    
    def _isnumber (self, number):
        number = number.replace('e','')
        number = number.replace('-','')
        number = number.replace('.','')
        test = number.isnumeric()
        return test
    
    def _setIdentation(self, word):
        identation = ''
        i = 0
        while i < 17 - len(word):
            identation += ' '
            i += 1
        return identation
    
    def _treatValue(self, word):
        if word.find('true') ==-1 and word.find('false') == -1:
            if self._isnumber(word) == False:
                word = "'"+word+"'"
        return word
    
    def _getPosition(self, content, word):
        for line in content:
            if line.find(word.upper() or word.lower()) > -1:
                position = content.index(line)
                break
        return position


if __name__ == '__main__':    
    
    class TopTest(AttrDisplay):        
        count = 0        
        def __init__(self):            
            self.attr1 = TopTest.count            
            self.attr2 = TopTest.count+1            
            TopTest.count += 2

    class SubTest(TopTest):        
            pass    

    X, Y = TopTest(), SubTest()    
    print(X)                         # Show all instance attrs    
    print(Y)                         # Show lowest class name

    class WriteTest(WriteInput):

        def __init__(self, section = 'input', f_input = 'input.in',
                    file_order = ['var1', 'var2', 'var3'],
                    parameters = {'var1':1, 'var2':2, 'var3':3}):
            self.section = section 
            self.f_input = f_input 
            self.file_order = file_order 
            self.parameters = parameters

        def write(self):
            self._writeInput(section = self.section, 
                            f_input = self.f_input, 
                            file_order = self.file_order, 
                            parameters = self.parameters)


    Z = WriteTest()
    Z.write()