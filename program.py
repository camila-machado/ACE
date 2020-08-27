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

This program provides an interface to calculate critical superconductivity 
temperature using espressotc module from shell and indicates the total 
calculation time.
"""

import time
from espressotc import *

def get_input_data():
    """ 

    This function read input for CriticalTemp instantiation from user shell
    @return: control (str): command for critical temperature calculation
             infile (str) : path of input file
    """
    user_input = sys.argv
    if len(user_input) > 2:
        control = user_input[2]
    else:
        control = 'n'

    infile = user_input[1]  

    return infile, control

################################################################################
##----------------------------------------------------------------------------##
################################################################################

start_time = time.time()

infile, control = get_input_data()
TC = CriticalTemp(infile= infile , control= control)
TC.calculate()

print('\nTotal Execution time:', time.time()-start_time)
"""
if __name__ == '__main__':

    #-------test CriticalTemp--------------
#    H3S_cif = '/home/camila/Documentos/EMA/Program-TC/cif_database/H3S.cif'
#    Nb_cif = '/home/camila/Documentos/EMA/Program-TC/cif_database/Nb.cif'
#    Hg_cif = '/home/camila/Documentos/EMA/Program-TC/cif_database/Hg.cif'

    for element in ['Cd','H3S','Hg','In','Nb','Ru','Sn','Ti','Zn']:
        start_time = time.time()
        print('\n----->CALC TC OF ELEMEMT: ', element)
        infile= '/home/camila/Documentos/EMA/Program-TC/{}/{}_config.db'.format(element, element)
        print(infile)
        TC = CriticalTemp(infile= infile , control= 'n')
        TC.calculate()
        print('\nTotal Execution time:', time.time()-start_time)
"""