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